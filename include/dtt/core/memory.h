#pragma once

#include "dtt/util.h"

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <type_traits>

namespace dtt::core {
/**
 * @file memory.h
 * @brief Definitions related to memory management utilities.
 * Provides functions for aligned memory allocation and deallocation,
 * for both host and device memory.
 * @author Yannik Rüfenacht
 */
enum class MemoryType { HOST, DEVICE };

template <typename T> class Memory {
  public:
    // RAII ownership: default ctors, disable copy, allow move
    Memory() = default;
    Memory(Memory &&other) noexcept {
        ptr_ = other.ptr_;
        size_ = other.size_;
        type_ = other.type_;
        alignment_ = other.alignment_;
        other.ptr_ = nullptr;
        other.size_ = 0;
    }

    Memory &operator=(Memory &&other) noexcept {
        if (this != &other) {
            release();
            ptr_ = other.ptr_;
            size_ = other.size_;
            type_ = other.type_;
            alignment_ = other.alignment_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    Memory(const Memory &other) = delete;
    Memory &operator=(const Memory &other) = delete;
    ~Memory() { release(); }

    /**
     * @brief Allocates memory with specified memory type and alignment.
     * @details For HOST memory, uses aligned_alloc for SIMD-friendly allocation.
     *          For DEVICE memory, uses cudaMalloc (requires DTT_ENABLE_CUDA).
     * @param count Number of elements to allocate.
     * @param type Memory type (HOST or DEVICE).
     * @param alignment Alignment in bytes (default is 64).
     * @return A Memory object managing the allocated memory.
     */
    static Memory allocate(std::size_t count, MemoryType type = MemoryType::HOST,
                           std::size_t alignment = 64) {
        Memory mem;
        mem.size_ = count;
        mem.type_ = type;
        if (count == 0) {
            mem.ptr_ = nullptr;
            return mem;
        }
        if (type == MemoryType::HOST) {
            if (alignment < alignof(T))
                alignment = alignof(T);
            // make alignment a power of two
            if ((alignment & (alignment - 1)) != 0) {
                // round up to next power of two
                std::size_t a = 1;
                while (a < alignment)
                    a <<= 1;
                alignment = a;
            }
            mem.alignment_ = alignment;
            void *p = nullptr;
            int rc = posix_memalign(&p, alignment, count * sizeof(T));
            if (rc != 0 || !p)
                throw makeError<MemoryError>("posix_memalign failed");
            mem.ptr_ = static_cast<T *>(p);
        }

#ifdef DTT_ENABLE_CUDA
        else if (type == MemoryType::DEVICE) {
            // device memory allocation
            DTT_CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&mem.ptr_), count * sizeof(T)));
        }
#endif
        return mem;
    }

    /**
     * @brief Extends the memory to a new size, preserving existing data.
     * @param new_size The new size in number of elements.
     */
    void extend(std::size_t new_size) {
        if (new_size <= size_)
            return;
        if (type_ == MemoryType::DEVICE) {
            throw makeError<MemoryError>("extend not supported for DEVICE memory");
        }
        Memory new_mem = allocate(new_size, type_);
        if (valid())
            std::memcpy(new_mem.ptr_, ptr_, size_ * sizeof(T));
        release();
        ptr_ = new_mem.ptr_;
        size_ = new_mem.size_;
        type_ = new_mem.type_;
        new_mem.ptr_ = nullptr;
        new_mem.size_ = 0;
        new_mem.alignment_ = 0;
    }

    /**
     * @brief Releases the allocated memory.
     */
    void release() {
        if (ptr_) {
            if (type_ == MemoryType::HOST) {
                std::free(ptr_);
            }
#ifdef DTT_ENABLE_CUDA
            else if (type_ == MemoryType::DEVICE) {
                DTT_CUDA_CHECK(cudaFree(ptr_));
            }
#endif
            ptr_ = nullptr;
            size_ = 0;
        }
    }

    /// Accessors
    T *data() { return ptr_; }
    const T *data() const { return ptr_; }
    std::size_t size() const { return size_; }
    bool valid() const { return size_ == 0 || ptr_ != nullptr; }

    T &operator[](std::size_t i) { return ptr_[i]; }
    const T &operator[](std::size_t i) const { return ptr_[i]; }

  private:
    T *ptr_{nullptr};
    std::size_t size_{0};
    MemoryType type_{MemoryType::HOST};
    std::size_t alignment_{alignof(T)};
};
} // namespace dtt::core