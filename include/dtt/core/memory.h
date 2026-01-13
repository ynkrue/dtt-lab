#pragma once

#include "dtt/util.h"

#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <type_traits>
#include <new>

namespace dtt::core
{
    /**
    * @file memory.h
    * @brief Definitions related to memory management utilities.
    * Provides functions for aligned memory allocation and deallocation,
    * for both host and device memory.
    * @author Yannik Rüfenacht
    */
    enum class MemoryType {
        HOST,
        DEVICE
    };

    template <typename T>
    class Memory {
    public:
        // RAII ownership: default ctors, disable copy, allow move
        Memory() = default;
        Memory(Memory&& other) noexcept = default;
        Memory& operator=(Memory&& other) noexcept = default;
        Memory(const Memory& other) = delete;
        Memory& operator=(const Memory& other) = delete;
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
        static Memory allocate(std::size_t count, MemoryType type = MemoryType::HOST, std::size_t alignment = 64) {
            Memory mem;
            mem.size_ = count;
            mem.type_ = type;
            if (type == MemoryType::HOST) {
#pragma GCC diagnostic push // buf in GCC 10+ complains about maybe uninitialized ptr_
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
                // aligned memory allocation (simd friendly)
                const auto bytes = static_cast<std::size_t>(count) * sizeof(T);
                mem.ptr_ = static_cast<T*>(std::aligned_alloc(alignment, bytes));
                if (!mem.ptr_)
                    throw makeError<MemoryError>("aligned malloc failed");
#pragma GCC diagnostic pop
            }
#ifdef DTT_ENABLE_CUDA
            else if (type == MemoryType::DEVICE) {
                // device memory allocation
                DTT_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&mem.ptr_), count * sizeof(T)));
            }
#endif
            return mem;
        }

        /**
         * @brief Extends the memory to a new size, preserving existing data.
         * @param mem The Memory object to extend.
         * @param new_size The new size in number of elements.
         * @return Reference to the extended Memory object.
         */
        static Memory& extend(Memory& mem, std::size_t new_size) {
            if (new_size <= mem.size_ || mem.type_ == MemoryType::DEVICE)
                return mem;
            Memory new_mem = allocate(new_size, mem.type_);
            if (mem.valid()) {
                std::memcpy(new_mem.data(), mem.data(), mem.size_ * sizeof(T));
            }
            mem.release();
            mem = std::move(new_mem);
            return mem;
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
        T* data() { return ptr_; }
        const T* data() const { return ptr_; }
        std::size_t size() const { return size_; }
        bool valid() const { return size_ == 0 || ptr_ != nullptr; }
    private:
        T* ptr_{nullptr};
        std::size_t size_{0};
        MemoryType type_{MemoryType::HOST};
    };
} // namespace dtt::core