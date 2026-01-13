#pragma once

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
    }

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
        void release();

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