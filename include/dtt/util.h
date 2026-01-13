#pragma once

#include <exception>
#include <string>
#include <type_traits>

#ifdef DTT_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace dtt {

/**
 * @file util.h
 * @brief Lightweight error helpers and CUDA check macro for DTT.
 */
class Error : public std::exception {
public:
    explicit Error(std::string msg = "Unknown error") : message_(std::move(msg)) {}
    const char *what() const noexcept override { return message_.c_str(); }
    virtual ~Error() = default;
    std::string message_;
};

/**
 * @brief Create an error of type T with a given message.
 */
template <typename T, typename... Args>
T makeError(std::string msg, Args &&...args)
    requires std::is_base_of_v<Error, T>
{
    T err(std::forward<Args>(args)...);
    err.message_ = std::move(msg);
    return err;
}

#ifdef DTT_ENABLE_CUDA
/**
 * @brief Wrapper for CUDA errors.
 */
struct CudaError : public Error {
    explicit CudaError(cudaError_t err) : Error(), cu_error(err) {}
    cudaError_t cu_error;
};
#endif

/**
 * @brief Generic memory allocation error.
 */
struct MemoryError : public Error {
    using Error::Error;
};

} // namespace dtt

#ifdef DTT_ENABLE_CUDA
/**
 * @brief Evaluates `call` and throws a CudaError on failure.
 */
#define DTT_CUDA_CHECK(call)                                                              \
    do {                                                                                  \
        cudaError_t _err = (call);                                                        \
        if (_err != cudaSuccess) {                                                        \
            throw dtt::makeError<dtt::CudaError>(                                         \
                "CUDA error at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) \
                    + " - " + std::string(cudaGetErrorString(_err)),                      \
                _err);                                                                    \
        }                                                                                 \
    } while (0)
#else
#define DTT_CUDA_CHECK(call) (call)
#endif

/**
 * @brief Mark unreachable code for the compiler.
 */
#if defined(__GNUC__) || defined(__clang__)
#define DTT_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#define DTT_UNREACHABLE() __assume(false)
#else
#include <utility>
#define DTT_UNREACHABLE() std::unreachable()
#endif
