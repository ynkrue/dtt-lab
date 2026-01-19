#include "dtt/core/forces.h"

namespace dtt::core {

ForcesBuffer ForcesBuffer::make(std::size_t count, dtt::core::MemoryType type,
                                std::size_t alignment) {
    ForcesBuffer buffer;
    buffer.count = count;

    buffer.fx = dtt::core::Memory<double>::allocate(count, type, alignment);
    buffer.fy = dtt::core::Memory<double>::allocate(count, type, alignment);

    return buffer;
}

ForcesView ForcesBuffer::view() {
    ForcesView v;
    v.fx = fx.data();
    v.fy = fy.data();
    v.count = count;
    return v;
}

ConstForcesView ForcesBuffer::const_view() const {
    ConstForcesView v;
    v.fx = fx.data();
    v.fy = fy.data();
    v.count = count;
    return v;
}

bool ConstForcesView::valid() const {
    return fx != nullptr && fy != nullptr && count > 0;
}

bool ForcesView::valid() const {
    return fx != nullptr && fy != nullptr && count > 0;
}

} // namespace dtt::core