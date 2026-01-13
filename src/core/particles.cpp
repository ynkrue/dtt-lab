#include "dtt/core/particles.h"

namespace dtt:core
{
    ParticlesBuffer ParticlesBuffer::make(std::size_t count, dtt::core::MemoryType type, std::size_t alignment) {
        ParticlesBuffer buffer;
        buffer.count = count;

        buffer.x = dtt::core::Memory<double>::allocate(count, type, alignment);
        buffer.y = dtt::core::Memory<double>::allocate(count, type, alignment);
        buffer.vx = dtt::core::Memory<double>::allocate(count, type, alignment);
        buffer.vy = dtt::core::Memory<double>::allocate(count, type, alignment);
        buffer.mass = dtt::core::Memory<double>::allocate(count, type, alignment);

        return buffer;
    }

    ParticlesView ParticlesBuffer::view() {
        ParticlesView v;
        v.x = x.data();
        v.y = y.data();
        v.vx = vx.data();
        v.vy = vy.data();
        v.mass = mass.data();
        v.count = count;
        return v;
    }

    ConstParticlesView ParticlesBuffer::const_view() const {
        ConstParticlesView v;
        v.x = x.data();
        v.y = y.data();
        v.vx = vx.data();
        v.vy = vy.data();
        v.mass = mass.data();
        v.count = count;
        return v;
    }
}