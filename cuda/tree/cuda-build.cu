#include "dtt/tree/build.h"
#include "dtt/tree/tree.h"

namespace dtt::tree::cu {

void build_tree_cuda(const dtt::core::ConstParticlesView &particles, const BuildParams &params) {
    (void)particles;
    (void)params;
}

} // namespace dtt::tree::cu