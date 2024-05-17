#include "BulletVolumeTree.h"

namespace World3d{
    std::unique_ptr<VolumeTreeBase> _CollisionThicknessCompressed_internal_get_volume_tree(){ return std::make_unique<BulletVolumeTree>(); }
}