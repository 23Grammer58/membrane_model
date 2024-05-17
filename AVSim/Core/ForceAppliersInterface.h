//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_FORCEAPPLIERSINTERFACE_H
#define AVSYM_FORCEAPPLIERSINTERFACE_H
#include "Object3D.h"

namespace World3d {
    typedef std::function<int(Object3D&)> ForceApplier;
};

#endif //AVSYM_FORCEAPPLIERSINTERFACE_H
