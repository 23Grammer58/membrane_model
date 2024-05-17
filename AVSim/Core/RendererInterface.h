//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_RENDERERINTERFACE_H
#define AVSYM_RENDERERINTERFACE_H

#include "Object3D.h"
namespace World3d {
    class RendererBase {
    public:
        virtual void set_renderer(ObjectID id, Object3D &obj, bool modifiable = true) { (void) id, (void) obj, (void) modifiable; };
        virtual void delete_renderer(ObjectID id, Object3D &obj) { (void) id, (void) obj; };
        virtual void render() {};
    };
};

#endif //AVSYM_RENDERERINTERFACE_H
