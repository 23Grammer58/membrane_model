//
// Created by alex on 30.06.2020.
//

#ifndef AORTIC_VALVE_RENDERER_H
#define AORTIC_VALVE_RENDERER_H

#include "GuiBackInterconnection.h"
#include "../World.h"

namespace World3d {
    class DefaultRenderer: public RendererBase {
        map<ObjectID, Processor> _renderers;
        map<int, pair<RenderedObject, InterconnectDataState>> _rendered_shapes;
        GuiBackInterconnection* interconnection = &g_interconnection;
    public:
        GuiBackInterconnection* get_interconnector() { return interconnection; }
        void set_interconnector(GuiBackInterconnection* back) { interconnection = back; }

        void delete_renderer(ObjectID id, Object3D& obj) override;
        void render() override final;
        void set_renderer(ObjectID id, Object3D& obj, bool modifiable = true) override final;
    };
}


#endif //AORTIC_VALVE_RENDERER_H
