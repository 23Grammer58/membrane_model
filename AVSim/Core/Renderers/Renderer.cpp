//
// Created by alex on 30.06.2020.
//

#include "Renderer.h"

void World3d::DefaultRenderer::delete_renderer(World3d::ObjectID id, Object3D &obj) {
    obj._renderer = nullptr;
    _renderers.erase(id);
    _rendered_shapes[id].second = InterconnectDataState::Deleted;
}

void World3d::DefaultRenderer::render() {
    g_interconnection.set_rendered_object_data(PartInterconnect::BackEnd, _rendered_shapes);
}

void World3d::DefaultRenderer::set_renderer(World3d::ObjectID id, Object3D &obj, bool modifiable) {
    auto it = _rendered_shapes.insert({id, pair<RenderedObject, InterconnectDataState>()});
    it.first->second.second = InterconnectDataState::New;
    auto& ro = it.first->second.first;
    ro.name = obj.name;
    ro.value.second = true;
    ro.vertex.second = true;
    ro.face.second = true;
    ro.vertex_lbl.second = true;
    ro.face.first.reserve(obj.m_mesh.num_faces()*3);
    ro.value.first.reserve(obj.m_mesh.num_vertices());
    ro.vertex.first.reserve(obj.m_mesh.num_vertices());
    for (auto f: obj.m_mesh.faces()){
        auto v = vert_around(obj.m_mesh, f);
        for (auto i: v)
            ro.face.first.push_back(i);
    }
    for (auto v: obj.m_mesh.vertices()){
        int lbl = obj.m_boundary[v];
        if (lbl > 0)
            ro.vertex_lbl.first.push_back({v, lbl});
    }
    for (auto v: obj.m_mesh.vertices()){
        ro.vertex.first.push_back(array<float, 3>{
                static_cast<float>(obj.m_x[v].x()),
                static_cast<float>(obj.m_x[v].y()),
                static_cast<float>(obj.m_x[v].z())});
        ro.value.first.push_back(static_cast<float>(sqrt(obj.m_F[v].squared_length())));
    }
    auto render = [it, mod = modifiable](Object3D& obj){
        auto& ro = it.first->second.first;
        auto& state = it.first->second.second;
        if (state == InterconnectDataState::Deleted) return -1;
        if (state == InterconnectDataState::NotModified && !mod) return 0;
        if (state == InterconnectDataState::NotModified)
            state = InterconnectDataState::Modified;
        ro.value.second = true;
        ro.vertex.second = true;
        ro.value.first.resize(obj.m_mesh.num_vertices());
        ro.vertex.first.resize(obj.m_mesh.num_vertices());
        for (auto v: obj.m_mesh.vertices()){
            ro.vertex.first[v] = array<float, 3>{
                    static_cast<float>(obj.m_x[v].x()),
                    static_cast<float>(obj.m_x[v].y()),
                    static_cast<float>(obj.m_x[v].z())};
            ro.value.first[v] = static_cast<float>(sqrt(obj.m_F[v].squared_length()));
        }
        return 0;
    };
    auto rit = _renderers.insert({id, move(render)});
    obj.set_renderer(rit.first->second);
}
