//
// Created by alex on 21.06.2020.
//

#include "GuiApplication.h"
#include "../TriangularMeshHelpers.h"

GuiApplication::GuiApplication(const Arguments& arguments):
        Platform::Application{arguments, Configuration{}
                .setTitle("Membranes Viewer")
                .setWindowFlags(Configuration::WindowFlag::Resizable)}{

    GL::Renderer::setClearColor(_controller._clearColor);
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);

    _controller.init();
    _shaders = CustomDrawable::init_shaders();

    _imgui = ImGuiIntegration::Context(windowSize());/*Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize());*/

    /* Set up proper blending to be used by ImGui. There's a great chance
       you'll need this exact behavior for the rest of your scene. If not, set
       this only for the drawFrame() call. */
    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
                                   GL::Renderer::BlendEquation::Add);
    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
                                   GL::Renderer::BlendFunction::OneMinusSourceAlpha);

    /* Loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
    _timeline.start();
}

Float GuiApplication::depthAt(const Vector2i& windowPosition) {
    /* First scale the position from being relative to window size to being
       relative to framebuffer size as those two can be different on HiDPI
       systems */
    const Vector2i position = windowPosition*Vector2{framebufferSize()}/Vector2{windowSize()};
    const Vector2i fbPosition{position.x(), GL::defaultFramebuffer.viewport().sizeY() - position.y() - 1};

    GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::Front);
    Image2D data = GL::defaultFramebuffer.read(
            Range2Di::fromSize(fbPosition, Vector2i{1}).padded(Vector2i{2}),
            {GL::PixelFormat::DepthComponent, GL::PixelType::Float});

    return Math::min<Float>(Containers::arrayCast<const Float>(data.data()));
}

Vector3 GuiApplication::unproject(const Vector2i& windowPosition, Float depth) const {
    /* We have to take window size, not framebuffer size, since the position is
       in window coordinates and the two can be different on HiDPI systems */
    const Vector2i viewSize = windowSize();
    const Vector2i viewPosition{windowPosition.x(), viewSize.y() - windowPosition.y() - 1};
    const Vector3 in{2*Vector2{viewPosition}/Vector2{viewSize} - Vector2{1.0f}, depth*2.0f - 1.0f};

    return _controller._camera->projectionMatrix().inverted().transformPoint(in);
}

void GuiApplication::keyPressEvent(KeyEvent& event) {
    if(_imgui.handleKeyPressEvent(event)) return;

    switch (event.key()){
        case KeyEvent::Key::W:
        case KeyEvent::Key::Up: {
            Float delta = _controller._key_move_velosity;
            const Vector3 p = Vector3(0, -delta, 0);
            _controller._cameraObject.translateLocal(p);
            _controller._translationPoint += p;
            redraw();
            break;
        }
        case KeyEvent::Key::S:
        case KeyEvent::Key::Down: {
            Float delta = _controller._key_move_velosity;
            const Vector3 p = Vector3(0, delta, 0);
            _controller._cameraObject.translateLocal(p);
            _controller._translationPoint += p;
            redraw();
            break;
        }
        case KeyEvent::Key::A:
        case KeyEvent::Key::Left: {
            Float delta = _controller._key_move_velosity;
            const Vector3 p = Vector3(delta, 0, 0);
            _controller._cameraObject.translateLocal(p);
            _controller._translationPoint += p;
            redraw();
            break;
        }
        case KeyEvent::Key::D:
        case KeyEvent::Key::Right:{
            Float delta = _controller._key_move_velosity;
            const Vector3 p = Vector3(-delta, 0, 0);
            _controller._cameraObject.translateLocal(p); /* is Z always 0? */
            _controller._translationPoint += p;
            redraw();
            break;
        }
        case KeyEvent::Key::Esc:
            exit(0);
        default:
            break;
    }

}

void GuiApplication::keyReleaseEvent(KeyEvent& event){
    if(_imgui.handleKeyReleaseEvent(event)) return;
}

void GuiApplication::mousePressEvent(MouseEvent& event) {
    if(_imgui.handleMousePressEvent(event)) return;
    /* Due to compatibility reasons, scroll is also reported as a press event,
       so filter that out. Could be removed once MouseEvent::Button::Wheel is
       gone from Magnum. */
    if(event.button() != MouseEvent::Button::Left &&
       event.button() != MouseEvent::Button::Middle)
        return;

    const Float currentDepth = depthAt(event.position());
    const Float depth = currentDepth == 1.0f ? _controller._lastDepth : currentDepth;
    _controller._translationPoint = unproject(event.position(), depth);
    /* Update the rotation point only if we're not zooming against infinite
       depth or if the original rotation point is not yet initialized */
    if(currentDepth != 1.0f || _controller._rotationPoint.isZero()) {
        _controller._rotationPoint = _controller._translationPoint;
        _controller._lastDepth = depth;
    }
}

void GuiApplication::mouseReleaseEvent(MouseEvent& event){
    if(_imgui.handleMouseReleaseEvent(event)) return;

}

void GuiApplication::mouseMoveEvent(MouseMoveEvent& event) {
    if(_imgui.handleMouseMoveEvent(event)) return;

    if(_controller._lastPosition == Vector2i{-1}) _controller._lastPosition = event.position();
    const Vector2i delta = event.position() - _controller._lastPosition;
    _controller._lastPosition = event.position();

    if(!event.buttons()) return;

    /* Translate */
    if(event.modifiers() & MouseMoveEvent::Modifier::Shift) {
        const Vector3 p = unproject(event.position(), _controller._lastDepth);
        _controller._cameraObject.translateLocal(_controller._translationPoint - p); /* is Z always 0? */
        _controller._translationPoint = p;

        /* Rotate around rotation point */
    } else _controller._cameraObject.transformLocal(
                Matrix4::translation(_controller._rotationPoint)*
                Matrix4::rotationX(-0.01_radf*delta.y())*
                Matrix4::rotationY(-0.01_radf*delta.x())*
                Matrix4::translation(-_controller._rotationPoint));

    redraw();
}

void GuiApplication::mouseScrollEvent(MouseScrollEvent& event) {
    if(_imgui.handleMouseScrollEvent(event)) {
        /* Prevent scrolling the page */
        event.setAccepted();
        return;
    }

    const Float currentDepth = depthAt(event.position());
    const Float depth = currentDepth == 1.0f ? _controller._lastDepth : currentDepth;
    const Vector3 p = unproject(event.position(), depth);
    /* Update the rotation point only if we're not zooming against infinite
       depth or if the original rotation point is not yet initialized */
    if(currentDepth != 1.0f || _controller._rotationPoint.isZero()) {
        _controller._rotationPoint = p;
        _controller._lastDepth = depth;
    }

    const Float direction = event.offset().y();
    if(!direction) return;

    /* Move towards/backwards the rotation point in cam coords */
    _controller._cameraObject.translateLocal(_controller._rotationPoint*direction*0.1f);

    event.setAccepted();
    redraw();
}

void GuiApplication::textInputEvent(TextInputEvent& event){
    if(_imgui.handleTextInputEvent(event)) return;
}

void GuiApplication::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});

    _imgui.relayout(event.windowSize());/*Vector2{event.windowSize()}/event.dpiScaling(),
                    event.windowSize(), event.framebufferSize());*/
}

void GuiApplication::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);
    readBackendData();
    setBackendData();

    _controller._camera->draw(_controller._drawables);
    //GL::Renderer::setPolygonMode(GL::Renderer::PolygonMode::Fill);
    _imgui.newFrame();

    /* Enable text input, if needed */
    if(ImGui::GetIO().WantTextInput && !isTextInputActive())
        startTextInput();
    else if(!ImGui::GetIO().WantTextInput && isTextInputActive())
        stopTextInput();

    if (_imguiState._showMeshStateWindow){
        //ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
        showGuiWindow();
    }

    /* Update application cursor */
    _imgui.updateApplicationCursor(*this);

    /* Set appropriate states. If you only draw ImGui, it is sufficient to
       just enable blending and scissor test in the constructor. */
    GL::Renderer::enable(GL::Renderer::Feature::Blending);
    GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
//    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::DepthTest);

    _imgui.drawFrame();

    /* Reset state. Only needed if you want to draw something else with
       different state after. */
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
//    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
    GL::Renderer::disable(GL::Renderer::Feature::Blending);

    swapBuffers();
    _timeline.nextFrame();

    redraw();
}

void GuiApplication::readBackendData() {
    g_interconnection.read_rendered_object_data(PartInterconnect::Gui, _interData);
}
void GuiApplication::setBackendData(){
    for (auto it = _interData.begin(); it != _interData.end(); ++it){
        auto& obj = *it;
        auto state = obj.second.second;
        switch (state){
            case InterconnectDataState::NotModified: {
                break;
            }
            case InterconnectDataState::Modified: {
                auto& self = obj.second.first;
                auto it = _meshes.find(obj.first);
                if (self.face.second)
                    it->second.updateIndexBuffer(self);
                if (self.vertex.second || self.value.second)
                    it->second.updateVertexBuffer(self);
                break;
            }
            case InterconnectDataState::New:{
                auto& self = obj.second.first;
                MeshInstance mesh(self, CustomDrawable::ShaderType::MeshLabel, CustomDrawable::RenderMode::WIREFRAME);
                mesh.addShaderType(self, CustomDrawable::ShaderType::MeshWireframe, CustomDrawable::RenderMode::WIREFRAME);
                mesh.addShaderType(self, CustomDrawable::ShaderType::PhongGL, CustomDrawable::RenderMode::FILL);

                auto p = _meshes.insert(pair<int, MeshInstance>{obj.first, move(mesh)});
                ((*p.first).second).drawable3D = make_unique<CustomDrawable>(_controller._manipulator, _shaders, (*p.first).second, _controller._drawables);
                break;
            }
            case InterconnectDataState::Deleted:{
                auto it = _meshes.find(obj.first);
                if (it == _meshes.end()) return;
                auto& mesh = it->second;
                mesh.drawable3D.reset(nullptr);
                _meshes.erase(it);
            }
        }
        if (state != InterconnectDataState::Deleted)
            obj.second.second = InterconnectDataState::NotModified;
    }
}

GuiApplication::~GuiApplication() {
    //The pointer referenced by this object will be deleted by the internal destructor.
    for (auto& i: _meshes)
        i.second.drawable3D.release();
}

void GuiApplication::addShaderType(MeshInstance& mesh, RenderedObject& ro, int shaderType, int renderType) {
    mesh.addShaderType(ro, shaderType, renderType);
}

void GuiApplication::showGuiWindow() {
    ImGui::Begin("Thin shell real-time models viewer", nullptr);
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0/Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    if (ImGui::CollapsingHeader("Main scene")){
        if (ImGui::ColorEdit4("Clear color", _controller._clearColor.data(), ImGuiColorEditFlags_NoInputs))
            GL::Renderer::setClearColor(_controller._clearColor);
        ImGui::DragFloat("Move velocity", &_controller._key_move_velosity, 0.005f);
    }
    if (ImGui::CollapsingHeader("Mesh control")){
        static int item_current_idx = -1;
        string label = "";
        if (item_current_idx > 0){
            auto it = _meshes.find(item_current_idx);
            if (it != _meshes.end()){
                auto nm = _interData[it->first].first.name;
                if (nm != "") label = nm;
                else label = to_string(it->first);
            }
        }
        if (ImGui::BeginCombo("objects", label.c_str())){
            for (auto& i: _interData){
                const bool is_selected = (item_current_idx == i.first);
                const string& name = (i.second.first.name == "") ? to_string(i.first) : i.second.first.name;
                if (ImGui::Selectable(name.c_str(), is_selected))
                    item_current_idx = i.first;

                // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (auto it = _meshes.find(item_current_idx); item_current_idx > 0 && it != _meshes.end()){
            auto& obj = it->second;
            ImGui::Text("Mesh info: ");
            ImGui::Text("\t vertices: %d", obj._nvertex);
            ImGui::Text("\t faces:    %d", obj._nfaces);
            ImGui::Text("\t edges:    %d", obj._nedges);
            ImGui::ColorEdit4("Main/diffuse color", obj.color.data(), ImGuiColorEditFlags_NoInputs);
            ImGui::SameLine();
            ImGui::ColorEdit4("Ambient color", obj.ambientColor.data(), ImGuiColorEditFlags_NoInputs);
            ImGui::ColorEdit4("Wireframe color", obj.wireframeColor.data(), ImGuiColorEditFlags_NoInputs);

            using CDRM = CustomDrawable::RenderMode;
            using CDST = CustomDrawable::ShaderType;
            const string shaders_item[CDST::NumShaderType] = {
                                            "MeshWireframe",
                                            "PhongGL",
                                            "PhongGLVColor",
                                            "MeshNormal",
                                            "MeshLabel"
            };
            const string render_items[CDRM::NumRENDERMODE] = {
                    "Fill",
                    "Wireframe"
            };
            const int dif = CustomDrawable::RenderMode::NumRENDERMODE;
            auto renderEditorWindow = [shaders_item, render_items, dif](int& shader_type, string name = "Shader"){
                ImGui::OpenPopup(name.c_str());
                char is_open = 1;
                if (ImGui::BeginPopupModal(name.c_str(), reinterpret_cast<bool*> (&is_open))) {
                    ImGui::Text("                               ");
                    ImGui::Columns(2, NULL, true);

                    int type = (shader_type >= 0) ? shader_type / dif : shader_type;
                    for (int i = 0; i < CustomDrawable::ShaderType::NumShaderType; i++) {
                        if (ImGui::Selectable(shaders_item[i].c_str(), type == i))
                            type = i;
                    }
                    ImGui::NextColumn();
                    int rend = (shader_type >= 0) ? shader_type % dif : shader_type;
                    for (int i = 0; i < CustomDrawable::RenderMode::NumRENDERMODE; i++) {
                        if (ImGui::Selectable(render_items[i].c_str(), rend == i))
                            rend = i;
                    }
                    ImGui::Columns(1);
                    if (type >= 0 && rend >= 0)
                        shader_type = type * dif + rend;

                    if (ImGui::Button("Save")) {
                        is_open = -1;
                        ImGui::CloseCurrentPopup();
                    }

                    ImGui::EndPopup();
                }
                return is_open;
            };
            ImGui::LabelText("", "Renderers");
            if (ImGui::BeginCombo("", "shaders included")){
                static int current_list_elem = -1;
                for (int ii = 0; ii < static_cast<int>(obj.shaderTypes.size()); ++ii){
                    int i = obj.shaderTypes[ii];
                    string item = shaders_item[i / dif] + string(" ") + render_items[i % dif];
                    ImGui::Selectable(item.c_str());

                    if (ImGui::IsItemActive() && !ImGui::IsItemHovered()){
                        int delta = ImGui::GetMouseDragDelta(0).y < 0.f ? -1 : 1;
                        int next = ii + delta;
                        if (next >= 0 && next < static_cast<int>(obj.shaderTypes.size())){
                            i = obj.shaderTypes[next];
                            swap(obj.shaderTypes[ii], obj.shaderTypes[next]);
                            item = shaders_item[i / dif] + string(" ") + render_items[i % dif];
                            ImGui::ResetMouseDragDelta();
                        }
                    }
                    if (ImGui::IsItemHovered()) {
                        current_list_elem = ii;
                    }
                    static bool edit_stat = false;
                    static int state_st;
                    if (ii == current_list_elem && ImGui::BeginPopupContextItem("object shader's context menu"))
                    {
                        if (ImGui::Selectable("Delete")) {
                            ii -= obj.deleteShaderType(_interData[it->first].first, i / dif, i % dif);
                        }
                        if (ImGui::Selectable("Edit")){
                            edit_stat  = true;
                            state_st = obj.shaderTypes[ii];
                        }
                        ImGui::EndPopup();
                    }
                    if (edit_stat && ii == current_list_elem){
                        char stat = renderEditorWindow(state_st);
                        edit_stat = (stat > 0);
                        if (!edit_stat && state_st != obj.shaderTypes[ii]) {
                            int type = state_st;
                            obj.deleteShaderType(_interData[it->first].first, obj.shaderTypes[ii] / dif,
                                                 obj.shaderTypes[ii] % dif);
                            obj.addShaderType(_interData[it->first].first, type / dif, type % dif);
                            int tmp = obj.shaderTypes[obj.shaderTypes.size() - 1];
                            for (int l = obj.shaderTypes.size() - 1; l >= ii; --l)
                                obj.shaderTypes[l] = obj.shaderTypes[l - 1];
                            obj.shaderTypes[ii] = tmp;
                        }
                    }
                }
                ImGui::EndCombo();
            } ImGui::SameLine();
            static bool clicked_add_render_type = false;
            static int new_type = -1;
            if (ImGui::Button("Add")) {
                clicked_add_render_type = true;
                new_type = 0;
            }
            if (clicked_add_render_type){
                int stat = renderEditorWindow(new_type, "New shader");
                clicked_add_render_type = (stat > 0);
                if (stat == -1 && new_type >= 0){
                    obj.addShaderType(_interData[it->first].first, new_type / dif, new_type % dif);
                }
            }

        }

    }
    ImGui::End();
}

void GuiApplication::MeshInstance::updateVertexRequire() {
    oldVertexRequire = _vertexRequire;
    for (auto i: shaderTypes){
        int shaderType = i / CustomDrawable::RenderMode::NumRENDERMODE;
        switch (shaderType) {
            case CustomDrawable::PhongGL: {
                _vertexRequire |= Vertex|Normal;
                break;
            }
            case CustomDrawable::PhongGLVColor:{
                _vertexRequire |= Vertex|Normal|Color;
                break;
            }
            case CustomDrawable::MeshWireframe:{
                _vertexRequire |= Vertex;
                break;
            }
            case CustomDrawable::MeshNormal:{
                _vertexRequire |= Vertex|Normal;
                break;
            }
            case CustomDrawable::MeshLabel:{
                break;
            }
            default: throw runtime_error("VertexRequire for such shader type = " + to_string(shaderType) + " is not specified");
        }
    }
}

void GuiApplication::MeshInstance::updateVertexBuffer(RenderedObject &ro, bool forced) {
    if (!_vertexRequire) return;
    _nvertex = ro.vertex.first.size();
    int stride = ((_vertexRequire & Vertex) > 0) * sizeof(Vector3) + ((_vertexRequire & Normal) > 0) * sizeof(Vector3) + ((_vertexRequire & Color) > 0) * sizeof(Color4);
    int shift = 0;
    auto& verts = ro.vertex.first;

    vector<char> vdata(verts.size()*stride);
    char* data = vdata.data();

    if (_vertexRequire & Vertex && (ro.vertex.second || forced)) {
        for (unsigned i = 0; i < verts.size(); ++i) {
            *(reinterpret_cast<Vector3 *> (data + i * stride + shift)) = Vector3(verts[i][0], verts[i][1], verts[i][2]);
        }
        shift += sizeof(Vector3);
    }
    if (_vertexRequire & Normal && (ro.vertex.second || forced)) {
        for (unsigned i = 0; i < verts.size(); ++i) {
            *(reinterpret_cast<Vector3 *> (data + i * stride + shift)) = Vector3(0, 0, 0);
        }
        for (unsigned i = 0; i < ro.face.first.size() / 3; ++i) {
            Vector3 p[3], d[3], n;
            Float w[3], a[3];
            for (unsigned j = 0; j < 3; ++j) {
                int ind = ro.face.first[3 * i + j];
                p[j] = Vector3(verts[ind][0], verts[ind][1], verts[ind][2]);
            }
            for (unsigned j = 0; j < 3; ++j) {
                d[j] = p[(j + 2) % 3] - p[(j + 1) % 3];
                a[j] = d[j].length();
            }
            n = Math::cross(d[2], -d[1]);
            float len = n.length();
            w[0] = len / (a[2] * a[1]);
            if (len >= 1.0e-8f) n /= len;
            w[1] = a[1] * w[0] / a[0];
            w[2] = a[2] * w[0] / a[0];
            for (unsigned j = 0; j < 3; ++j) {
                if (w[j] > 1.0f) w[j] = 1.0f;
                if (w[j] < -1.0f) w[j] = -1.0f;
                w[j] = float(Math::asin(w[j]));
                *(reinterpret_cast<Vector3 *> ( data + ro.face.first[3 * i + j] * stride + shift)) += w[j] * n;
            }

        }
        shift += sizeof(Vector3);
    }
    if (_vertexRequire & Color && (ro.value.second || forced)){
        pair<Float, Float> range;
        if (valueRange.first != valueRange.second)
            range = valueRange;
        else {
            auto p = minmax_element(ro.value.first.begin(), ro.value.first.end());
            if (p.first.base() && p.second.base())
                range = pair<Float, Float>{*p.first, *p.second};
            else
                range = pair<Float, Float>{*ro.value.first.begin(), *ro.value.first.begin() + 1};
        }

        Float len = range.second - range.first;
        if (len < FLT_EPSILON) len = 1.0f;
        for (unsigned i = 0; i < ro.vertex.first.size(); ++i){
            float approx = max(min((1 - (ro.value.first[i] - range.first) / len), 1.0f), 0.0f);
            ColorHsv c(approx * 240.0_degf, 1.0, 1.0);
            *(reinterpret_cast<Color4*> (data+i*stride + shift)) = Color4::fromHsv(c);
        }
        shift += sizeof(Color4);
    }
    vertexBuffer.setData(move(vdata), Magnum::GL::BufferUsage::DynamicDraw);
    vdata.clear();

    if (ro.vertex_lbl.second && ro.vertex_lbl.first.size()){
        struct LabeledVertex{
            Vector3 v;
            Color4 c;
        };
        vector<LabeledVertex> lv(ro.vertex_lbl.first.size());
        for (unsigned i = 0; i < ro.vertex_lbl.first.size(); ++i){
            int j = ro.vertex_lbl.first[i].first;
            lv[i] = {Vector3(ro.vertex.first[j][0], ro.vertex.first[j][1], ro.vertex.first[j][2]),
                     Color4::fromHsv(ColorHsv(ro.vertex_lbl.first[i].second * 1.0f / _lbl_spliter.second * 240.0_degf, 1.0, 1.0))
                     };
        }
        lbl_vertexBuffer.setData(lv);
    }
}

MeshIndexType GuiApplication::MeshInstance::updateIndexBuffer(RenderedObject &ro, bool forced) {
    _nfaces = ro.face.first.size() / 3;
    if (!ro.face.second && !forced) return static_cast<MeshIndexType>(mesh.indexType());
    mesh.setCount(ro.face.first.size());
    std::pair<Containers::Array<char>, MeshIndexType> compressed =
            MeshTools::compressIndices(ro.face.first, 0);
    indexBuffer.setData(compressed.first);

    if (forced || mesh.indexTypeSize() != meshIndexTypeSize(compressed.second)) {
        mesh.setIndexBuffer(indexBuffer, 0, compressed.second);
    }

    set<pair<int, int>> edges;
    for (unsigned i = 0; i < ro.face.first.size() / 3; ++i){
        edges.insert(minmax(ro.face.first[3*i], ro.face.first[3*i+1]));
        edges.insert(minmax(ro.face.first[3*i+1], ro.face.first[3*i+2]));
        edges.insert(minmax(ro.face.first[3*i+2], ro.face.first[3*i]));
    }
    _nedges = edges.size();
    edges.clear();

    if (ro.vertex_lbl.second && ro.vertex_lbl.first.size()){
        lbl_mesh.second = true;
        set<pair<int, int>> edges;
        for (unsigned i = 0; i < ro.face.first.size() / 3; ++i){
            edges.insert(minmax(ro.face.first[3*i], ro.face.first[3*i+1]));
            edges.insert(minmax(ro.face.first[3*i+1], ro.face.first[3*i+2]));
            edges.insert(minmax(ro.face.first[3*i+2], ro.face.first[3*i]));
        }
        vector<int> verts(ro.vertex.first.size(), 0);
        pair<int, int>& lbl_spliter = _lbl_spliter;
        for (auto& i: ro.vertex_lbl.first) {
            verts[i.first] = i.second;
            if (i.second > lbl_spliter.second)
                lbl_spliter.second = i.second;
        }
        vector<UnsignedInt> edge_dt;
        edge_dt.reserve(ro.vertex.first.size() + 1);
        for (auto& i: edges){
            if (verts[i.first] && verts[i.second]) {
                edge_dt.push_back(i.first);
                edge_dt.push_back(i.second);
            }
        }
        edges.clear();
        vector<int> loc_index(ro.vertex.first.size(), -1);
        for (unsigned i = 0; i < ro.vertex_lbl.first.size(); ++i)
            loc_index[ro.vertex_lbl.first[i].first] = i;
        for (auto& i: edge_dt) i = loc_index[i];
        int indexCount = edge_dt.size();
        compressed = MeshTools::compressIndices(move(edge_dt), 0);
        lbl_indexBuffer.setData(compressed.first);
        lbl_mesh.first.setCount(indexCount)
                .setIndexBuffer(lbl_indexBuffer, 0, compressed.second)
                .addVertexBuffer(lbl_vertexBuffer, 0, Shaders::GenericGL3D::Position{}, Shaders::GenericGL3D::Color4{});
    }

    return compressed.second;
}

void GuiApplication::MeshInstance::addShaderType(RenderedObject& ro, int shaderType, int renderMode, bool forceIndex){
    int type = shaderType * CustomDrawable::RenderMode::NumRENDERMODE + renderMode;
    if (auto it = shaderTypes.push_back_unique(type); it == shaderTypes.end()) return;
    updateVertexRequire();
    updateVertexBuffer(ro, true);
    updateIndexBuffer(ro, forceIndex);
    addVertexBuffer();
}

GuiApplication::MeshInstance::MeshInstance(RenderedObject &ro, int shaderType, int renderMode) {
    mesh.setPrimitive(MeshPrimitive::Triangles);
    lbl_mesh.second = false;
    lbl_mesh.first.setPrimitive(MeshPrimitive::Lines);
    addShaderType(ro, shaderType, renderMode, true);
}

int GuiApplication::MeshInstance::deleteShaderType(RenderedObject &ro, int shaderType, int renderMode) {
    int ret = 0;
    int type = shaderType * CustomDrawable::RenderMode::NumRENDERMODE + renderMode;
    if (auto it = find(shaderTypes.begin(), shaderTypes.end(), type); it != shaderTypes.end()) {
        shaderTypes.erase(it);
        updateVertexRequire();
        ++ret;
    }
    if (oldVertexRequire != _vertexRequire){
        updateVertexBuffer(ro, true);
        addVertexBuffer();
    }
    return ret;
}

void GuiApplication::MeshInstance::addVertexBuffer() {
    if (oldVertexRequire != _vertexRequire || !_vertexRequire) {
        switch (_vertexRequire) {
            case 0x7:
                mesh.addVertexBuffer(vertexBuffer, 0, Shaders::GenericGL3D::Position{}, Shaders::GenericGL3D::Normal{},
                                     Shaders::GenericGL3D::Color4{});
                break;
            case 0x5:
                mesh.addVertexBuffer(vertexBuffer, 0, Shaders::GenericGL3D::Position{}, Shaders::GenericGL3D::Normal{});
                break;
            case 0x1:
                mesh.addVertexBuffer(vertexBuffer, 0, Shaders::MeshVisualizerGL3D::Position{});
                break;
            case 0x3:
                mesh.addVertexBuffer(vertexBuffer, 0, Shaders::GenericGL3D::Position{}, Shaders::GenericGL3D::Color4{});
                break;
            case 0x0:
                mesh.addVertexBuffer(vertexBuffer, 0);
                break;
            default:
                throw runtime_error("addVertex is not specified for such _vertexRequire yet");
        };
        oldVertexRequire = _vertexRequire;
    }
}

void GuiApplication::SceneController::init() {
    _cameraObject
            .setParent(&_scene)
            .translate(Vector3::zAxis(5.0f));
    (*(_camera = new SceneGraph::Camera3D{_cameraObject}))
            .setAspectRatioPolicy(SceneGraph::AspectRatioPolicy::Extend)
            .setProjectionMatrix(Matrix4::perspectiveProjection(35.0_degf, 1.0f, 0.01f, 1000.0f))
            .setViewport(GL::defaultFramebuffer.viewport().size());
    /* Base object, parent of all (for easy manipulation) */
    _manipulator.setParent(&_scene);

    _lastDepth = ((_camera->projectionMatrix()*_camera->cameraMatrix()).transformPoint({}).z() + 1.0f)*0.5f;
}

vector<unique_ptr<GL::AbstractShaderProgram>> CustomDrawable::init_shaders() {
    vector<unique_ptr<GL::AbstractShaderProgram>> res;
    Color4 color = 0x2f83cc_rgbf;
    res.push_back(make_unique<Shaders::MeshVisualizerGL3D>(Shaders::MeshVisualizerGL3D::Flag::Wireframe));
    static_cast<Shaders::MeshVisualizerGL3D*>(res[0].get())->setColor(0x2f83cc_rgbf)
            .setWireframeColor(0xdcdcdc_rgbf)
            .setViewportSize(Vector2{GL::defaultFramebuffer.viewport().size()});
    res.push_back(make_unique<Shaders::PhongGL>());
    static_cast<Shaders::PhongGL*>(res[1].get())->setLightPositions({{7.0f, 5.0f, 2.5f, 0.0f}})
            .setLightColors({Color3{1.0f}})
            .setDiffuseColor(0x000000_rgbf)
            .setAmbientColor(color)
            .setShininess(80.0f);
    res.push_back(make_unique<Shaders::PhongGL>(Shaders::PhongGL::Flag::VertexColor));
    static_cast<Shaders::PhongGL*>(res[2].get())->setLightPositions({{7.0f, 5.0f, 2.5f, 0.0f}})
            .setLightColors({Color3{1.0f}})
            .setDiffuseColor(color)
            .setAmbientColor(Color3::fromHsv({color.hue(), 1.0f, 0.3f}))
            .setShininess(80.0f);
    res.push_back(make_unique<Shaders::MeshVisualizerGL3D>(Shaders::MeshVisualizerGL3D::Flag::Wireframe|Shaders::MeshVisualizerGL3D::Flag::NormalDirection));
    static_cast<Shaders::MeshVisualizerGL3D*>(res[3].get())->setViewportSize(Vector2{GL::defaultFramebuffer.viewport().size()})
            .setLineLength(0.3f);

    res.push_back(make_unique<Shaders::VertexColorGL3D>());
    static_cast<Shaders::VertexColorGL3D*>(res[4].get());
    return res;
}

void CustomDrawable::draw(const Matrix4 &transformationMatrix, SceneGraph::Camera3D &camera) {
    for (auto type: _mesh.shaderTypes){
        int renderMode = type % NumRENDERMODE,
                shaderType = type / NumRENDERMODE;
        switch (renderMode) {
            case FILL: {
                GL::Renderer::setPolygonMode(GL::Renderer::PolygonMode::Fill);
                break;
            }
            case WIREFRAME:{
                GL::Renderer::setPolygonMode(GL::Renderer::PolygonMode::Line);
                break;
            }
            default:
                throw runtime_error("such renderMode = " + to_string(renderMode) + " is not implemented yet");
        }
        switch (shaderType) {
            case PhongGL:{
                (*static_cast<Shaders::PhongGL*>(_shaders[PhongGL].get()))
                        .setAmbientColor(_mesh.ambientColor)
                        .setDiffuseColor((renderMode == WIREFRAME) ? _mesh.wireframeColor : _mesh.color)
                        .setTransformationMatrix(transformationMatrix)
                        .setNormalMatrix(transformationMatrix.normalMatrix())
                        .setProjectionMatrix(camera.projectionMatrix())
                        .draw(_mesh.mesh);
                break;
            }
            case PhongGLVColor:{
                (*static_cast<Shaders::PhongGL*>(_shaders[PhongGLVColor].get()))
                        .setTransformationMatrix(transformationMatrix)
                        .setNormalMatrix(transformationMatrix.normalMatrix())
                        .setProjectionMatrix(camera.projectionMatrix())
                        .draw(_mesh.mesh);
                break;
            }
            case MeshWireframe:{
                (*static_cast<Shaders::MeshVisualizerGL3D*>(_shaders[MeshWireframe].get()))
                        .setColor(_mesh.color)
                        .setWireframeColor(_mesh.wireframeColor)
                        .setTransformationMatrix(transformationMatrix)
                        .setProjectionMatrix(camera.projectionMatrix())
                        .draw(_mesh.mesh);
                break;
            }
            case MeshNormal:{
                (*static_cast<Shaders::MeshVisualizerGL3D*>(_shaders[MeshNormal].get()))
                        .setColor(_mesh.color)
                        .setWireframeColor(_mesh.wireframeColor)
                        .setTransformationMatrix(transformationMatrix)
                        .setProjectionMatrix(camera.projectionMatrix())
                        .setNormalMatrix(transformationMatrix.normalMatrix())
                        .draw(_mesh.mesh);
                break;
            }
            case MeshLabel:{
                if (_mesh.lbl_mesh.second) {
                    GL::Renderer::setDepthFunction(GL::Renderer::DepthFunction::LessOrEqual);
                    (*static_cast<Shaders::VertexColorGL3D *>(_shaders[MeshLabel].get()))
                            .setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
                            .draw(_mesh.lbl_mesh.first);
                    GL::Renderer::setDepthFunction(GL::Renderer::DepthFunction::Less);
                }
                break;
            }
            default:
                throw runtime_error("such shaderType = " + to_string(shaderType) + " is not implemented yet");
        }
        GL::Renderer::setPolygonMode(GL::Renderer::PolygonMode::Fill);

    }
}