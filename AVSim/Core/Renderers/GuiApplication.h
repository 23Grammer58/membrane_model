//
// Created by alex on 21.06.2020.
//

#ifndef AORTIC_VALVE_GUIAPPLICATION_H
#define AORTIC_VALVE_GUIAPPLICATION_H

#ifdef USE_MAGNUM_GUI
#include <vector>
#include <memory>
#include <algorithm>
#include <set>
#include <Corrade/Containers/ArrayViewStl.h>
#include <Corrade/Containers/Array.h>
//#include <Corrade/Containers/Optional.h>
//#include <Corrade/PluginManager/Manager.h>
#include <Corrade/Utility/Arguments.h>
//#include <Corrade/Utility/DebugStl.h>
//#include <Corrade/Containers/GrowableArray.h>
//#include <Corrade/Containers/Pointer.h>
#include <Magnum/Image.h>
//#include <Magnum/ImageView.h>
//#include <Magnum/Mesh.h>
#include <Magnum/GL/DefaultFramebuffer.h>
//#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/PixelFormat.h>
//#include <Magnum/PixelFormat.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Math/Color.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <Magnum/GL/Buffer.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Capsule.h>
#include <Magnum/Primitives/Cone.h>
#include <Magnum/GL/TextureFormat.h>

#include "GuiBackInterconnection.h"


using namespace Magnum;
using namespace Math::Literals;
using namespace std;

template <class T>
class unique_vector: public vector<T>{
public:
    typename vector<T>::iterator push_back_unique( const T& value ){
        auto it = find(vector<T>::begin(), vector<T>::end(), value);
        if (it == vector<T>::end()) {
            vector<T>::push_back(value);
            return --vector<T>::end();
        }
        return vector<T>::end();
    }
};


typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> MagnumObject3D;
typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D> Scene3D;

class GuiApplication: public Platform::Application {
public:
    struct MeshInstance {
        enum VertexRequire{
            None = 0,
            Vertex = 1,
            Color = 2,
            Normal = 4
        };

        int _nedges = 0;
        int _nvertex = 0;
        int _nfaces = 0;

        GL::Mesh mesh;
        GL::Buffer vertexBuffer;
        GL::Buffer indexBuffer;
        pair<GL::Mesh, bool> lbl_mesh;
        pair<int, int> _lbl_spliter = {0, 1};
        GL::Buffer lbl_vertexBuffer;
        GL::Buffer lbl_indexBuffer;

        Color4 color = 0x065091_rgbf;
        Color4 ambientColor = 0x536f87_rgbf;
        Color4 wireframeColor = 0x11ebde_rgbf;
        unique_ptr<SceneGraph::Drawable3D> drawable3D{nullptr};
        unique_vector<int> shaderTypes;
        int _vertexRequire = None;
        int oldVertexRequire = None;
        pair<Float, Float> valueRange = {0.0f, 0.0f};

        void updateVertexRequire();
        void updateVertexBuffer(RenderedObject& ro, bool forced = false);
        MeshIndexType updateIndexBuffer(RenderedObject& ro, bool forced = false);

        void addShaderType(RenderedObject& ro, int shaderType, int renderMode = 0, bool forceIndex = false);
        int deleteShaderType(RenderedObject& ro, int shaderType, int renderMode);
        MeshInstance(RenderedObject& ro, int shaderType = 0, int renderMode = 0);
        MeshInstance() {};
    private:
        void addVertexBuffer();
    };

    struct SceneController{
        Color4 _clearColor = 0x72909aff_rgbaf;

        Scene3D _scene;
        MagnumObject3D _manipulator, _cameraObject;
        SceneGraph::Camera3D* _camera;
        SceneGraph::DrawableGroup3D _drawables;
        //Vector3 _previousPosition;

        Float _lastDepth;
        Vector2i _lastPosition{-1};
        Vector3 _rotationPoint, _translationPoint;

        Float _key_move_velosity = 0.05;

        void init();
    };

    struct ImGuiState{
        bool _showMeshStateWindow = true;
        bool _showDemoWindow = true;
        bool _showAnotherWindow = false;
        Color4 _clearColor = 0x72909aff_rgbaf;
        Float _floatValue = 0.0f;
    };

public:
    explicit GuiApplication(const Arguments& arguments);

    void drawEvent() override;

    void viewportEvent(ViewportEvent& event) override;

    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;

    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void mouseMoveEvent(MouseMoveEvent& event) override;
    void mouseScrollEvent(MouseScrollEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;

    virtual ~GuiApplication();

private:
    Float depthAt(const Vector2i& windowPosition);
    Vector3 unproject(const Vector2i& windowPosition, Float depth) const;
    void showGuiWindow();

    void readBackendData();
    void setBackendData();
    map<int, pair<RenderedObject, InterconnectDataState>> _interData;


    ImGuiIntegration::Context _imgui{NoCreate};
    ImGuiState _imguiState;

    vector<unique_ptr<GL::AbstractShaderProgram>> _shaders;

    map<int, MeshInstance> _meshes;
    void addShaderType(MeshInstance& mesh, RenderedObject& ro, int shaderType, int renderType);

    SceneController _controller;

    Timeline _timeline;
};

class CustomDrawable: public SceneGraph::Drawable3D {
public:
    enum ShaderType{
        MeshWireframe,
        PhongGL,
        PhongGLVColor,
        MeshNormal,
        MeshLabel,
        NumShaderType
    };
    enum RenderMode{
        FILL,
        WIREFRAME,
        NumRENDERMODE
    };

    static vector<unique_ptr<GL::AbstractShaderProgram>> init_shaders();

public:
    explicit CustomDrawable(MagnumObject3D& object, vector<unique_ptr<GL::AbstractShaderProgram>>& shaders, GuiApplication::MeshInstance& mesh, SceneGraph::DrawableGroup3D& group):
            SceneGraph::Drawable3D{object, &group}, _shaders{shaders}, _mesh{mesh} {}

private:
    void draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) override;

    vector<unique_ptr<GL::AbstractShaderProgram>>& _shaders;
    GuiApplication::MeshInstance& _mesh;
};

#endif

#endif //AORTIC_VALVE_GUIAPPLICATION_H
