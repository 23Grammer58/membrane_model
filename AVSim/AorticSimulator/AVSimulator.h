//
// Created by alex on 17.11.2020.
//

#ifndef AORTIC_VALVE_AVSIMULATOR_H
#define AORTIC_VALVE_AVSIMULATOR_H
#include <cstdint>
#include <string>
#include <thread>
#include "Configuration.h"
#include <fstream>
#include "AVSim/Core/World.h"
#include "AVSim/Core/Renderers/Renderer.h"
#include "AVSim/Core/Renderers/GuiApplication.h"
#include "AVSim/Core/Renderers/GuiBackInterconnection.h"
#include "Messurer.h"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif


//TODO: надо повернуть и центрировать входные данные

struct AorticSurface{
    uint32_t nT; //число треугольников на триангуляции поверхности
    uint32_t nP; //число вершин в триангуляции
    uint32_t*  tri;//[nT*3] - массив, где каждые следующие 3 числа задают индексы (начиная с нуля) вершин треугольника
    //треугольники должны быть ориентированы по часовой стрелке при взгляде на аорту снаружи (т.е. нормаль внутрь аорты)
    float* pnt; //[nP*3] - массив, типа {{x0, y0, z0}, {x1, y1, z1}, ... {xi, yi, zi} ...},
    // где каждые следующие 3 числа задают координаты новой вершины
};

AorticSurface fromMesh(const World3d::Mesh& m);
void freeAorticSurface(AorticSurface* as);

struct SewLines{
    struct Line{
        uint32_t nP; //число точек в линии
        float* pnt; //[nP*3] - последовательные вершины линии
    };
    uint8_t nL; //число линий
    Line* lines; //[nL] - набор из nL линий
    //в массиве pnt каждой линии первая и последняя точки это точки комиссуры
    //линии должны быть ориентированы так, чтобы выполнялось lines[i%3].pnt[lines[i%3].nP-1] = lines[(i + 1)%3].pnt[0]
    //кроме того, точки lines[0].pnt[0], lines[1].pnt[0], lines[2].pnt[0] должны быть расположены против часовой стрелки
    // при взгляде на треугольник, который они образуют, сверху аорты (из восходящей аорты)
};
void freeSewLines(SewLines* sl);
SewLines readSewLines(std::string file);

struct CoaptationResult{
    using Point = float[3];
    struct Billowing{
        Point billowPlate[3]; //плоскость, относительно которой измеряется прогиб
        //здесь 3 это число лепестков
        Point billowPnt[3]; //нижняя точка прогиба для каждого лепестка
        float billowing[3]; //величина прогиба для каждого лепестка
        //эти данные можно отображать рисуя по точкам из billowPlate плоскость и рисуя точки billowPnt на лепестках
        //и отрисовывая линию расстояния от точки billowPnt[i] до плоскости billowPlate с подписью,
        //что это расстояние равно billowing[i]
    };
    struct CoaptDistribution{
        uint32_t N; //число точек распределения (первая точка задаёт один край лепестка, а последняя - второй край)
        float width; //ширина области распределения
        float* bottom; //[N] - нижний край распределения коаптации
        float* top; //[N] - верхний край распределения
        //для демонстрации этих данных следует строить диаграмму с изображёнными на ней линиями top и bottom,
        //разница между top[i] и bottom[i] задаёт величину коптации в i-ом сечении по длине лепестка
    };
    Billowing bill;
    CoaptDistribution cd[3];
    float centralCoaptation; //величина центральной коаптации
    float maxCoaptationHeight; //высота коаптации
    uint8_t coaptStatus;
    //if (coaptStatus & (1 << 0)) than cusp[0] coapting with cusp[1]
    //if (coaptStatus & (1 << 1)) than cusp[0] coapting with cusp[2]
    //if (coaptStatus & (1 << 2)) than cusp[1] coapting with cusp[0]
    //if (coaptStatus & (1 << 3)) than cusp[1] coapting with cusp[2]
    //if (coaptStatus & (1 << 4)) than cusp[2] coapting with cusp[0]
    //if (coaptStatus & (1 << 5)) than cusp[2] coapting with cusp[1]
    //if (coaptStatus & (1 << 6)) than valve is closed (has points of triple coaptation)
};

CoaptationResult createCoaptationResult(int N);
void freeCoaptationResult(CoaptationResult* cr);

struct IntermediateResults{
    struct Leaflet {
        uint8_t state; //текущее состояние моделируемого лепестка:
        //    NotModified = 0 - данные лепестка не изменились с момента последнего считывания данных
        //    Modified = 1 - данные лепестка изменились с момента последнего считывания данных
        //    New = 2 - этот лепесток считывается впервые
        //    Deleted = 3 - этот лепесток более не участвует в моделировании и должен быть удалён из отображающихся
        uint8_t nm_sz; char* name; //[nm_sz] - имя лепестка аорты
        //далее в массивах этой структуры длина (_sz) 0 означает, что данный массив не изменён с момента последнего считывания данных
        uint32_t v_sz;   float*  vertex; //[v_sz*3] - массив координат всех вершин лепестка
        uint32_t val_sz; float*  value;  //[val_sz] - массив, где i-ое значение задаёт какое-то число на i-ой вершине
        uint32_t vlb_sz; int32_t*  v_label;//[vlb_sz] - массив, целочисленных значений лейблов вершин
        uint32_t f_sz;   uint32_t* faces;  //[f_sz*3] - массив, где каждые 3 числа задают индексы вершин отдельного треугольника
        struct String{
            uint32_t sz; //размер строки
            char* str; //строка
        };
        uint32_t av_sz;  String* availableValues; //[av_sz] - массив имён величин, которые можно запросить для сохранения в массив value
    };
    uint8_t nL; //количество лепестков
    Leaflet* leaflets;//[nL]
};

//g_interconnection from "GuiBackInterconnection.h" responsible for interconnection between my model and yours
struct InterConnector{
public:
    map<int, pair<RenderedObject, InterconnectDataState>> _interData;
    GuiBackInterconnection* g_ic = &g_interconnection;

    void readBackendData() {
        g_ic->read_rendered_object_data(PartInterconnect::Gui, _interData);
    }
    void sendBackEndData(std::string msg) {
        g_ic->send_msg(PartInterconnect::Gui, msg);
    }
    void update(IntermediateResults& ir);
    static void destructIntermediateResults(IntermediateResults& ir);
};

class AVSimulator {
public:
    using Mesh = World3d::Mesh;
    using Object3D = World3d::Object3D;
    using Vector = World3d::Vector;
    using Point = World3d::Point;
    using Point_2 = World3d::Point_2;

    struct SewLineData{
        using Node = Vector;
        using Direction = Vector;
        using Line = vector<Node>;
        using ConstrictionCDF = std::function<bool(std::vector<double>& val, double J0, double l_l)>;
        /*
         * val - input-output vector,
         *       on input contains sequentially increasing set of numbers ranging from 0 to 1
         *       on output contains sequentially increasing set of values of CDF function on values of vector
         * J0 - relation l_s/l_l
         * l_l - length of sewed boundary on template
         * l_s - length of sewed line on aorta (= J0 * l_l)
         */

        vector<Line> lines;
        vector<ConstrictionCDF> cdfs;
        std::array<Point, 3> commissure_points;
        Direction orientation; //Direction from the left ventricle to the ascending aorta along which the coaptation height will be measured
        static ConstrictionCDF getDefaultCDF(){ return [](std::vector<double>& val, double J0, double l_l){ return true; }; }
    };
    struct AorticSurfaceData{
        typedef std::vector<uint32_t> CGAL_Polygon;
        std::vector<array<float, 3>> points;
        std::vector<CGAL_Polygon> polygons;
    };

    SewLineData sld;
    AorticSurfaceData asd;

    World m_w;
    ObjectID m_aid;
    std::array<World3d::ObjectID, 3> m_lid;
    Configuration conf;

public:
    void setSewLines(const SewLines& sw);
    void setAorticSurface(const AorticSurface& as);
    void sewLeaflets(const Object3D& aorta, vector<Object3D*>& leafs);

    Object3D make_aorta_object(double long_extend_coef_up = 1.1, double long_extend_coef_down = 1.2, double radial_extend_coef = 1.2);
    void set_aorta_object(){
        Object3D aorta = make_aorta_object();
        m_aid = m_w.addObject3D(move(aorta), 0);
    }
    void read_leafs();
    vector<Object3D*> getLeafs(){ return vector<Object3D*>{&m_w.obj(m_lid[0]), &m_w.obj(m_lid[1]), &m_w.obj(m_lid[2])}; }
    Object3D& getAorta(){ return m_w.obj(m_aid); }

    void _run_computations(GuiBackInterconnection* ic = &g_interconnection);

    void _make_postsavings();

    void OccludeAorticValve(int argc, const char* argv[], const AorticSurface& as, const SewLines& sw, GuiBackInterconnection* ic = &g_interconnection);

    //cr is input-out variable,
    //  on input cr.cd[i].N should have nonnegative value and memory for cr.cd[i].bottom and cr.cd[i].top should be allocated
    void fillCoaptationMessures(CoaptationResult& cr);
    /*here will be getters of coaptation measures*/

    struct SewLinesMeta{
        typedef std::vector<Point>                                    Polyline_type;
        static void save_polyline_csv(Polyline_type& line, string fname);;

        struct SewLineMeassures{
            double len = -1; //length of sew line
            double linDist = -1;
            std::array<Point, 2> comissPnts;
            double maxDist1 = -1; //max distance between halfs of sewline according commissure plane
            std::array<Point, 2> pnts1;
            double maxDist2 = -1; //max distance between halfs of sewline according two commissure and lowest point
            std::array<Point, 2> pnts2;
            Point lowPoint; double dist_low_plane = -1; //lowest point and its distance to commissure plane
            Polyline_type ozakiLine; double ozakiLineLen = -1;
            std::string to_str(){
                std::stringstream ss;
                ss << "len=" << len << ", linD=" << linDist << ", ComissPointsCrd[2]{ " << comissPnts[0] << ", " << comissPnts[1] << " }, "
                        << "maxDist1=" << maxDist1 << ", MaxDist1Pnts[2]{ " << pnts1[0] << ", " << pnts1[1] << " }, "
                        << "maxDist2=" << maxDist2 << ", MaxDist2Pnts[2]{ " << pnts2[0] << ", " << pnts2[1] << " }, "
                        << "distLow="  << dist_low_plane << ", LowestPnt{" << lowPoint << "}";
                return ss.str();
            }
        };
        std::array<Point, 3> commissure_points;
        Vector orientation; //vector orthogonal to the plane of commissure points
        Point circle_center; double circle_R; // center of circle circumscribed around a triangle of three commissure points
        std::array<SewLineMeassures, 3> meassures;

        static SewLinesMeta Init(AVSimulator& m);
        bool prepareOzakiLine(const Object3D& aorta);
    };


    void transformLeafsIntoAorta(const vector<Object3D*>& leafs, std::array<Vector, 3> orients);//common function
    void transformLeafsIntoAorta(const Object3D& aorta, const vector<Object3D*>& leafs);
    void sewLeafBndToSewBnd(const vector<Object3D*>& leafs);
    void sewLeafletsCoarse(const Object3D& aorta, const vector<Object3D*>& leafs);
    void sewLeafletsMid(const Object3D& aorta, const vector<Object3D*>& leafs);
};

static int Listener(int argc, char* argv[]){
    InterConnector ic;
    /*do something with ic*/

#ifdef USE_MAGNUM_GUI
    GuiApplication app({argc, argv});
    return app.exec();
#else
    return 0;
#endif

    return 0;
}
//
//static int AorticValveTest(int argc, char* argv[]){
//    using namespace World3d;
//    Object3D aorticSurface("../data/Shwets/aorta.stl");
//    AorticSurface as = fromMesh(aorticSurface.m_mesh);
//    SewLines sl = readSewLines("../data/Shwets/aorta_rec_bnd.bnd");
//    CoaptationResult cr = createCoaptationResult(1000);
//    int new_argc = 3;
//    const char* new_argv[] = {"AorticValveTest", "-c", "../data/configuration.json"};
//
//    AVSimulator model;
//    std::thread front_end(Listener, argc, argv);
//    model.OccludeAorticValve(new_argc, new_argv, as, sl, &g_interconnection);
//    model.fillCoaptationMessures(cr);
//    front_end.join();
//
//    freeCoaptationResult(&cr);
//    freeSewLines(&sl);
//    freeAorticSurface(&as);
//    return 0;
//}
//int AorticValvePigsTest(int argc, char* argv[]);
//int MeasureTest();


#endif //AORTIC_VALVE_AVSIMULATOR_H
