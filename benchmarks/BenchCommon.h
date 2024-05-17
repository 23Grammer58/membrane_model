//
// Created by Liogky Alexey on 20.05.2022.
//

#ifndef AVSYM_BENCHCOMMON_H
#define AVSYM_BENCHCOMMON_H

#include <iostream>
#include <thread>
#include "AVSim/Core/Renderers/GuiApplication.h"
#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/DomainCollection/AniMeshCollection.h"
#include "AVSim/Core/Renderers/GuiBackInterconnection.h"
#include "AVSim/Core/World.h"
#include "AVSim/Core/ForceAppliers/ForceAppliers.h"
#include "AVSim/Core/Forces/Forces.h"
#include "AVSim/Core/Renderers/Renderer.h"
#include <Eigen/Dense>

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

using namespace std;

static int start_debug_gui(int argc, char** argv){
#ifdef USE_MAGNUM_GUI
    GuiApplication app({argc, argv});
    return app.exec();
#else
    return 0;
#endif
}

#endif //AVSYM_BENCHCOMMON_H
