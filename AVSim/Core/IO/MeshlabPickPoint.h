//
// Created by alex on 30.08.2021.
//

#ifndef AORTIC_VALVE_MESHLABPICKPOINT_H
#define AORTIC_VALVE_MESHLABPICKPOINT_H

#include <string>
#include <vector>
#include <array>
#include <iostream>

struct MeshlabPickPoint {
    struct DocumentData{
        std::string data, time;
        std::string user;
        std::string dataFName;
        std::string templateName;
    };
    DocumentData meta;
    std::vector<std::array<double, 3>> pp;

    std::ostream& save(std::ostream& out) const;
    MeshlabPickPoint& read(std::string fname);
};

std::ostream& operator<<(std::ostream& out, const MeshlabPickPoint& ppd);
MeshlabPickPoint read_pp_format(std::string fname);


#endif //AORTIC_VALVE_MESHLABPICKPOINT_H
