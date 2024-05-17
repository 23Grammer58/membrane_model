//
// Created by alex on 30.08.2021.
//

#include "MeshlabPickPoint.h"
#include "inmost.h"
using namespace INMOST;

std::ostream& MeshlabPickPoint::save(std::ostream& out) const {
    auto& ppd = *this;
    out << "<!DOCTYPE PickedPoints>\n";
    out << "<PickedPoints>\n";
    out << " <DocumentData>\n"
           "  <DateTime date=\"" << ppd.meta.data << "\" time=\"" << ppd.meta.time << "\"/>\n"
           "  <User name=\"" << ppd.meta.user << "\"/>\n"
           "  <DataFileName name=\"" << ppd.meta.dataFName << "\"/>\n"
           "  <templateName name=\"" << ppd.meta.templateName << "\"/>\n"
           " </DocumentData>\n";
    for (auto it = ppd.pp.begin(); it != ppd.pp.end(); ++it) {
        auto& p = *it;
        out << std::setprecision(8) << "<point x=\"" << p[0] << "\" y=\"" << p[1] << "\" z=\"" << p[2] << "\" name=\""
            << std::distance(ppd.pp.begin(), it) << "\" active=\"1\"/>\n";
    }
    return out << "</PickedPoints>";
}

MeshlabPickPoint &MeshlabPickPoint::read(std::string fname) {
    std::ifstream f(fname);
    if (!f.is_open()) throw std::runtime_error("Can't open file \"" + fname + "\"");
    std::string type;
    std::getline(f, type);
    if (type != "<!DOCTYPE PickedPoints>")
        std::cout << R"(WARNING: Wrong type: expected "<!DOCTYPE PickedPoints>" instead of ")" + type + "\" at row:1" << std::endl;
    XMLReader reader(fname, f);
    auto tree = reader.ReadXML();
    f.close();

    MeshlabPickPoint& ppd = *this;
    if (tree.GetName() == "PickedPoints"){
        for (auto& c: tree.children){
            if (c.GetName() == "point"){
                ppd.pp.push_back({stod(c.GetAttrib("x")), stod(c.GetAttrib("y")), stod(c.GetAttrib("z"))});
            } else if (c.GetName() == "DocumentData"){
                for (auto& dd: c.children){
                    if (dd.GetName() == "DateTime"){
                        ppd.meta.data = dd.GetAttrib("date");
                        ppd.meta.time = dd.GetAttrib("time");
                    } else if (dd.GetName() == "User")
                        ppd.meta.user = dd.GetAttrib("name");
                    else if (dd.GetName() == "DataFileName")
                        ppd.meta.dataFName = dd.GetAttrib("name");
                    else if (dd.GetName() == "templateName")
                        ppd.meta.templateName = dd.GetAttrib("name");
                    else
                        std::cout << "WARNING: Encountered unknown tag \"" << c.GetName() << "\"" << std::endl;
                }
            } else
                std::cout << "WARNING: Encountered unknown tag \"" << c.GetName() << "\"" << std::endl;
        }
    }

    return ppd;
}

std::ostream &operator<<(std::ostream &out, const MeshlabPickPoint &ppd) { return ppd.save(out); }

MeshlabPickPoint read_pp_format(std::string fname) {
    MeshlabPickPoint ppd;
    ppd.read(std::move(fname));
    return ppd;
}
