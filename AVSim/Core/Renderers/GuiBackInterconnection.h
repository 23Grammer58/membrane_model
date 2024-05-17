//
// Created by alex on 23.06.2020.
//

#ifndef AORTIC_VALVE_GUIBACKINTERCONNECTION_H
#define AORTIC_VALVE_GUIBACKINTERCONNECTION_H

#include <vector>
#include <map>
#include <string>
#include <mutex>
#include <iostream>

using namespace std;


struct RenderedObject{
    string name;

    pair<vector<array<float, 3>>, bool> vertex; //pair<array, array_is_changed>
    pair<vector<float>, bool> value;
    pair<vector<std::uint32_t>, bool> face;
    pair<vector<pair<std::uint32_t, std::uint32_t>>, bool> vertex_lbl; //pair<pair<vertex_index, vertex_label>, array_is_changed>
};

enum class PartInterconnect{
    BackEnd,
    Gui
};

enum class InterconnectDataState {
    NotModified = 0,
    Modified,
    New,
    Deleted
};

class GuiBackInterconnection {
    mutable mutex _objs_lk;
    map<int, pair<RenderedObject, InterconnectDataState>> _objs;

    mutable mutex _gui_to_back_msgs_lk;
    vector<string> _gui_to_back_msgs;
    mutable mutex _back_to_gui_msgs_lk;
    vector<string> _back_to_gui_msgs;
public:
    void send_msg(PartInterconnect com, string msg);
    std::vector<string> get_msg(PartInterconnect com);

    void read_rendered_object_data(PartInterconnect com, map<int, pair<RenderedObject, InterconnectDataState>>& objs);
    void set_rendered_object_data(PartInterconnect com, map<int, pair<RenderedObject, InterconnectDataState>>& set_objs);
    void lock_rendered_objs(PartInterconnect com);
    map<int, pair<RenderedObject, InterconnectDataState>>& get_for_modify(PartInterconnect com);
    void unlock_rendered_objs(PartInterconnect com);
};

extern GuiBackInterconnection g_interconnection;


#endif //AORTIC_VALVE_GUIBACKINTERCONNECTION_H
