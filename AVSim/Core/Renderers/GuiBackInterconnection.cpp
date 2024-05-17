//
// Created by alex on 23.06.2020.
//

#include "GuiBackInterconnection.h"


GuiBackInterconnection g_interconnection;

void GuiBackInterconnection::send_msg(PartInterconnect com, string msg) {
    switch (com) {
        case PartInterconnect::BackEnd:{
            std::lock_guard<std::mutex> guard(_back_to_gui_msgs_lk);
            _back_to_gui_msgs.push_back(msg);
            break;
        }
        case PartInterconnect::Gui:{
            std::lock_guard<std::mutex> guard(_gui_to_back_msgs_lk);
            _gui_to_back_msgs.push_back(msg);
            break;
        }
    }
}

std::vector<string> GuiBackInterconnection::get_msg(PartInterconnect com) {
    std::vector<string> result;
    switch (com) {
        case PartInterconnect::BackEnd:{
            std::lock_guard<std::mutex> guard(_gui_to_back_msgs_lk);
            result = std::move(_gui_to_back_msgs);
            _gui_to_back_msgs.clear();
            break;
        }
        case PartInterconnect::Gui:{
            std::lock_guard<std::mutex> guard(_back_to_gui_msgs_lk);
            result = std::move(_back_to_gui_msgs);
            _back_to_gui_msgs.clear();
            break;
        }
    }
    return result;
}

void GuiBackInterconnection::read_rendered_object_data(PartInterconnect com,
                                                       map<int, pair<RenderedObject, InterconnectDataState>> &objs) {
    (void) com;
    std::lock_guard<std::mutex> guard(_objs_lk);

    for (auto it = _objs.begin(); it != _objs.end(); ){
        auto& obj = *it;
        decltype(objs.find(obj.first)) q;
        bool next = true;
        if (obj.second.second == InterconnectDataState::New || (q = objs.find(obj.first)) == objs.end()){
            auto p = objs.insert(pair<int, pair<RenderedObject, InterconnectDataState>>{obj.first, pair<RenderedObject, InterconnectDataState>()});
            std::swap(p.first->second, obj.second);
            p.first->second.second = InterconnectDataState::New;
            obj.second.first.face.second = false;
            obj.second.first.vertex.second = false;
            obj.second.first.value.second = false;
            obj.second.first.vertex_lbl.second = false;
            obj.second.second = InterconnectDataState::NotModified;
        }else if (obj.second.second == InterconnectDataState::Modified){
            auto& p = q->second;
            p.second = InterconnectDataState::Modified;
            auto& inob = p.first;
            auto& fromob = obj.second.first;
            if (fromob.vertex.second) {
                swap(fromob.vertex, inob.vertex);
                fromob.vertex.second = false;
                fromob.vertex_lbl.second = true;
            } else inob.vertex.second = false;
            if (fromob.face.second) {
                swap(fromob.face, inob.face);
                fromob.face.second = false;
            } else inob.face.second = false;
            if (fromob.value.second) {
                swap(fromob.value, inob.value);
                fromob.value.second = false;
            } else inob.value.second = false;
            if (fromob.vertex_lbl.second) {
                swap(fromob.vertex_lbl, inob.vertex_lbl);
                fromob.vertex_lbl.second = false;
            } else inob.vertex_lbl.second = false;
            obj.second.second = InterconnectDataState::NotModified;
        } else if (obj.second.second == InterconnectDataState::Deleted){
            auto& p = q->second;
            p.second = InterconnectDataState::Deleted;
            it = _objs.erase(it);
            next = false;
        }
        if (next) ++it;
    }
}

void GuiBackInterconnection::set_rendered_object_data(PartInterconnect com,
                                                      map<int, pair<RenderedObject, InterconnectDataState>> &set_objs) {
    (void) com;
    std::lock_guard<std::mutex> guard(_objs_lk);
    map<int, pair<RenderedObject, InterconnectDataState>>& mesh = _objs;
    for (auto it = set_objs.begin(); it != set_objs.end(); ) {
        RenderedObject& obj = (*it).second.first;
        bool flag = true;
        switch ((*it).second.second){
            case InterconnectDataState::New:{
                mesh.insert(*it);
                break;
            }
            case InterconnectDataState::Modified:{
                auto& mop = mesh[(*it).first];
                mop.second = InterconnectDataState::Modified;
                auto& mo = mop.first;
                if (obj.vertex.second) {
                    swap(mo.vertex, obj.vertex);
                    obj.vertex.second = false;
                }
                if (obj.value.second) {
                    swap(mo.value, obj.value);
                    obj.value.second  = false;
                }
                if (obj.face.second) {
                    swap(mo.face, obj.face);
                    obj.face.second  = false;
                }
                if (obj.vertex_lbl.second) {
                    swap(mo.vertex_lbl, obj.vertex_lbl);
                    obj.vertex_lbl.second  = false;
                }
                break;
            }
            case InterconnectDataState::NotModified:{
                break;
            }
            case InterconnectDataState::Deleted:{
                mesh[(*it).first].second = InterconnectDataState::Deleted;
                it = set_objs.erase(it);
                flag = false;
                break;
            }
        }
        if (flag) {
            (*it).second.second = InterconnectDataState::NotModified;
            ++it;
        }
    }
}

void GuiBackInterconnection::lock_rendered_objs(PartInterconnect com) {
    (void) com;
    _objs_lk.lock();
}

map<int, pair<RenderedObject, InterconnectDataState>> &GuiBackInterconnection::get_for_modify(PartInterconnect com) {
    (void) com;
    if (_objs_lk.try_lock())
        throw runtime_error("For modifying RenderedObject you should initially call lock_rendered_objs(...)");
    return _objs;
}

void GuiBackInterconnection::unlock_rendered_objs(PartInterconnect com) {
    (void) com;
    _objs_lk.unlock();
}
