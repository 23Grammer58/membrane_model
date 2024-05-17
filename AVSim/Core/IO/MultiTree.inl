//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_MULTITREE_INL
#define AORTIC_VALVE_MULTITREE_INL

#include "MultiTree.h"

namespace World3d{

template<typename T>
typename MultiTree<T>::NodeID MultiTree<T>::createFreeNode() {
    NodeID id = std::numeric_limits<NodeID>::max();
    if (!m_removed.empty()){
        id = m_removed.back();
        m_removed.resize(m_removed.size() - 1);
    } else {
        id = m_nodes.size();
        m_nodes.resize(m_nodes.size() + 1);
    }
    return id;
}
template<typename T>
typename MultiTree<T>::NodeID MultiTree<T>::createChild(NodeID id) { 
    auto res = createFreeNode();
    m_nodes[id].add_child(res);
    m_nodes[res].add_parent(id);
    return res;
}
template<typename T>
MultiTree<T>& MultiTree<T>::addLink(NodeID parent, NodeID child){
    m_nodes[parent].add_child(child);
    m_nodes[child].add_parent(parent);
    return *this;
}
template<typename T>
MultiTree<T>& MultiTree<T>::removeLink(NodeID parent, NodeID child){
    m_nodes[parent].del_child(child);
    m_nodes[child].del_parent(parent);
    return *this;
}
template<typename T>
template<typename iterator>
MultiTree<T>& MultiTree<T>::removeUnreacheableNodes(iterator node_id_begin, iterator node_id_end){
    std::vector<bool> discovered(m_nodes.size(), false);
    std::stack<NodeID> to_visit; 
    for (auto it = node_id_begin; it != node_id_end; ++it) to_visit.push(*it);
    for (auto v: m_removed) discovered[v] = true;
    while(!to_visit.empty()){
        auto u = to_visit.top(); to_visit.pop();
        if (discovered[u]) continue;
        discovered[u] = true;
        for (auto v: m_nodes[u].childs())
            if (!discovered[v]) to_visit.push(v);
    }
    for (std::size_t v = 0; v < m_nodes.size(); ++v) 
        if (!discovered[v]){ 
            m_nodes[v].setEmpty();
            m_removed.push_back(v);
        }
    return *this;    
}

template<typename T>
bool MultiTree<T>::isHaveCycles(NodeID start) const {
    std::vector<bool> discovered(m_nodes.size(), false), finished(m_nodes.size(), false);
    std::stack<unsigned> to_visit;
    to_visit.push(start);
    while(!to_visit.empty()){
        auto u = to_visit.top(); to_visit.pop();
        if (u == std::numeric_limits<unsigned>::max()){
            u = to_visit.top(); to_visit.pop();
            finished[u] = true;
            discovered[u] = false;
            continue;
        }
        discovered[u] = true;
        to_visit.push(u); to_visit.push(std::numeric_limits<unsigned>::max());
        for (auto v: m_nodes[u].childs()){
            if (discovered[v]) return true;
            if (!finished[v]) to_visit.push(v);
        }
    }
    return false;
}
template<typename T>
MultiTree<T>& MultiTree<T>::removeNode(NodeID id){
    std::sort(m_removed.begin(), m_removed.end());
    for (unsigned i = 0, j = 0; i < m_nodes.size(); ++i){
        if (j < m_removed.size() && i == m_removed[j]) { ++j; continue; }
        if (i == id) continue;
        m_nodes[i].del_parent(id);
        m_nodes[i].del_child(id);
    }
    m_nodes[id].setEmpty();
    m_removed.push_back(id);
    return *this;
}
template<typename T>
MultiTree<T>& MultiTree<T>::removeBranch(NodeID id){
    std::set<unsigned> to_remove;
    {
        std::stack<unsigned> not_visited;
        not_visited.push(id);
        while(!not_visited.empty()){
            auto lid = not_visited.top(); not_visited.pop();
            to_remove.insert(lid);
            for (auto v: m_nodes[lid].childs())
                if (to_remove.insert(v).second) 
                    not_visited.push(v);
        }
    }
    for (unsigned i = 0; i < m_nodes.size(); ++i)
        for (auto lid: to_remove){
            m_nodes[i].del_parent(lid);
            m_nodes[i].del_child(lid);
        }

    for (auto lid: to_remove){
        m_nodes[lid].setEmpty();
        m_removed.push_back(lid);
    }
    
    return *this;
}

}

#endif //MULTITREE