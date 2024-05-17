//
// Created by alex on 30.09.2023.
//

#ifndef AORTIC_VALVE_MULTITREE_H
#define AORTIC_VALVE_MULTITREE_H

#include <vector>
#include <stack>
#include <set>
#include <map>
#include <numeric>

namespace World3d{
template<typename T>
struct MultiTree{
    using NodeID = unsigned;
    struct Node{
        T m_val = T();
        std::vector<NodeID> m_parents;
        std::vector<NodeID> m_childs;
        Node() = default;
        Node(T val): m_val{std::move(val)} {}
        Node(T val, NodeID parent) : m_val{std::move(val)} { m_parents.push_back(parent); }
        
        bool contain_parent(NodeID parent) const { return std::find(m_parents.begin(), m_parents.end(), parent) != m_parents.end(); }
        bool contain_child(NodeID child) const { return std::find(m_childs.begin(), m_childs.end(), child) != m_childs.end(); }
        std::size_t number_of_parents() const { return m_parents.size(); }
        std::size_t number_of_childs() const { return m_childs.size(); }
        const std::vector<NodeID>& parents() const { return m_parents; }
        const std::vector<NodeID>& childs() const { return m_childs; }
        const T& value() const { return m_val; }
        std::vector<NodeID>& parents(){ return m_parents; }
        std::vector<NodeID>& childs(){ return m_childs; }
        NodeID parent(unsigned i) const { return m_parents[i]; }
        NodeID child(unsigned i) const { return m_childs[i]; }
        T& value(){ return m_val; }
        bool have_parent() const { return !m_parents.empty(); }
        bool have_child() const { return !m_childs.empty(); }
        Node& add_parent(NodeID parent) { if (!contain_parent(parent)) m_parents.push_back(parent); return *this; }
        Node& add_child(NodeID child) { if (!contain_child(child)) m_childs.push_back(child); return *this; }
        Node& del_parent(NodeID parent) { auto it = std::find(m_parents.begin(), m_parents.end(), parent); if (it != m_parents.end()) m_parents.erase(it); return *this; }
        Node& del_child(NodeID child) { auto it = std::find(m_childs.begin(), m_childs.end(), child); if (it != m_childs.end()) m_childs.erase(it); return *this; }
        void setEmpty(){ m_parents.resize(0), m_childs.resize(0); m_val = T(); }
        void clear(){ m_parents.clear(), m_childs.clear(); m_val = T(); }
    };

    struct SubTreeIterator{
    protected: 
        MultiTree* back_ptr;
        NodeID start;
        unsigned niter = 0;
        std::vector<bool> discovered;
        std::stack<NodeID> to_visit;

        SubTreeIterator(MultiTree* back_ptr, NodeID start, bool set_end): start{start}, back_ptr{back_ptr} { 
            if (!set_end) {
                discovered.resize(back_ptr->m_nodes.size(), false);
                discovered[start] = true;
                to_visit.push(start);  
            }  else
                niter = std::numeric_limits<unsigned>::max();
        }
    public:
        typedef std::forward_iterator_tag  iterator_category;
        typedef Node value_type;
        typedef Node& reference;
        typedef Node* pointer;
        bool operator ==(const SubTreeIterator & other) { 
            return start == other.start && ((to_visit.empty() && other.to_visit.empty()) || niter == other.niter); 
        }
        bool operator !=(const SubTreeIterator & other) { 
            if (start != other.start || to_visit.empty() != other.to_visit.empty()) return true;
            if (to_visit.empty() && other.to_visit.empty()) return false;
            return niter != other.niter;
        }
        bool operator < (const SubTreeIterator & other) const {
            if (!to_visit.empty() && other.to_visit.empty()) return true;
            if (to_visit.empty() && !other.to_visit.empty()) return false;
            return std::accumulate(discovered.begin(), discovered.end(), 0) < std::accumulate(other.discovered.begin(), other.discovered.end(), 0); 
        }
        NodeID getID() const { return to_visit.top(); }
        reference operator*() const { return back_ptr->operator[](to_visit.top()); }
        pointer operator->() const { return &(back_ptr->operator[](to_visit.top())); }
        SubTreeIterator& operator ++() { 
            auto id = to_visit.top(); to_visit.pop();
            for (std::size_t vi = back_ptr->operator[](id).childs().size(); vi > 0; --vi){
                auto v = back_ptr->operator[](id).childs()[vi - 1];
                if (!discovered[v]) {
                    discovered[v] = true;
                    to_visit.push(v);
                }
            }
            ++niter;
            return *this;
        }
        friend class MultiTree<T>;
    };

    std::size_t size() const { return m_nodes.size() - m_removed.size(); }
    std::size_t maxID() const { return m_nodes.size(); }
    MultiTree() = default;
    MultiTree(T head_value){ m_nodes.emplace_back(); }
    NodeID createFreeNode();
    NodeID createChild(NodeID id);
    MultiTree& addLink(NodeID parent, NodeID child);
    MultiTree& removeLink(NodeID parent, NodeID child);
    bool isHaveCycles(NodeID start) const;
    MultiTree& removeNode(NodeID id);
    MultiTree& removeBranch(NodeID id);
    template<typename iterator>
    MultiTree& removeUnreacheableNodes(iterator head_node_id_begin, iterator head_node_id_end);
    MultiTree& removeUnreacheableNodes(std::initializer_list<NodeID> head_node_ids) { return removeUnreacheableNodes(head_node_ids.begin(), head_node_ids.end()); }
    MultiTree& removeUnreacheableNodes(NodeID head_id) { return removeUnreacheableNodes(&head_id, (&head_id) + 1); }
    void clear() { m_nodes.clear(); m_removed.clear(); }
    SubTreeIterator begin(NodeID start) { return SubTreeIterator(this, start, false); }
    SubTreeIterator end(NodeID start) { return SubTreeIterator(this, start, true); }

    const Node& operator[](NodeID id) const { return m_nodes[id]; }
    Node& operator[](NodeID id) { return m_nodes[id]; }
    NodeID getID(const Node& n) const { assert(&n >= m_nodes.data() && &n < m_nodes.data() + m_nodes.size() && "Wrong node"); return &n - m_nodes.data(); } 
    bool isValid(NodeID id) const { return id < m_nodes.size() && std::find(m_removed.begin(), m_removed.end(), id) == m_removed.end(); }

    std::vector<Node> m_nodes;
    std::vector<NodeID> m_removed;
};

}

#include "MultiTree.inl"

#endif //AORTIC_VALVE_MULTITREE_H