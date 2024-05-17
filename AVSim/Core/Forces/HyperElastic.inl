//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_HYPERELASTIC_INL
#define AORTIC_VALVE_HYPERELASTIC_INL

namespace World3d {
    namespace HyperElasticHelpers {
        template<class Nameable>
        auto VarContainer<Nameable>::find(string name) {
            auto it = streq.find(name);
            if (it != streq.end())
                return seq.begin() + it->second;
            else
                return seq.end();
        }

        template<class Nameable>
        const auto VarContainer<Nameable>::find(string name) const {
            auto it = streq.find(name);
            if (it != streq.end())
                return seq.begin() + it->second;
            else
                return seq.end();
        }

        template<class Nameable>
        Nameable &VarContainer<Nameable>::at(string name) {
            auto it = streq.find(name);
            if (it != streq.end())
                return seq[it->second];
            else
                return seq.at(seq.size());
        }

        template<class Nameable>
        const Nameable &VarContainer<Nameable>::at(string name) const{
            auto it = streq.find(name);
            if (it != streq.end())
                return seq[it->second];
            else
                return seq.at(seq.size());
        }

        template<class Nameable>
        void VarContainer<Nameable>::push_back(Nameable var) {
            seq.push_back(var);
            update();
        }

        template<class Nameable>
        void VarContainer<Nameable>::push_back(std::initializer_list <Nameable> vars) {
            for (auto &var: vars) {
                seq.push_back(var);
            }
            update();
        }

        template<class Nameable>
        void VarContainer<Nameable>::remove(string name) {
            int id = streq[name] - seq.data();
            streq.erase(name);
            seq.erase(seq.begin() + id);
            update();
        }

        template<class Nameable>
        void VarContainer<Nameable>::remove(int id) {
            streq.erase(seq[id].getName());
            seq.erase(seq.begin() + id);
            update();
        }

        template<class Nameable>
        void VarContainer<Nameable>::update() {
            for (int i = 0; i < seq.size(); ++i) {
                streq[seq[i].getName()] = i;
            }
        }

        template<class Nameable>
        void VarContainer<Nameable>::test() {
            for (auto &i: streq) {
                std::cout << i.first << " <-> " << seq[i.second].getName() << " < - > " << seq[i.second].sym
                          << std::endl;
            }
        }
    }
}

#endif //AORTIC_VALVE_HYPERELASTIC_INL
