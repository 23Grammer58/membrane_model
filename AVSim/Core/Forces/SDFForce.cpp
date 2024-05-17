#include "SDFForce.h"
using namespace World3d;

static int sdf_property_unique_n = 0;
void SDFForce::create_sdf_property(){
    if (!_obj) return;
    int n = sdf_property_unique_n;
    _sdf_data = _obj->m_mesh.add_property_map<V_ind, SignedDistanceField::SDF>("v:SDFForce::sdf_data_" + std::to_string(n)).first;
    sdf_property_unique_n++;
}
void SDFForce::clear_sdf_property(){
    if (!_obj) return;
    int n = sdf_property_unique_n;
    _obj->m_mesh.remove_property_map(_sdf_data);
}
void SDFForce::fill_sdf_data(){
    if (!_obj || !m_field) return;
    for (auto v: _obj->m_mesh.vertices()){
        if (!_obj->is_movable(v)) continue;
        _sdf_data[v] = (*m_field)(_obj->m_x[v] - CGAL::ORIGIN);
    }
}

SDFForce::SDFForce(std::shared_ptr<SignedDistanceField> field, DReal P, DReal half_penetration_depth, DReal dist_shift): m_field(std::move(field)){
    set_dist_funcs(P, half_penetration_depth, dist_shift);
    type = "SDFForce";
}
SDFForce::SDFForce(std::shared_ptr<SignedDistanceField> field, Func1D dist_f, Func1D dist_df): 
    m_field(std::move(field)), m_dist_func{dist_f}, m_dist_deriv{dist_df} {
    type = "SDFForce";
}
void SDFForce::set_dist_funcs(DReal P, DReal half_penetration_depth, DReal dist_shift){
    auto r = get_default_dist_func(P, half_penetration_depth, dist_shift);
    m_dist_func = r.first;
    m_dist_deriv = r.second;
}
std::pair<SDFForce::Func1D, SDFForce::Func1D> SDFForce::get_default_dist_func(DReal P, DReal w, DReal d_shift){
    Func1D f, df;
    f = [P, width = w, d_shift](double dist) -> double { 
        double d = dist + d_shift;
        return (2*d < width) ? P*(1 - 2*d / width) : 0; 
    };
    df = [P, width = w, d_shift](double dist) -> double { 
        double d = dist + d_shift;
        return (2*d < width) ? -2*P / width : 0; 
    };
    return {f, df};
}
std::pair<SDFForce::Func1D, SDFForce::Func1D> SDFForce::get_default_dist_func(DReal* P, DReal* w, DReal* d_shift){
    Func1D f, df;
    f = [_P = P, _w = w, _ds = d_shift](double dist) -> double { 
        DReal P = *_P, width = *_w;
        DReal d = dist;
        if (_ds) d += *_ds;
        return (2*d < width) ? P*(1 - 2*d / width) : 0; 
    };
    df = [_P = P, _w = w, _ds = d_shift](double dist) -> double { 
        DReal P = *_P, width = *_w; 
        DReal d = dist;
        if (_ds) d += *_ds;
        return (2*d < width) ? -2*P / width : 0; 
    };
    return {f, df};
}
void SDFForce::set_dist_funcs(std::function<double(double d)> dist_func, std::function<double(double d)> dist_deriv){
    m_dist_func = std::move(dist_func), m_dist_deriv = std::move(dist_deriv);
}
void SDFForce::registerObj(Object3D *obj) {
    if (_obj) clear_sdf_property();
    _obj = obj;
    _S = set_S(*obj);
    create_sdf_property();
}
int SDFForce::operator()(Object3D &obj) {
    if (&obj != _obj) registerObj(&obj);
    fill_sdf_data();
    std::vector<F_ind> ff;
    for (auto v: obj.m_mesh.vertices()){
        if (!obj.is_movable(v)) continue;
        face_around(obj.m_mesh, v, ff);
        double Aq = 0;
        for (auto f: ff) Aq += sqrt(_S[f].squared_length());
        auto sdf = _sdf_data[v];
        DReal F = m_dist_func(sdf.sdf);
        if (F != 0.0) obj.m_F[v] += obj.withBCmask(v, Aq/3 * F * sdf.grad_sdf);
    }

    return 0;
}

int SDFForce::element_matrix(Object3D &obj, ForceBase::LocMatrix &matr, World3d::F_ind element) {
    matr.resize(3 * 3, 3 * 3);
    auto v = vert_around(obj.m_mesh, element);
    auto S = _S[element];
    auto Aq = sqrt(S.squared_length());
    auto s = S / Aq;
    std::array<SignedDistanceField::SDF, 3> sdf = {_sdf_data[v[0]], _sdf_data[v[1]], _sdf_data[v[2]]};
    for (int n = 0; n < 3; ++n) {
        auto dist  =  sdf[n].sdf;
        double P = m_dist_func(dist);
        if (P == 0.0) continue;
        for (int l = 0; l < 3; ++l) {
            auto x = P / 6 * CGAL::cross_product(s, obj.m_x[v[(l + 2) % 3]] - obj.m_x[v[(l + 1) % 3]]);
            for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                matr[3*n + i][3*l+j] = sdf[n].grad_sdf[i] * x[j];
        }

        double dPddist = m_dist_deriv(dist);
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j){
            matr[3 * n + i][3 * n + j] += (sdf[n].grad_sdf[i] * sdf[n].grad_sdf[j] *dPddist)* Aq / 3;
            //matr[3 * n + i][3 * n + j] += (P*sdf[n].grad_grad_sdf(i,j))* Aq / 3;//TODO: recomment this
        }
    }

    for (auto i = 0; i < 3; ++i) {
        auto bc = obj.getBC(v[i]);
        for (int l = 0; l < 3; ++l){
            if (!(bc & (1 << l)))
                matr.ApplyDirichlet(3 * i + l);
        }
    }

    return 0;
}
