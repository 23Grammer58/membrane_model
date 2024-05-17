//
// Created by alex on 03.11.2020.
//

#include "LinearSolverInterface.h"

void SparseMatrix::LinSum(SparseMatrix& res, unsigned nmtx, const double* c/*[nmtx]*/, const SparseMatrix** A/*[nmtx]*/){
    if (nmtx == 0) return;
    if (!c || !A || !A[0]) throw std::runtime_error("Empty input coefficients or matrices");
    int N = A[0]->rows();
    if (nmtx == 1){
        if (c[0] != 0){
            res = *(A[0]);
            for (int i = 0; i < N; ++i){
                for (auto it = res[i].begin(); it != res[i].end(); ++it)
                    it->second *= c[0];
            }
        } else {
            res.resize(N);
            for (int i = 0; i < N; ++i)
                res[i].clear();
        }
        return;        
    }
    for (int k = 1; k < nmtx; ++k)
        if (A[k]->rows() != N) throw std::runtime_error("Matrix " + std::to_string(k) + " has wrong dimension");
    res.resize(N);
    std::vector<std::map<int, double>::const_iterator> Aits(nmtx);
    for (int i = 0; i < N; ++i){
        int j = N;
        for (int k = 0; k < nmtx; ++k) {
            Aits[k] = A[k]->operator[](i).begin();
            if (Aits[k] != A[k]->operator[](i).end() && Aits[k]->first < j) j = Aits[k]->first;
        }
        std::map<int, double>::iterator it = res.operator[](i).begin();
        while (it != res[i].end() && it->first < j) it = res[i].erase(it);
        while (j < N){
            double val = 0;
            int next_j = N;
            for (int k = 0; k < nmtx; ++k){
                if (Aits[k] == A[k]->operator[](i).end()) continue;
                if (Aits[k]->first == j){
                    val += c[k]*Aits[k]->second;
                    ++Aits[k];
                    if (Aits[k] != A[k]->operator[](i).end() && Aits[k]->first < next_j) next_j = Aits[k]->first;
                } else {
                    if (Aits[k]->first < next_j) next_j = Aits[k]->first;
                }
            }
            if (val != 0){
                if (it->first == j){
                    it->second = val;
                    ++it;
                } else 
                    res[i].emplace_hint(it, j, val);
            }
            
            j = next_j;
            while (it != res[i].end() && it->first < j) it = res[i].erase(it);
        }
    }
    return;
}
