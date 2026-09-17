#pragma once

#include "matrix.h"
#include <random>
#include <set>

spMtx<int> generate_mask(const size_t n, const size_t min_deg, const size_t max_deg) {
    std::set<int> positions;
    std::uniform_int_distribution<int> deg_distr(min_deg, max_deg);
    std::uniform_int_distribution<int> col_distr(0, n - 1);
    std::mt19937 generator{std::random_device{}()};
    spMtx<int> mask;
    mask.Rst = new int[n+1];
    mask.m = mask.n = n;

    mask.Rst[0] = 0;
    for (size_t i = 0; i < n; ++i)
        mask.Rst[i+1] = mask.Rst[i] + deg_distr(generator);
    mask.nz = mask.Rst[n];
    mask.Col = new int[mask.nz];

    size_t j = 0;
    for (size_t i = 0; i < mask.n; ++i) {
        size_t deg = mask.Rst[i+1] - mask.Rst[i];
        while (positions.size() < deg)
            positions.insert(col_distr(generator));
        for (const int column : positions)
            mask.Col[j++] = column;
        positions.clear();
    }

    return mask;
}


spMtx<int> generate_adjacency_matrix(const size_t n, const size_t min_deg, const size_t max_deg) {
    spMtx<int> Res = generate_mask(n, min_deg, max_deg);

    Res.Val = new int[Res.nz];
    for (size_t j = 0; j < Res.nz; ++j)
        Res.Val[j] = 1;

    return build_symm_from_lower(extract_lower_triangle(Res));
}

template <typename T>
spMtx<T> transpose(const spMtx<T> &A) {
    spMtx<T> AT(A.n, A.m, A.nz);

    // filling the column indices array and current column positions array
    for (size_t i = 0; i < A.nz; ++i)
        ++AT.Rst[A.Col[i]+1];
    for (size_t i = 0; i < AT.m; ++i)
        AT.Rst[i+1] += AT.Rst[i];

    // transposing
    for (size_t i = 0; i < A.m; ++i) {
        for (int j = A.Rst[i]; j < A.Rst[i+1]; ++j) {
            AT.Val[AT.Rst[A.Col[j]]] = std::move(A.Val[j]);
            AT.Col[AT.Rst[A.Col[j]]++] = i;
        }
    }
    // set Rst indices to normal state
    // AT.Rst[AT.m] already has the correct value
    for (int i = AT.m - 1; i > 0; --i)
        AT.Rst[i] = AT.Rst[i-1];
    AT.Rst[0] = 0;

    return AT;
}