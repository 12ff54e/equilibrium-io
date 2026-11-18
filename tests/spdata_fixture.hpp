#pragma once
#include "equilibrium/spdata.hpp"

#include <array>

inline equilibrium::Spdata spdata_fixture(std::size_t lst) {
    equilibrium::Spdata data;
    data.meta_data = "spdata coefficient fixture; units supplied by caller";
    data.lsp = 3;
    data.lst = lst;
    data.psi_wall = 0.25;
    data.psi_sep = 0.5;
    std::array fields{&data.b, &data.r, &data.z, &data.j};
    for (std::size_t f = 0; f < fields.size(); ++f) {
        auto& values = *fields[f];
        values.resize(data.lsp * 9 * lst);
        for (std::size_t r = 0; r < data.lsp; ++r)
            for (std::size_t c = 0; c < 9; ++c)
                for (std::size_t t = 0; t < lst; ++t) {
                    double v = static_cast<double>(10000 * (f + 1) + 1000 * r +
                                                   100 * c + t);
                    if (r == 0 && c == 0) v = static_cast<double>(f + 1);
                    if (r == 0 && (c == 3 || c == 6)) v = 0.0;
                    values[data.field_index(r, c, t)] = v;
                }
    }
    std::array profiles{&data.q, &data.f,       &data.i,
                        &data.p, &data.r_minor, &data.psi_t};
    for (std::size_t f = 0; f < profiles.size(); ++f) {
        profiles[f]->resize(data.lsp * 3);
        for (std::size_t k = 0; k < data.lsp * 3; ++k)
            (*profiles[f])[k] = -static_cast<double>(100 * (f + 1) + k);
    }
    return data;
}
