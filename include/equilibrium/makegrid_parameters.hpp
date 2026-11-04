#pragma once

namespace equilibrium {
struct MakegridParametersSpec {
    bool normalize_by_currents = false;
    bool assume_stellarator_symmetry = false;
    int number_of_field_periods = 0;
    double r_grid_minimum = 0.0;
    double r_grid_maximum = 0.0;
    int number_of_r_grid_points = 0;
    double z_grid_minimum = 0.0;
    double z_grid_maximum = 0.0;
    int number_of_z_grid_points = 0;
    int number_of_phi_grid_points = 0;
};

inline bool operator==(const MakegridParametersSpec& a,
                       const MakegridParametersSpec& b) {
    return a.normalize_by_currents == b.normalize_by_currents &&
           a.assume_stellarator_symmetry == b.assume_stellarator_symmetry &&
           a.number_of_field_periods == b.number_of_field_periods &&
           a.r_grid_minimum == b.r_grid_minimum &&
           a.r_grid_maximum == b.r_grid_maximum &&
           a.number_of_r_grid_points == b.number_of_r_grid_points &&
           a.z_grid_minimum == b.z_grid_minimum &&
           a.z_grid_maximum == b.z_grid_maximum &&
           a.number_of_z_grid_points == b.number_of_z_grid_points &&
           a.number_of_phi_grid_points == b.number_of_phi_grid_points;
}

}  // namespace equilibrium
