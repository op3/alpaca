/*
    This file is part of alpaca.

    alpaca is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    alpaca is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

    Copyright (C) 2021 Udo Friman-Gayer
*/

#include <cmath>
#include <numbers>

#include <enoki/array.h>

#include "alpaca/EulerAngleRotation.hh"

namespace alpaca {

array<double, 2>
EulerAngleRotation::get_theta_phi(const array<double, 3> x_y_z_norm) const {

  return {acos(x_y_z_norm[2]),
          fmod(atan2(x_y_z_norm[1], x_y_z_norm[0]) + 2. * std::numbers::pi,
               2. * std::numbers::pi)};
}

array<double, 3>
EulerAngleRotation::get_x_y_z_norm(const array<double, 2> theta_phi) const {

  double cos_the{cos(theta_phi[0])}, sin_the{sin(theta_phi[0])},
      cos_phi{cos(theta_phi[1])}, sin_phi{sin(theta_phi[1])};

  return array<double, 3>{sin_the * cos_phi, sin_the * sin_phi, cos_the};
}

array<array<double, 3>, 3> EulerAngleRotation::rotation_matrix(
    const array<double, 3> Phi_Theta_Psi) const {

  auto [sin_phi, cos_phi] = enoki::sincos(Phi_Theta_Psi[0]);
  auto [sin_the, cos_the] = enoki::sincos(Phi_Theta_Psi[1]);
  auto [sin_psi, cos_psi] = enoki::sincos(Phi_Theta_Psi[2]);

  return array<array<double, 3>, 3>{
      array<double, 3>{cos_psi * cos_phi - cos_the * sin_phi * sin_psi,
                       cos_psi * sin_phi + cos_the * cos_phi * sin_psi,
                       sin_psi * sin_the},
      array<double, 3>{-sin_psi * cos_phi - cos_the * sin_phi * cos_psi,
                       -sin_psi * sin_phi + cos_the * cos_phi * cos_psi,
                       cos_psi * sin_the},
      array<double, 3>{sin_the * sin_phi, -sin_the * cos_phi, cos_the}};
}

} // namespace alpaca
