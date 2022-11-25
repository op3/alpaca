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

#include <array>
#include <numbers>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/TestUtilities.hh"

using namespace alpaca;

/**
 * \brief Test rotations of 3D vectors using Euler angles.
 *
 * Test by rotating the three canonical Cartesian axes into each other.
 */
int main() {

  const double epsilon = 1e-8;

  EulerAngleRotation<double> eul_ang_rot;

  EulerAngles<double> x_axis{1., 0., 0.};
  CoordDir<double> x_axis_sph{0.5 * std::numbers::pi, 0.};
  EulerAngles<double> y_axis{0., 1., 0.};
  CoordDir<double> y_axis_sph{0.5 * std::numbers::pi, 0.5 * std::numbers::pi};
  EulerAngles<double> z_axis{0., 0., 1.};
  // Warning: On the z axis, the angle phi in spherical coordinates is actually
  // undefined. Therefore, a test in which, after a rotation into the z axis,
  // Cartesian coordinates are converted back into spherical coordinates, may
  // not result in the coordinates theta = 0, phi =0. This is taken into account
  // in the tests by only requiring that the value of theta is (numerically
  // close to) zero.
  CoordDir<double> z_axis_sph{0., 0.};

  EulerAngles<double> euler_angles{0., 0., 0.};
  EulerAngles<double> xp_yp_zp{0., 0., 0.};
  CoordDir<double> thetap_phip{0., 0.};

  // Rotate x axis into y axis
  // Phi   = -pi/2
  // Theta = 0
  // Psi   = 0
  euler_angles = EulerAngles<double>{-0.5 * std::numbers::pi, 0., 0.};

  xp_yp_zp = eul_ang_rot.rotate(x_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), y_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), x_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(x_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), y_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      2, eul_ang_rot.rotate_back(thetap_phip, euler_angles).data(),
      x_axis_sph.data(), epsilon);

  // Rotate x axis into z axis
  // Phi   = pi/2
  // Theta = pi/2
  // Psi   = 0
  euler_angles =
      EulerAngles<double>{0.5 * std::numbers::pi, 0.5 * std::numbers::pi, 0.};

  xp_yp_zp = eul_ang_rot.rotate(x_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), x_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(x_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), z_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      2, eul_ang_rot.rotate_back(thetap_phip, euler_angles).data(),
      x_axis_sph.data(), epsilon);

  // Rotate y axis into x axis
  // Phi   = pi/2
  // Theta = 0
  // Psi   = 0
  euler_angles = EulerAngles<double>{0.5 * std::numbers::pi, 0., 0.};

  xp_yp_zp = eul_ang_rot.rotate(y_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), x_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), y_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(y_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), x_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      2, eul_ang_rot.rotate_back(thetap_phip, euler_angles).data(),
      y_axis_sph.data(), epsilon);

  // Rotate y axis into z axis
  // Phi   = 0
  // Theta = -pi/2
  // Psi   = 0
  euler_angles = EulerAngles<double>{0., -0.5 * std::numbers::pi, 0.};

  xp_yp_zp = eul_ang_rot.rotate(y_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), y_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(y_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), z_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      2, eul_ang_rot.rotate_back(thetap_phip, euler_angles).data(),
      y_axis_sph.data(), epsilon);

  // From here on, note the special role of the z axis in spherical coordinates
  // which leads to arbitrary values for phi.

  // Rotate z axis into x axis
  // Phi   = 0
  // Theta = pi/2
  // Psi   = pi/2
  euler_angles =
      EulerAngles<double>{0., 0.5 * std::numbers::pi, 0.5 * std::numbers::pi};

  xp_yp_zp = eul_ang_rot.rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), x_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), z_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(z_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), x_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      eul_ang_rot.rotate_back(thetap_phip, euler_angles)[0], z_axis_sph[0],
      epsilon);

  // Rotate z axis into y axis
  // Phi   = 0
  // Theta = pi/2
  // Psi   = 0
  euler_angles = EulerAngles<double>{0., 0.5 * std::numbers::pi, 0.};

  xp_yp_zp = eul_ang_rot.rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), y_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), z_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(z_axis_sph, euler_angles);
  test_numerical_equality<double>(2, thetap_phip.data(), y_axis_sph.data(),
                                  epsilon);
  test_numerical_equality<double>(
      eul_ang_rot.rotate_back(thetap_phip, euler_angles)[0], z_axis_sph[0],
      epsilon);

  // Rotate z axis into z axis (trivial)
  // Phi   = 0
  // Theta = 0
  // Psi   = 0
  euler_angles = EulerAngles<double>{0., 0., 0.};

  xp_yp_zp = eul_ang_rot.rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, eul_ang_rot.rotate_back(xp_yp_zp, euler_angles).data(), z_axis.data(),
      epsilon);

  thetap_phip = eul_ang_rot.rotate(z_axis_sph, euler_angles);
  test_numerical_equality<double>(thetap_phip[0], z_axis_sph[0], epsilon);
  test_numerical_equality<double>(
      eul_ang_rot.rotate_back(thetap_phip, euler_angles)[0], z_axis_sph[0],
      epsilon);

  // Test the get_theta_phi method which calculates the corresponding spherical
  // coordinates theta and phi for a given normalized Cartesian vector.
  CoordDir<double> theta_phi;

  // Test the x, y, z, -x, -y, and -z axis.
  theta_phi = eul_ang_rot.get_theta_phi({1., 0., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({-1., 0., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], std::numbers::pi, epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({0., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 0.5 * std::numbers::pi,
                                  epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({0., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 3. * 0.5 * std::numbers::pi,
                                  epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({0., 0., 1.});
  test_numerical_equality<double>(theta_phi[0], 0., epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({0., 0., -1.});
  test_numerical_equality<double>(theta_phi[0], std::numbers::pi, epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  // Test more vectors in the xy plane to see whether phi is set in the correct
  // quadrant.
  theta_phi = eul_ang_rot.get_theta_phi({1., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 0.25 * std::numbers::pi,
                                  epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({-1., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 3. * 0.25 * std::numbers::pi,
                                  epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({-1., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 5. * 0.25 * std::numbers::pi,
                                  epsilon);

  theta_phi = eul_ang_rot.get_theta_phi({1., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], 0.5 * std::numbers::pi,
                                  epsilon);
  test_numerical_equality<double>(theta_phi[1], 7. * 0.25 * std::numbers::pi,
                                  epsilon);
}
