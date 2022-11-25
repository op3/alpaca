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

#pragma once

#include <enoki/array.h>
#include <numbers>

namespace alpaca {

/** Euler Angles */
template <typename T> using EulerAngles = enoki::Array<T, 3>;

/** Directional coordinates (theta, phi) */
template <typename T> using CoordDir = enoki::Array<T, 2>;

/**
 * \brief Class to perform arbitrary rotations of 3D vectors using Euler angles.
 *
 * Any rotation in three-dimensional space can be expressed in terms of three
 * Euler angles, which are denoted as \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$
 * here \cite Weisstein2020. Here, an arbitrary orientation of a coordinate
 * system with the axes \f$x\f$, \f$y\f$, and \f$z\f$ into a coordinate system
 * with the new axes \f$x^\prime\f$, \f$y^\prime\f$, and \f$z^\prime\f$ is
 * achieved by the 'zxz' scheme or the 'x-convention': First, the original
 * vector \f$v\f$ in the \f$xyz\f$ system is rotated around the original \f$z\f$
 * axis by an angle \f$\Phi\f$. The rotation matrix is given by:
 *
 * \f[
 *      D\left( \Phi \right) = \begin{pmatrix}
 *           \cos \left( \Phi \right) &  \sin \left( \Phi \right) & 0 \\
 *          -\sin \left( \Phi \right) &  \cos \left( \Phi \right) & 0 \\
 *           0                        &  0                        & 1 \\
 *      \end{pmatrix}
 * \f]
 *
 * Then, the resulting vector is rotated around the original \f$x\f$ axis by an
 * angle \f$\Theta\f$:
 *
 * \f[
 *      C\left( \Theta \right) = \begin{pmatrix}
 *           1                        &  0                          & 0 \\
 *           0                        &  \cos \left( \Theta \right) &  \sin
 * \left( \Theta \right) \\
 *           0                        & -\sin \left( \Theta \right) &  \cos
 * \left( \Theta \right) \\ \end{pmatrix} \f]
 *
 * At last, the vector is rotated around the original \f$z\f$ axis by an angle
 * \f$\Psi\f$.
 *
 * \f[
 *      B\left( \Psi \right) = \begin{pmatrix}
 *           \cos \left( \Psi \right) &  \sin \left( \Psi \right) & 0 \\
 *          -\sin \left( \Psi \right) &  \cos \left( \Psi \right) & 0 \\
 *           0                        &  0                        & 1 \\
 *      \end{pmatrix}
 * \f]
 *
 * The matrices are named in the same way as Eq. (1) in Ref. \cite
 * Weisstein2020, but they have explicit angle arguments. The action of the
 * three rotation matrices on the vector \f$v\f$ results in a new vector
 * \f$v^\prime\f$:
 *
 * \f[
 *      v^\prime = \underbrace{B \left( \Psi \right) C \left( \Theta \right) D
 * \left( \Phi \right)}_{\equiv A \left( \Phi, \Theta, \Psi \right)} v. \f]
 *
 * In the equation above, the symbol \f$A\f$ for the total rotation matrix has
 * been introduced. In order to obtain \f$v\f$ from \f$v^\prime\f$, the
 * rotations can be reversed:
 *
 * \f[
 *      v = \underbrace{D \left( -\Phi \right) C \left( -\Theta \right) B \left(
 * -\Psi \right)}_{\equiv A^{-1} \left( \Phi, \Theta, \Psi \right)} v^\prime.
 * \f]
 *
 * This is equivalent to calculating the inverse matrix \f$A^{-1}\f$ of \f$A\f$.
 *
 * For the present application, it is instructive to consider the action of the
 * three matrices on a unit vector along the positive \f$z\f$ axis, since the
 * default orientation of all angular correlations is along the \f$z\f$ axis (in
 * other words: the first photon is supposed to be emitted/propagating along the
 * \f$z\f$ axis):
 *
 * \f[
 *      B \left( \Psi \right) C \left( \Theta \right) D \left( \Phi \right)
 * \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = B \left( \Psi \right) C \left(
 * \Theta \right) \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix}
 * \sin \left( \Theta \right) \sin \left( \Psi \right) \\ \sin \left( \Theta
 * \right) \cos \left( \Psi \right) \\ \cos \left( \Theta \right) \end{pmatrix}.
 * \f]
 *
 * For the convention used in this code, an arbitrary vector in spherical
 * coordinates is identified by its polar angle \f$\theta\f$ and its azimuthal
 * angle \f$\varphi\f$:
 *
 * \f[
 *      \begin{pmatrix} \sin \left( \theta \right) \cos \left( \varphi \right)
 * \\ \sin \left( \theta \right) \sin \left( \varphi \right) \\ \cos \left(
 * \theta \right) \end{pmatrix}. \f]
 *
 * By comparing the two equations above, it can be seen that a given vector in
 * spherical coordinates can be obtained by applying the matrix \f$B C\f$ to a
 * unit vector along the positive \f$z\f$ axis with the following choice of the
 * Euler angles:
 *
 * \f[
 *      \Theta = \theta
 * \f]
 * \f[
 *      \Psi = -\varphi + \frac{\pi}{2}.
 * \f]
 *
 * For the representation of 2D and 3D vectors, this class uses the std::array
 * container class.
 *
 */
template <typename T> class EulerAngleRotation {

public:
  /**
   * \brief Rotate a 3D vector.
   *
   * \param x_y_z \f$v\f$, 3D vector
   * \param phi_theta_psi Euler angles in radians
   *
   * \return \f$v^\prime\f$, 3D vector
   */
  inline EulerAngles<T> rotate(const EulerAngles<T> x_y_z,
                               const EulerAngles<T> Phi_Theta_Psi) const {
    if (no_rotation_required(Phi_Theta_Psi)) {
      return x_y_z;
    }

    const enoki::Array<EulerAngles<T>, 3> A = rotation_matrix(Phi_Theta_Psi);
    return EulerAngles<T>{
        A[0][0] * x_y_z[0] + A[0][1] * x_y_z[1] + A[0][2] * x_y_z[2],
        A[1][0] * x_y_z[0] + A[1][1] * x_y_z[1] + A[1][2] * x_y_z[2],
        A[2][0] * x_y_z[0] + A[2][1] * x_y_z[1] + A[2][2] * x_y_z[2]};
  }

  /**
   * \brief Rotate a 3D vector.
   *
   * \param theta_phi \f$v\f$, spherical coordinates \f$\theta\f$ and
   * \f$\varphi\f$ in radians. \param Phi_Theta_Psi Euler angles in radians
   *
   * \return \f$v^\prime\f$, spherical coordinates \f$\theta\f$ and
   * \f$\varphi\f$ in radians.
   */
  inline CoordDir<T> rotate(const CoordDir<T> theta_phi,
                            const EulerAngles<T> Phi_Theta_Psi) const {
    if (no_rotation_required(Phi_Theta_Psi)) {
      return theta_phi;
    }
    return get_theta_phi(rotate(get_x_y_z_norm(theta_phi), Phi_Theta_Psi));
  }

  /**
   * \brief Rotate a 3D vector back.
   *
   * Performs the same action on a 3D vector which would be performed by
   * EulerAngleRotation::rotate() if the angles \f$\Phi\f$ and \f$\Psi\f$ are
   * switched, and the negative value of each angle is used.
   *
   * \param x_y_z \f$v^\prime\f$, 3D vector
   * \param Phi_Theta_Psi Euler angles in radians
   *
   * \return \f$v\f$, 3D vector
   */
  inline EulerAngles<T> rotate_back(const EulerAngles<T> xp_yp_zp,
                                    const EulerAngles<T> Phi_Theta_Psi) const {
    if (no_rotation_required(Phi_Theta_Psi)) {
      return xp_yp_zp;
    }

    return rotate(xp_yp_zp,
                  {-Phi_Theta_Psi[2], -Phi_Theta_Psi[1], -Phi_Theta_Psi[0]});
  }

  /**
   * \brief Rotate a 3D vector back.
   *
   * See also the implementation of rotate_back for Cartesian vectors.
   *
   * \param theta_phi \f$v^\prime\f$, spherical coordinates \f$\theta\f$ and
   * \f$\varphi\f$ in radians. \param Phi_Theta_Psi Euler angles in radians
   *
   * \return \f$v\f$, spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in
   * radians.
   */
  inline CoordDir<T> rotate_back(const CoordDir<T> thetap_phip,
                                 const EulerAngles<T> Phi_Theta_Psi) const {
    if (no_rotation_required(Phi_Theta_Psi)) {
      return thetap_phip;
    }

    return get_theta_phi(
        rotate_back(get_x_y_z_norm(thetap_phip), Phi_Theta_Psi));
  }

  /**
   * \brief Convert Cartesian to spherical coordinates.
   *
   * Given a normalized Cartesian vector with three coordinates \f$x\f$,
   * \f$y\f$, and \f$z\f$, this function calculates the corresponding angles
   * \f$\theta\f$ and \f$\varphi\f$ in spherical coordinates.
   *
   * At the moment, the input vector is not tested for normalization.
   *
   * \param x_y_z_norm Normalized Cartesian vector.
   *
   * \return Spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in radians.
   */
  inline CoordDir<T> get_theta_phi(const EulerAngles<T> x_y_z_norm) const {

    return {enoki::acos(x_y_z_norm[2]),
            enoki::fmod(enoki::atan2(x_y_z_norm[1], x_y_z_norm[0]) +
                            2. * std::numbers::pi,
                        2. * std::numbers::pi)};
  }

protected:
  /**
   * \brief Convert spherical to Cartesian coordinates.
   *
   * Given spherical coordinates \f$\theta\f$ and \f$\varphi\f$, this function
   * calculates the corresponding normalized Cartesian three-component vector.
   *
   * \param theta_phi Spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in
   * radians.
   *
   * \return Normalized Cartesian vector.
   */
  inline EulerAngles<T> get_x_y_z_norm(const CoordDir<T> theta_phi) const {

    auto [sin_the, cos_the] = enoki::sincos(theta_phi[0]);
    auto [sin_phi, cos_phi] = enoki::sincos(theta_phi[1]);

    return EulerAngles<T>{sin_the * cos_phi, sin_the * sin_phi, cos_the};
  }

  /**
   * \brief Calculate rotation matrix for the three Euler angles.
   *
   * \param Phi_Theta_Psi Euler angles in radians
   *
   * \return \f$3 \times 3\f$ matrix \f$A\f$
   */
  enoki::Array<EulerAngles<T>, 3>
  rotation_matrix(const EulerAngles<T> Phi_Theta_Psi) const {

    auto [sin_phi, cos_phi] = enoki::sincos(Phi_Theta_Psi[0]);
    auto [sin_the, cos_the] = enoki::sincos(Phi_Theta_Psi[1]);
    auto [sin_psi, cos_psi] = enoki::sincos(Phi_Theta_Psi[2]);

    return enoki::Array<EulerAngles<T>, 3>{
        EulerAngles<T>{cos_psi * cos_phi - cos_the * sin_phi * sin_psi,
                       cos_psi * sin_phi + cos_the * cos_phi * sin_psi,
                       sin_psi * sin_the},
        EulerAngles<T>{-sin_psi * cos_phi - cos_the * sin_phi * cos_psi,
                       -sin_psi * sin_phi + cos_the * cos_phi * cos_psi,
                       cos_psi * sin_the},
        EulerAngles<T>{sin_the * sin_phi, -sin_the * cos_phi, cos_the}};
  }

  /**
   * \brief Check if all three Euler angles are zero.
   *
   * This function is used in the code to decide whether a calculation needs to
   * be performed or whether the input value can simply be returned.
   *
   * \param Phi_Theta_Psi Euler angles in radians
   *
   * \return true, if all Euler angles are zero, else false.
   */
  inline bool no_rotation_required(const EulerAngles<T> Phi_Theta_Psi) const {
    return (Phi_Theta_Psi[0] == 0. && Phi_Theta_Psi[1] == 0. &&
            Phi_Theta_Psi[2] == 0.);
  }
};

} // namespace alpaca
