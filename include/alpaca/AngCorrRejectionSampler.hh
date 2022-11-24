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

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/SphereRejectionSampler.hh"

namespace alpaca {

/**
 * \brief Sample directions in spherical coordinates from an angular
 * correlation.
 *
 * Compared to its base class, SphereRejectionSampler, this class provides a
 * member variable to store an AngularCorrelation object that acts as the
 * probability distribution \f$W\f$. Although SphereRejectionSampler already
 * accepts an arbitrary function of \f$\theta\f$ and \f$\varphi\f$, a function
 * object of the class AngularCorrelation cannot be passed this way. Therefore,
 * the present class was derived. For more information, see the base class.
 */
template <typename T> class AngCorrRejectionSampler {

public:
  /**
   * \brief Constructor
   *
   * In contrast to the base class, the upper limit \f$W_\mathrm{max}\f$ does
   * not have to be provided explicitly. The member function
   * AngularCorrelation::get_upper_limit() is called instead.
   *
   * \param w \f$W \left( \theta, \varphi \right)\f$, angular correlation
   * \param seed Random number seed.
   * \param max_tri Maximum number of sampled points
   * \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$
   * before the algorithm terminates without success and returns \f$\left( 0, 0
   * \right)\f$.
   */
  AngCorrRejectionSampler(AngularCorrelation<T> &w, const int seed,
                          const unsigned int max_tri = 1000)
      : angular_correlation(w.get_initial_state(), w.get_cascade_steps()),
        sampler(angular_correlation, w.get_upper_limit(), seed, max_tri) {}

  /**
   * \brief Sample random vector from a probability distribution and record the
   * number of tries.
   *
   * \return std::pair which contains \f$N\f$, the number of tries that were
   * needed to find a valid vector, and the accepted vector \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand}\right)\f$. Returns a std::pair
   * of \f$N_\mathrm{max}\f$ and \f$\left(0, 0 \right)\f$, if the maximum number
   * of trials \f$N_\mathrm{max}\f$ is reached by the algorithm and no random
   * vector was accepted.
   */
  inline pair<unsigned int, array<double, 2>> sample() {
    return sampler.sample();
  }

  inline array<double, 2> operator()() { return sampler(); }

  /**
   * \brief Sample a random vector from an arbitrarily rotated probability
   * distribution.
   *
   * This function allows to rotate the probability distribution using Euler
   * angles in the x convention. In the present code, this is an important
   * feature in gamma-ray cascades, where the direction of emission/propagation
   * and the polarization axis of the initial photon define the reference frame
   * for the emission of a subsequent photon. [Usually the analytical
   * expressions for these angular correlations are implemented in a fixed
   * reference frame (here, the direction of emission/propagation is along the
   * positive \f$z\f$ axis, and the polarization along the \f$x\f$ axis) for
   * simplicity.]
   *
   * \param euler_angles Euler angles \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$
   * in radians which define an arbitrary rotation in 3D space in the x
   * convention.
   *
   * \return Accepted vector \f$\left( \theta_\mathrm{rand},
   * \varphi_\mathrm{rand}\right)\f$. Returns \f$\left( 0, 0 \right)\f$ if the
   * maximum number of trials \f$N_\mathrm{max}\f$ is reached by the algorithm
   * and no random vector was accepted.
   */
  inline array<double, 2> operator()(const array<double, 3> euler_angles) {
    return sampler(euler_angles);
  }

protected:
  AngularCorrelation<T>
      angular_correlation; /**< Gamma-gamma angular correlation */
  SphereRejectionSampler<T, AngularCorrelation<T>>
      sampler; /**< Spherical sampler */
};

} // namespace alpaca
