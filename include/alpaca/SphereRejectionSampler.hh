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

#include <enoki/random.h>
#include <functional>
#include <numbers>
#include <numeric>
#include <utility>
#include <vector>

using std::pair;
using std::vector;

#include "alpaca/EulerAngleRotation.hh"

namespace alpaca {

using Distribution = std::function<double(const double, const double)>;

/**
 * \brief Sample from a probability distribution in spherical coordinates using
 * rejection sampling.
 *
 * Given a probability distribution
 *
 * \f[
 * W\left( \Omega \right) \mathrm{d}\Omega = W \left( \theta, \varphi\right)
 * \sin \left( \theta \right) \mathrm{d} \theta \mathrm{d} \varphi, \f]
 *
 * which determines the probability of finding a vector in the solid angle
 * element \f$\mathrm{d}\Omega\f$ (expressed as \f$\sin \left( \theta \right)
 * \mathrm{d} \theta \mathrm{d} \varphi\f$ in spherical coordinates, with the
 * polar angle \f$\theta\f$ and the aximuthal angle \f$\varphi\f$) this class
 * samples random pairs \f$\left( \theta, \varphi \right)\f$ distributed
 * according to \f$W\f$.
 *
 * The sampling algorithm used here is 'rejection sampling' (see, e.g. Sec. 2.3
 * in Ref. \cite RobertCasella1999). It requires an upper limit
 * \f$W_\mathrm{max}\f$ for the maximum value of the distribution:
 *
 * \f[
 *      W_\mathrm{max} \geq \mathrm{max}_{\theta \in \left[ 0, \pi \right],
 * \varphi \in \left[ 0, 2\pi \right]} W \left( \theta, \varphi \right). \f]
 *
 * It starts by sampling a point \f$\left( \theta_\mathrm{rand},
 * \varphi_\mathrm{rand} \right)\f$ from a uniform distribution on a sphere
 * surface (see, e.g. Ref. \cite Weisstein20202):
 *
 * \f[
 *      \theta_\mathrm{rand} = \mathrm{arccos} \left( 2u - 1 \right)
 * \f]
 * \f[
 *      \varphi_\mathrm{rand} = 2 v \pi.
 * \f]
 *
 * Here, \f$u\f$ and \f$v\f$ denote independent uniform random numbers in the
 * range \f$\left[ 0, 1 \right]\f$. After that, a uniform random number
 * \f$W_\mathrm{rand}\f$ in the range \f$\left[ 0, W_\mathrm{max} \right]\f$ is
 * sampled. If
 *
 * \f[
 *      W_\mathrm{rand} \leq W \left( \theta_\mathrm{rand},
 * \varphi_\mathrm{rand} \right), \f]
 *
 * then the vector \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand}
 * \right)\f$ is accepted. If not, the vector is rejected and the algorithm
 * start from the beginning by sampling a new point on the sphere surface. Note
 * that the last inequality ensures that the probability of sampling a vector
 * \f$\left( \theta, \varphi \right)\f$ is proportional to \f$W \left( \theta,
 * \varphi \right)\f$. This algorithm obviously works best if the probability of
 * a vector being rejected is small, i.e. if \f$W_\mathrm{max}\f$ is a good
 * approximation of the real maximum of \f$W\f$, and if the distribution is not
 * too different from a uniform distribution \f$W = \left( 4 \pi
 * \right)^{-1}\f$. Note that rejection sampling does not require the
 * distribution to be normalized.
 *
 * In principle, rejection sampling proceeds in an infinite loop until a valid
 * vector has been found. Here, the loop is limited to a maximum number of trial
 * vectors \f$N_\mathrm{max}\f$. If this maximum number is reached, the point
 * \f$\left( 0, 0 \right)\f$ is returned.
 *
 * As a measure of the efficiency of the algorithm for a given distribution,
 * define the efficiency
 *
 * \f[
 *      \epsilon = \frac{1}{\langle N \rangle},
 * \f]
 *
 * which is the inverse of the average number of trial vectors \f$\langle N
 * \rangle\f$ that have to be sampled before a vector is accepted.
 */
template <typename T, typename Dist> class SphereRejectionSampler {

public:
  // virtual ~SphereRejectionSampler() = default;
  // SphereRejectionSampler(const SphereRejectionSampler &) = default;
  // SphereRejectionSampler &operator=(const SphereRejectionSampler &) =
  // default; SphereRejectionSampler(SphereRejectionSampler &&) = default;
  // SphereRejectionSampler &operator=(SphereRejectionSampler &&) = default;

  /**
   * \brief Constructor
   *
   * \param dis \f$W \left( \theta, \varphi \right)\f$, probability distribution
   * in spherical coordinates. \param dis_max \f$W_\mathrm{max}\f$, upper limit
   * for the maximum of the probability distribution. \param seed Random number
   * seed. \param max_tri Maximum number of sampled points \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$ before the algorithm
   * terminates without success and returns \f$\left( 0, 0 \right)\f$ (default:
   * 1000).
   */
  SphereRejectionSampler(Dist dis, const double dis_max, const size_t seed,
                         const unsigned int max_tri = 1000)
      : distribution(dis), distribution_maximum(dis_max), max_tries(max_tri),
        rng(seed) {}

  /**
   * \brief Sample a random vector from probability distribution and record the
   * number of tries.
   *
   * \return std::pair which contains \f$N\f$, the number of tries that were
   * needed to find a valid vector, and the accepted vector \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand}\right)\f$. Returns a std::pair
   * of \f$N_\mathrm{max}\f$ and \f$\left(0, 0 \right)\f$, if the maximum number
   * of trials \f$N_\mathrm{max}\f$ is reached by the algorithm and no random
   * vector was accepted.
   */
  pair<unsigned int, array<double, 2>> sample() {
    array<double, 2> theta_phi;
    double dis_val;

    for (unsigned int i = 0; i < max_tries; ++i) {

      theta_phi = sample_theta_phi();
      dis_val = sample_val();

      if (dis_val <= distribution(theta_phi[0], theta_phi[1])) {
        return {i + 1, {theta_phi[0], theta_phi[1]}};
      }
    }

    return {max_tries, {0., 0.}};
  }

  array<double, 2> operator()() {
    pair<unsigned int, array<double, 2>> sampled_theta_phi = sample();

    return {sampled_theta_phi.second[0], sampled_theta_phi.second[1]};
  }

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
  array<double, 2> operator()(const array<double, 3> euler_angles) {
    return euler_angle_rotation.rotate(operator()(), euler_angles);
  }

  /**
   * \brief Estimate the efficiency of rejection sampling for the given
   * distribution.
   *
   * The efficiency \f$\epsilon\f$ is estimated by sampling \f$n\f$ vectors from
   * the distribution and calculating the average number of trials \f$\langle N
   * \rangle\f$ for this set.
   *
   * \param n_tries \f$n\f$, number of sampled vectors.
   *
   * \return Estimate for \f$\epsilon\f$ from the \f$n\f$ samples.
   */
  double estimate_efficiency(const unsigned int n_tries) {
    vector<unsigned int> required_tries(n_tries);

    pair<unsigned int, array<double, 2>> sampled_theta_phi;

    for (unsigned int i = 0; i < n_tries; ++i) {
      sampled_theta_phi = sample();
      required_tries[i] = sampled_theta_phi.first;
    }

    return static_cast<double>(n_tries) /
           static_cast<double>(std::accumulate(required_tries.begin(),
                                               required_tries.end(), 0));
  }

protected:
  /**
   * \brief Sample random value for W between [0, distribution_maximum).
   *
   * \return Random (uniform) value for W.
   */
  inline T sample_val() { return rng.next_float64() * distribution_maximum; }

  /**
   * \brief Sample polar angle of a uniformly randomly distributed point on a
   * sphere surface.
   *
   * \return \f$\theta_\mathrm{rand}\f$, random polar angle.
   */
  inline T sample_theta() { return enoki::acos(2. * rng.next_float64() - 1.); }

  /**
   * \brief Sample azimuthal angle of a uniformly randomly distributed point on
   * a sphere surface.
   *
   * \return \f$\varphi_\mathrm{rand}\f$, random azimuthal angle.
   */
  inline T sample_phi() { return 2. * std::numbers::pi * rng.next_float64(); };

  /**
   * \brief Sample uniformly randomly distributed point on a sphere surface.
   *
   * \return \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$,
   * random point on sphere surface.
   */
  inline array<T, 2> sample_theta_phi() {
    return {sample_theta(), sample_phi()};
  }

  Dist distribution;           /**< \f$W \left( \theta, \varphi \right)\f$,
                                          (unnormalized)   probability distribution. */
  double distribution_maximum; /**< \f$W_\mathrm{max}\f$, maximum of
                                        probability distribution. */
  unsigned int max_tries;      /**< \f$N_\mathrm{max}\f$, maximum number of
                                        tries to find a random vector. */

  enoki::PCG32<T> rng; /**< Deterministic random number engine. */
  EulerAngleRotation<T>
      euler_angle_rotation; /**< Instance of the EulerAngleRotation class */
};

} // namespace alpaca
