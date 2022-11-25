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

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/State.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"
#include "alpaca/W_pol_dir.hh"

using std::invalid_argument;
using std::pair;
using std::unique_ptr;
using std::vector;

namespace alpaca {

inline vector<CascadeStep> gen_cascade_steps(State ini_sta,
                                             const vector<State> cas_sta) {
  vector<CascadeStep> cascade_steps;
  cascade_steps.reserve(cas_sta.size());

  cascade_steps.push_back({Transition{ini_sta, cas_sta[0]}, cas_sta[0]});

  for (size_t i = 1; i < cas_sta.size(); ++i) {
    cascade_steps.push_back(
        {Transition{cas_sta[i - 1], cas_sta[i]}, cas_sta[i]});
  }
  return cascade_steps;
}

/**
 * \brief Class for a gamma-gamma correlation.
 *
 * Calculates the angular correlation \f$W_{\gamma \gamma} \left( \theta,
 * \varphi \right)\f$ between the first and the last photon in a sequence of
 * \f$n-1\f$ (\f$n > 2\f$) electromagnetic (EM) transitions between \f$n\f$
 * states of a quantized system (a 'cascade' of EM transitions with \f$n-1\f$
 * steps). The angles \f$\theta \in \left[ 0, \pi \right]\f$ and \f$\varphi \in
 * \left[ 0, 2 \pi \right]\f$ are the polar and azimuthal angles in a spherical
 * coordinate system, respectively.
 * States of the system are identified by their total angular momentum quantum
 * numbers \f$J\f$ (\f$2J \in \mathcal{N} \f$, \f$J \geq 0\f$) and their
 * parities \f$\pi \in \lbrace -1, 1 \rbrace \f$. They are labeled by indices
 * \f$0 < i \geq n\f$, where \f$i = 1\f$ denotes the first, and \f$i = n\f$ the
 * last state of a cascade. EM transitions are identified by their multipolarity
 * \f$L\f$ (\f$L \in \mathcal{N} \f$, \f$L \geq 0\f$) and their EM character
 * \f$\lambda \in \lbrace \mathrm{E}, \mathrm{M} \rbrace\f$ which can be either
 * electric (E) or magnetic (M). They are labeled by indices \f$1 \geq i < n\f$,
 * where the \f$i\f$-th transition is assumed to be the transition that connects
 * the states labeled \f$i\f$ and \f$i+1\f$. A transition between two states may
 * be a mixture of up to two multipolarities \f$L_i\f$ and \f$L_i^\prime\f$,
 * whose relative contributions are quantified by the corresponding multipole
 * mixing ratio \f$\delta_i\f$ (see below about the convention for \f$\delta\f$
 * which is used in the present implementation).
 *
 * The minimum possible cascade with \f$n = 3\f$ states is denoted symbolically
 * as:
 *
 * \f[
 *      J_1^{\pi_1}
 *      \left( \begin{array}{ccc} L_1 \\ L_1^\prime \end{array} \right)
 *      J_2^{\pi_2}
 *      \left( \begin{array}{ccc} L_2 \\ L_2^\prime \end{array} \right)
 *      J_3^{\pi_3}.
 * \f]
 *
 * A cascade between \f$n > 3\f$ states is denoted as:
 *
 * \f[
 *      J_1^{\pi_1}
 *      \left( \begin{array}{ccc} L_1 \\ L_1^\prime \end{array} \right)
 *      J_2^{\pi_2}
 *      \left( \begin{array}{ccc} L_2 \\ L_2^\prime \end{array} \right)
 *      ...
 *      J_{n-1}^{\pi_{n-1}}
 *      \left( \begin{array}{ccc} L_{n-1} \\ L_{n-1}^\prime \end{array} \right)
 *      J_n^{\pi_n}.
 * \f]
 *
 * The transitions with \f$ 1 < i < n-1\f$ are assumed to be unobserved.
 * The photon which corresponds to the first (see below about the notion of
 * 'first' and 'second') transition is assumed to be observed along the positive
 * z direction (\f$\theta = 0\f$). It may also be interpreted as a beam of
 * photons which travels in that direction and causes the excitation of the
 * state \f$J_2^{\pi_2}\f$ from the initial state \f$J_1^{\pi_1}\f$. Apart from
 * a correlation between the directions of motion of the two photons
 * [direction-direction (dir-dir) correlation], the present implementation can
 * take into account that the polarization of the first photon is observed in
 * addition [a polarization-direction (pol-dir) correlation]. Here, it is
 * assumed that the polarization is along the x axis
 * (\f$\theta = \pi/2\f$, \f$\varphi = 0\f$ denotes the positive x axis).
 * Only the observation of polarization information, which introduces a
 * dependence of the angular correlation on the azimuthal angle \f$\varphi\f$,
 * allows for a distinction between different EM characters, which are also
 * related to the parities of the corresponding states. For this reason, the
 * code will assume a pol-dir correlation if parities and EM characters
 * associated with the first (the 'first' transition is well-defined in the user
 * interface of the code) transition are given, and a dir-dir correlation
 * otherwise.
 *
 * The formalism of angular correlations which was used in the present
 * implementation is mainly based on a review article by Fagg and Hanna \cite
 * FaggHanna1959 and on a book chapter by Biedenharn \cite AjzenbergSelove1960.
 * Consequently, Biedenharn's convention for the multipole mixing ratio is used.
 * For a comparison to the other popular conventions of Rose and Brink and
 * Krane, Steffen, and Wheeler, see Refs. \cite RoseBrink1967 and \cite
 * KraneSteffenWheeler1973, respectively. When a two-step cascade is considered
 * in which the first and the last state are identical, Biedenharn's convention
 * has the advantage that \f$\delta_1 = \delta_2\f$.
 *
 * In the angular correlation formalism, the expansion coefficients of
 * \f$W_{\gamma \gamma}\f$ in terms of Legendre polynomials are separable into
 * contributions by the different transitions. Therefore, a 'first' and 'last'
 * transition of the cascade actually need not be identified. In other words, it
 * does not matter whether the photon observed in z direction with an
 * polarization in x direction is the first or the last one of the cascade.
 * However, the unobserved transitions must take place in between the two,
 * otherwise they would not have to be considered anyway.
 *
 * In order to describe an observation of the 'first' photon in an arbitrary
 * direction with an arbitrary orientation of the polarization in the plane
 * perpendicular to the direction of propagation, the code allows to rotate the
 * angular correlation by three Euler angles \f$\Phi\f$, \f$\Theta\f$, and
 * \f$\Psi\f$. They are defined to be the rotation angles around the z-, x', and
 * z' axes, respectively, as described, for example, in Ref. \cite
 * Weisstein2020.
 *
 * Like many other quantum mechanical computer codes, this one also uses \f$2
 * J\f$ instead of \f$J\f$ and \f$2L\f$ instead of \f$L\f$ to be able to
 * represent both integer and half-integer angular momentum quantum numbers
 * internally as integers.
 *
 * The angular correlation is normalized to \f$4\pi\f$, i.e.:
 *
 * \f[
 *      \int_0^{2\pi} \int_0^\pi W_{\gamma \gamma} \left( \theta, \varphi
 * \right) \sin \left( \theta \right) \mathrm{d} \theta \mathrm{d} \varphi = 4
 * \pi \f]
 */
template <typename T> class AngularCorrelation {

public:
  /**
   * \brief Constructor
   *
   * An EM cascade is specified entirely by the arguments of the constructor.
   * The first argument is the initial state of the cascade, which is sometimes
   * denoted as the 'oriented' state in the literature, because its
   * decay/excitation defines the coordinate system. All other cascade steps are
   * given as a list of pairs of transitions and the states which they populate.
   * The transition between the initial state and the second state, and the one
   * between the \f$n-1\f$-th and the \f$n\f$-th state, are assumed to be the
   * two observed transitions. All other transitions are treated as unobserved.
   *
   * If the parities of the initial state and the second state, and both EM
   * characters of the first transition are given, the code will assume a
   * pol-dir correlation. If none of these data is given, the code will assume a
   * dir-dir correlation.
   *
   * The constructor checks the input data for consistency in terms of angular
   * momentum coupling and selection rules for EM transitions. Further checks
   * are performed by the constructors of the State and Transition classes.
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_ste Cascade steps, given as a list of arbitrary length which
   * contains Transition-State pairs. The first and the last transition of this
   * list are assumed to be observed.
   *
   * \throw invalid_argument if at least one of the following consistency checks
   * fails:
   * 1. The number of cascade steps is larger than one, because two transitions
   * are needed for a gamma-gamma correlation.
   * 2. At least one multipolarity of a transition fulfils the triangle
   * inequality with the angular momenta of the initial and final states of that
   * transition.
   * 3. The EM characters of all multipolarities of a single transition fulfil
   * the parity selection rules.
   * In particular, if at least one parity or EM character is specified, all
   * others must also be specified for the given transition and the states which
   * it links. Note that this also applies to so-called 'pure' transitions where
   * only a single multipolarity is allowed.
   */
  AngularCorrelation(const State ini_sta,
                     const vector<pair<Transition, State>> cas_ste)
      : euler_angle_rotation(EulerAngleRotation<T>()), w_gamma_gamma(nullptr) {
    check_cascade(ini_sta, cas_ste);

    w_gamma_gamma = std::make_unique<W_pol_dir>(ini_sta, cas_ste);
  }

  /**
   * \brief Constructor with transition inference
   *
   * Simplified version of the general constructor which takes a vector of State
   * objects as cascade steps. See the general constructor for more information.
   *
   * For the multipolarities, electromagnetic characters, and mixing ratios of
   * the transitions, the most likely values will be inferred according to the
   * following rules:
   *
   * \enum The multipolarities are assumed to be the lowest possible ones
   * permitted by the triangle inequality for each transitions. \enum The EM
   * character for the transition between the initial state and the second state
   * is set to 'unknown' if any of the parities of the two states is unknown. In
   * this case, a dir-dir correlation will be calculated. If both parities are
   * known, the EM characters corresponding to the lowest possible multipole
   * orders will be chosen. The EM characters for all other transitions are
   * assumed in the same way(Note that the pol-dir correlation only depends on
   * the EM character of the transition for which the polarization information
   * is available.). \enum All multipole-mixing ratios are assumed to be zero.
   *
   * All the aforementioned rules are based on the single-particle 'Weisskopf'
   * estimates (see, e.g., Sec. 6.A in Ref. \cite BlattWeisskopf1979), which
   * favor the lowest multipole order.
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_sta Cascade states, given as a list of arbitrary length which
   * contains State objects. The transitions between the first and the seconds,
   * and between the next-to last and the last state of this list are assumed to
   * be observed.
   *
   * \throw invalid_argument if the number of cascade steps is smaller or equal
   * to one, because two transitions are needed for a gamma-gamma correlation.
   */
  AngularCorrelation(const State ini_sta, const vector<State> cas_sta)
      : AngularCorrelation(ini_sta, gen_cascade_steps(ini_sta, cas_sta)) {}

  /**
   * Copy constructor
   */
  AngularCorrelation(const AngularCorrelation &ang_cor)
      : AngularCorrelation(ang_cor.get_initial_state(),
                           ang_cor.get_cascade_steps()){};

  /**
   * \brief Return the angular correlation for given spherical coordinates.
   *
   * This function assumes that the direction of propagation of the first photon
   * is in the positive z direction. If the correlation is a pol-dir
   * correlation, the function assumes that the polarization axis is the x axis.
   *
   * \param theta Polar angle in spherical coordinates in radians
   * (\f$\theta \in \left[ 0, \pi \right]\f$).
   * \param phi Azimuthal angle in spherical coordinates in radians
   * (\f$\varphi \in \left[ 0, 2 \pi \right]\f$).
   *
   * \return \f$W_{\gamma \gamma} \left( \theta, \varphi \right)\f$
   */
  inline T operator()(const T theta, const T phi) const {
    return w_gamma_gamma->operator()(theta, phi);
  }

  /**
   * \brief Return the angular correlation for an arbitrary coordinate system.
   *
   * This function takes an array of three Euler angles as an addition
   * parameter, to rotate the direction of propagation and the polarization axis
   * (if defined) of the first photon. The 'zxz' convention or 'x' convention is
   * used for the order of the rotations \cite Weisstein2020. If all Euler
   * angles are set to zero, the direction of propagation is in the positive z
   * direction, and the polarization axis (if defined) is the x axis.
   * As implied by the notation, the angles \f$\theta\f$ and \f$\varphi\f$ are
   * still defined in the original coordinate system, i.e. \f$\theta = 0\f$ is
   * still the z axis.
   *
   * \param theta Polar angle in spherical coordinates in radians
   * (\f$\theta \in \left[ 0, \pi \right]\f$).
   * \param phi Azimuthal angle in spherical coordinates in radians
   * (\f$\varphi \in \left[ 0, 2 \pi \right]\f$).
   * \param euler_angles Euler angles \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$
   * in radians which define an arbitrary rotation in 3D space in the x
   * convention.
   *
   * \return \f$W_{\gamma \gamma} \left( \theta, \varphi \right)\f$
   */
  T operator()(const T theta, const T phi,
               const EulerAngles<T> euler_angles) const {
    CoordDir<T> thetap_phip =
        euler_angle_rotation.rotate_back(CoordDir<T>{theta, phi}, euler_angles);

    return (*this)(thetap_phip[0], thetap_phip[1]);
  }

  /**
   * \brief Return the initial state of the angular correlation.
   *
   * \return Initial state.
   */
  inline State get_initial_state() const {
    return w_gamma_gamma->get_initial_state();
  }

  /**
   * \brief Return the cascade steps.
   *
   * \return vector of Transition-State pairs.
   */
  inline vector<pair<Transition, State>> get_cascade_steps() const {
    return w_gamma_gamma->get_cascade_steps();
  }

  /**
   * \brief Return an upper limit for possible values of the gamma-gamma angular
   * correlation.
   *
   * Some applications, for example the rejection-sampling (or 'accept-reject')
   * algorithm (see, e.g. Sec. 2.3 in \cite RobertCasella1999), which can be
   * used to sample random directions that are distributed according to a given
   * angular correlation, require an expression, or at least an estimate, for
   * the maximum absolute value of \f$W \left( \theta, \varphi \right)\f$, i.e.:
   *
   * \f[
   * 		\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in
   * \left[ 0, 2\pi \right]} | W \left( \theta, \varphi \right) |. \f]
   *
   * If a useful upper limit estimate exists for a given angular correlation,
   * this function will return it. If no useful upper limit exists, or the
   * absolute value of \f$W\f$ does not have a limit, this function returns a
   * negative number.
   *
   * This method calls the equivalent method of the W_gamma_gamma member object.
   *
   * \return \f$\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in
   * \left[ 0, 2\pi \right]} | W \left( \theta, \varphi \right) | \f$, or an
   * upper limit for this quantity. If no useful upper limit can be given or if
   * there is no limit, a negative number is returned.
   */
  inline double get_upper_limit() const {
    return w_gamma_gamma->get_upper_limit();
  }

protected:
  /**
   * \brief Check consistency of the input.
   *
   * Checks whether the cascade has more than one step.
   * After that, calls the functions
   * AngularCorrelation::check_triangle_inequalities() and
   * AngularCorrelation::check_em_transitions().
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_ste Cascade steps, given as a list of arbitrary length which
   * contains Transition-State pairs.
   *
   * \throw invalid_argument if the cascade has less than two steps.
   */
  static void check_cascade(const State ini_sta,
                            const vector<CascadeStep> cas_ste) {
    if (cas_ste.size() < 2) {
      throw invalid_argument(
          "Cascade must have at least two transition - state pairs.");
    }

    check_angular_momenta(ini_sta, cas_ste);
    check_triangle_inequalities(ini_sta, cas_ste);
    check_em_transitions(ini_sta, cas_ste);
  }

  /**
   * \brief Check whether angular momenta are either all half integer or all
   * integer.
   *
   * Since photons are particles with a helicity ('spin') of \f$1 \hbar\f$,
   * half-integer momentum transfers in transitions are not possible (see, e.g.,
   * Sec. XXI.IV.26 in \cite Messiah19622 or any textbook on the quantization of
   * electromagnetic radiation). Therefore, a series of states must have
   * uniformly half-integer or integer angular momentum quantum numbers.
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_ste Cascade steps, given as a list of arbitrary length which
   * contains Transition-State pairs.
   *
   * \throw invalid_argument if a mixed use of half-integer and integer angular
   * momentum quantum numbers is detected.
   */
  static void
  check_angular_momenta(const State ini_sta,
                        const vector<pair<Transition, State>> cas_ste) {
    const int even_odd = ini_sta.two_J % 2;

    for (size_t i = 0; i < cas_ste.size(); ++i) {
      if (cas_ste[i].second.two_J % 2 != even_odd) {
        throw invalid_argument(
            "Unphysical mixing of half-integer and integer spins in cascade.");
      }
    }
  }
  /**
   * \brief Check triangle inequality for all cascade steps.
   *
   * Checks whether at least one of
   *
   * \f[
   *      \left| J_i - J_{i+1} \right| \leq L_i \leq J_i + J_{i+1}
   * \f]
   *
   * and
   *
   * \f[
   *      \left| J_i - J_{i+1} \right| \leq L_i^\prime \leq J_i + J_{i+1}
   * \f]
   *
   * is fulfilled for each cascade step.
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_ste Cascade steps, given as a list of arbitrary length which
   * contains Transition-State pairs.
   *
   * \throw invalid_argument if none of the two multipolarities fulfils the
   * triangle inequality for any cascade step.
   */
  static void
  check_triangle_inequalities(const State ini_sta,
                              const vector<pair<Transition, State>> cas_ste) {

    if (!fulfils_triangle_inequality<int>(
            ini_sta.two_J, cas_ste[0].second.two_J, cas_ste[0].first.two_L) &&
        !fulfils_triangle_inequality<int>(
            ini_sta.two_J, cas_ste[0].second.two_J, cas_ste[0].first.two_Lp)) {
      throw invalid_argument(
          "Triangle inequality selection rule not fulfilled for any "
          "multipolarity of transition #1: " +
          cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
    }

    for (size_t i = 1; i < cas_ste.size(); ++i) {
      if (!fulfils_triangle_inequality<int>(cas_ste[i - 1].second.two_J,
                                            cas_ste[i].second.two_J,
                                            cas_ste[i].first.two_L) &&
          !fulfils_triangle_inequality<int>(cas_ste[i - 1].second.two_J,
                                            cas_ste[i].second.two_J,
                                            cas_ste[i].first.two_Lp)) {
        throw invalid_argument(
            "Triangle inequality selection rule not fulfilled for any "
            "multipolarity of transition #" +
            std::to_string(i + 1) + ": " +
            cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
      }
    }
  }
  /**
   * \brief Check parity selections rules for all cascade steps.
   *
   * For each cascade step, first checks whether any of \f$\pi_i\f$,
   * \f$\pi_{i+1}\f$, \f$\lambda_i\f$, \f$\lambda_i^\prime\f$ is given. If
   * yes, checks whether all of them are given. If yes, checks whether the EM
   * selections rules (see, e.g. Appendix 'Electromagnetic Transitions and
   * Moments' in Ref. \cite deShalitTalmi2004)
   *
   * \f[
   *      \pi_i \pi_{i+1} = \left( -1 \right)^{L_i} ~~~ \to ~~~ \lambda_i =
   * \mathrm{E} \f]
   *
   * and
   *
   * \f[
   *      \pi_i \pi_{i+1} = \left( -1 \right)^{L_i+1} ~~~ \to ~~~ \lambda_i =
   * \mathrm{M} \f]
   *
   * are fulfilled.
   * The same are also supposed to hold for \f$L_i^\prime\f$.
   *
   * \param ini_sta Initial state of the cascade.
   * \param cas_ste Cascade steps, given as a list of arbitrary length which
   * contains Transition-State pairs.
   *
   * \throw invalid_argument if at least one of the following conditions is
   * not fulfilled for any cascade step:
   * 1. Either all information about parities and EM characters for a cascade
   * step is given, or none.
   * 2. The parity selection rules apply for both multipolarities of the
   * transition.
   */
  static void
  check_em_transitions(const State ini_sta,
                       const vector<pair<Transition, State>> cas_ste) {
    if (ini_sta.parity != Parity::unknown &&
        cas_ste[0].second.parity != Parity::unknown) {
      if (cas_ste[0].first.em_char != EMCharacter::unknown) {
        if (!valid_em_character(ini_sta.parity, cas_ste[0].second.parity,
                                cas_ste[0].first.two_L,
                                cas_ste[0].first.em_char)) {
          throw invalid_argument(
              "Incorrect electroEMCharacter::magnetic character '" +
              Transition::em_str_rep(cas_ste[0].first.em_char) +
              "' for transition #1: " +
              cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
        }

        if (cas_ste[0].first.em_charp != EMCharacter::unknown) {
          if (!valid_em_character(ini_sta.parity, cas_ste[0].second.parity,
                                  cas_ste[0].first.two_Lp,
                                  cas_ste[0].first.em_charp)) {
            throw invalid_argument(
                "Incorrect electroEMCharacter::magnetic character '" +
                Transition::em_str_rep(cas_ste[0].first.em_charp) +
                "' for transition #1: " +
                cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
          }
        } else {
          throw invalid_argument(
              "Only one electroEMCharacter::magnetic character defined for "
              "transition #1: " +
              cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
        }
      }

      if (cas_ste[0].first.em_char == EMCharacter::unknown &&
          cas_ste[0].first.em_charp != EMCharacter::unknown) {
        throw invalid_argument(
            "Only one electroEMCharacter::magnetic character defined for "
            "transition #1: " +
            cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
      }
    } else if (cas_ste[0].first.em_char != EMCharacter::unknown ||
               cas_ste[0].first.em_charp != EMCharacter::unknown) {
      throw invalid_argument(
          "ElectroEMCharacter::magnetic character defined, but one or both "
          "parities missing "
          "for transition #1: " +
          cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
    }

    for (size_t i = 1; i < cas_ste.size(); ++i) {
      if (cas_ste[i - 1].second.parity != Parity::unknown &&
          cas_ste[i].second.parity != Parity::unknown) {
        if (cas_ste[i].first.em_char != EMCharacter::unknown) {
          if (!valid_em_character(
                  cas_ste[i - 1].second.parity, cas_ste[i].second.parity,
                  cas_ste[i].first.two_L, cas_ste[i].first.em_char)) {
            throw invalid_argument(
                "Incorrect electroEMCharacter::magnetic character '" +
                Transition::em_str_rep(cas_ste[0].first.em_char) +
                "' for transition #" + to_string(i + 1) + ": " +
                cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                         cas_ste[i].second));
          }
        }

        if (cas_ste[i].first.em_charp != EMCharacter::unknown) {
          if (!valid_em_character(
                  cas_ste[i - 1].second.parity, cas_ste[i].second.parity,
                  cas_ste[i].first.two_Lp, cas_ste[i].first.em_charp)) {
            throw invalid_argument(
                "Incorrect electroEMCharacter::magnetic character '" +
                Transition::em_str_rep(cas_ste[0].first.em_charp) +
                "' for transition #" + to_string(i + 1) + ": " +
                cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                         cas_ste[i].second));
          }
        } else {
          throw invalid_argument(
              "Only one electroEMCharacter::magnetic character defined for "
              "transition #" +
              to_string(i + 1) + ": " +
              cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                       cas_ste[i].second));
        }

        if (cas_ste[i].first.em_char == EMCharacter::unknown &&
            cas_ste[i].first.em_charp != EMCharacter::unknown) {
          throw invalid_argument(
              "Only one electroEMCharacter::magnetic character defined for "
              "transition #" +
              to_string(i + 1) + ": " +
              cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                       cas_ste[i].second));
        }
      } else if (cas_ste[i].first.em_char != EMCharacter::unknown ||
                 cas_ste[i].first.em_charp != EMCharacter::unknown) {
        throw invalid_argument(
            "ElectroEMCharacter::magnetic character defined, but one or both "
            "parities missing "
            "for transition #" +
            to_string(i + 1) + ": " +
            cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
      }
    }
  }

  /**
   * \brief Checks the parity selection rule for a single transition.
   *
   * The transition of interest is assumed to have the label \f$i\f$.
   *
   * \param p0 \f$p_i\f$Parity of initial state
   * \param p1 \f$p_{i+1}\f$Parity of final state
   * \param two_L \f$2L_i\f$, two times the multipolarity of the transition.
   * \param em \f$\lambda_i\f$ EM character of the transition with
   * multipolarity \f$L\f$.
   *
   * \return true, if the parity selection rule is fulfilled, false otherwise.
   */
  static bool valid_em_character(const Parity p0, const Parity p1,
                                 const int two_L, const EMCharacter em) {
    if (p0 == p1) {
      if ((two_L / 2) % 2 == 0) {
        if (em != EMCharacter::electric) {
          return false;
        }
      } else if (em != EMCharacter::magnetic) {
        return false;
      }

      return true;
    }

    if ((two_L / 2) % 2 == 0) {
      if (em != EMCharacter::magnetic) {
        return false;
      }
    } else if (em != EMCharacter::electric) {
      return false;
    }

    return true;
  }
  /**
   * \brief Infer the most likely transition that connects two given states.
   *
   * This method is used in connection with the simplified constructor of the
   * AngularCorrelation class. Given two states of a cascade, it infers the
   * most likely transition that connects them. See the corresponding
   * constructor for more information.
   *
   * \param states Two states for which the most likely transition should be
   * inferred.
   *
   * \return Transition object.
   *
   * \throw invalid_argument, if both states have spin 0. In this case, no EM
   * transition is possible.
   */
  Transition infer_transition(const pair<State, State> states) const;

  EulerAngleRotation<T>
      euler_angle_rotation; /**< Instance of the EulerAngleRotation class */

  /**
   * Distribution that is being correlated
   */
  unique_ptr<W_pol_dir> w_gamma_gamma;
};

} // namespace alpaca
