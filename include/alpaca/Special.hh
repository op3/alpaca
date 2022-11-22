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

    Copyright (C) 2021-2022 Udo Friman-Gayer, Oliver Papst
*/

#include <numbers>
#include <type_traits>

using std::numbers::pi;

#include <enoki/array.h>
#include <enoki/special.h>

template <typename T, typename U> struct legendre_;

template <typename T, int L>
struct legendre_<T, std::integral_constant<int, L>> {
  constexpr static T val(T x) {
    return ((2 * L - 1) * x *
                legendre_<T, std::integral_constant<int, L - 1>>::val(x) -
            (L - 1) *
                legendre_<T, std::integral_constant<int, L - 2>>::val(x)) /
           L;
  }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 6>> {
  constexpr static T val(T x) {
    return 0.0625 * (231. * enoki::pow(x, 6) - 315. * enoki::pow(x, 4) +
                     105. * enoki::pow(x, 2) - 5.);
  }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 5>> {
  constexpr static T val(T x) {
    return 0.125 * (63. * enoki::pow(x, 5) - 70. * enoki::pow(x, 3) + 15 * x);
  }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 4>> {
  constexpr static T val(T x) {
    return 0.125 * (35. * enoki::pow(x, 4) - 30. * enoki::pow(x, 2) + 3);
  }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 3>> {
  constexpr static T val(T x) { return 0.5 * (5 * enoki::pow(x, 3) - 3 * x); }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 2>> {
  constexpr static T val(T x) { return 0.5 * (3 * enoki::pow(x, 2) - 1); }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 1>> {
  constexpr static T val(T x) { return x; }
};

template <typename T> struct legendre_<T, std::integral_constant<int, 0>> {
  constexpr static T val([[maybe_unused]] T x) { return 1; }
};

template <int L, typename T> T legendre(T x) {
  return legendre_<T, std::integral_constant<int, L>>::val(x);
}

/**
 * \brief Calculate legendre polynomial
 *
 * This templated function allows for the calculation of legendre polynomials
 * using any type of data. The low-degree functions up to l=6 are hard-coded in
 * a look-up-table, starting from l=5, recursion is used for calculation.
 * The Condon-Shortley phase term is omitted, just like in the STL and GSL
 * implementations.
 *
 * \tparam T data type used for the calculation.
 * \param l degree of the polynomial
 * \param x evaluate at this position
 *
 * \return Value of the legendre polynomial at x
 */
template <typename T> T legendre(unsigned int l, T x) {
  switch (l) {
  case 0:
    return legendre<0>(x);
  case 1:
    return legendre<1>(x);
  case 2:
    return legendre<2>(x);
  case 3:
    return legendre<3>(x);
  case 4:
    return legendre<4>(x);
  case 5:
    return legendre<5>(x);
  case 6:
    return legendre<6>(x);
  }
  return ((2 * l - 1) * x * legendre(l - 1, x) - (l - 1) * legendre(l - 2, x)) /
         l;
}

template <typename T, typename U, typename V> struct assoc_legendre_;

template <typename T, int L, int M>
struct assoc_legendre_<T, std::integral_constant<int, L>,
                       std::integral_constant<int, M>> {
  constexpr static T val(T x) { return 0 * x; }
};

#define DEF_ASSOC_LEGENDRE(N, M, func)                                         \
  template <typename T>                                                        \
  struct assoc_legendre_<T, std::integral_constant<int, N>,                    \
                         std::integral_constant<int, M>> {                     \
    constexpr static T val([[maybe_unused]] T x) { return func; }              \
  }

DEF_ASSOC_LEGENDRE(0, 0, 1.);
DEF_ASSOC_LEGENDRE(1, 0, x);
DEF_ASSOC_LEGENDRE(1, 1, enoki::sqrt(1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(2, 0, 0.5 * (3. * enoki::sqr(x) - 1.));
DEF_ASSOC_LEGENDRE(2, 1, 3. * x * enoki::sqrt(1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(2, 2, 3 * (1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(3, 0, 0.5 * (5. * enoki::pow(x, 3) - 3. * x));
DEF_ASSOC_LEGENDRE(3, 1,
                   -1.5 * (1. - 5. * enoki::sqr(x)) *
                       enoki::sqrt(1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(3, 2, 15. * x * (1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(3, 3, 15. * enoki::pow(1 - enoki::sqr(x), 1.5));
DEF_ASSOC_LEGENDRE(4, 0,
                   0.125 * (35. * enoki::pow(x, 4) - 30. * enoki::sqr(x) + 3.));
DEF_ASSOC_LEGENDRE(4, 1,
                   2.5 * (7. * enoki::pow(x, 3) - 3 * x) *
                       enoki::sqrt(1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(4, 2, 7.5 * (7 * enoki::sqr(x) - 1.) * (1. - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(4, 3, 105. * x * enoki::pow(1. - enoki::sqr(x), 1.5));
DEF_ASSOC_LEGENDRE(4, 4, 105. * enoki::sqr(1 - enoki::sqr(x)));
DEF_ASSOC_LEGENDRE(5, 0,
                   0.125 * (63. * enoki::pow(x, 5) - 70. * enoki::pow(x, 3) +
                            15. * x));
DEF_ASSOC_LEGENDRE(5, 2,
                   -52.5 * (3. * enoki::pow(x, 5) - 4. * enoki::pow(x, 3) + x));
DEF_ASSOC_LEGENDRE(6, 0,
                   0.0625 * (231. * enoki::pow(x, 6) - 315. * enoki::pow(x, 4) +
                             105. * enoki::sqr(x) - 5));
DEF_ASSOC_LEGENDRE(6, 2,
                   -13.125 * (33. * enoki::pow(x, 6) - 51. * enoki::pow(x, 4) +
                              19. * enoki::sqr(x) - 1.));

template <int L, int M, typename T> T assoc_legendre(T x) {
  return assoc_legendre_<T, std::integral_constant<int, L>,
                         std::integral_constant<int, M>>::val(x);
}

constexpr unsigned int switch_pair(unsigned int l, unsigned int m) {
  return (l << 8) + m;
}

/**
 * \brief Calculate associated legendre polynomial
 *
 * Only implemented for l < 4 and 0 ≤ m ≤ 4.
 * This templated function allows for the calculation of legendre polynomials
 * using any type of data. The functions are hard-coded in a look-up-table, no
 * recursions take place.
 * The Condon-Shortley phase term is omitted, just like in the STL and GSL
 * implementations.
 *
 * \tparam T data type used for the calculation.
 * \param l degree of the polynomial
 * \param m order of the polynomial
 * \param x evaluate at this position
 *
 * \return Value of the associated legendre polynomial at x
 */
template <typename T> T assoc_legendre(unsigned int l, unsigned int m, T x) {
  switch (switch_pair(l, m)) {
  case switch_pair(0, 0):
    return assoc_legendre<0, 0>(x);
  case switch_pair(1, 0):
    return assoc_legendre<1, 0>(x);
  case switch_pair(1, 1):
    return assoc_legendre<1, 1>(x);
  case switch_pair(2, 0):
    return assoc_legendre<2, 0>(x);
  case switch_pair(2, 1):
    return assoc_legendre<2, 1>(x);
  case switch_pair(2, 2):
    return assoc_legendre<2, 2>(x);
  case switch_pair(3, 0):
    return assoc_legendre<3, 0>(x);
  case switch_pair(3, 1):
    return assoc_legendre<3, 1>(x);
  case switch_pair(3, 2):
    return assoc_legendre<3, 2>(x);
  case switch_pair(3, 3):
    return assoc_legendre<3, 3>(x);
  case switch_pair(4, 0):
    return assoc_legendre<4, 0>(x);
  case switch_pair(4, 1):
    return assoc_legendre<4, 1>(x);
  case switch_pair(4, 2):
    return assoc_legendre<4, 2>(x);
  case switch_pair(4, 3):
    return assoc_legendre<4, 3>(x);
  case switch_pair(4, 4):
    return assoc_legendre<4, 4>(x);
  case switch_pair(5, 0):
    return assoc_legendre<5, 0>(x);
  case switch_pair(5, 2):
    return assoc_legendre<5, 2>(x);
  case switch_pair(6, 0):
    return assoc_legendre<6, 0>(x);
  case switch_pair(6, 2):
    return assoc_legendre<6, 2>(x);
  }
  return std::numeric_limits<double>::quiet_NaN() * x;
}

/**
 * \brief Elliptic integral of the second kind \f$E\left( \varphi | m
 * \right)\f$ for arbitrary real parameters
 *
 * The length of the spiral trajectory, \f$S\left( \Theta \right)\f$, on which
 * points are sampled in Ref. \cite Koay2011 [Eq. (4) therein] can be
 * expressed by incomplete elliptic integrals of the second kind \f$E\left(
 * \varphi | m \right)\f$ {Eq. (5) in Ref. \cite Koay2011}:
 *
 * \f[
 *      E\left( \varphi | m \right) = \int_0^{\varphi} \sqrt{ 1 - m \left[
 * \sin \left( \theta \right) \right]^2 } \mathrm{d} \theta. \f]
 *
 * where \f$\varphi\f$ is an angle, and \f$m\f$ is assumed to be a negative
 * real number
 * [\f$m = -c^2\f$, see Eq. (4)] in Ref. \cite Koay2011.
 * A definition of the elliptic integrals can be found in Sec. 17 of the
 * 'Handbook of Mathematical Functions' by Abramowitz and Stegun \cite
 * AbramowitzStegun1974. It indicates that \f$m\f$ may even be complex [see,
 * e.g. the discussion below Eq. (17.2.18) therein]. The functions \f$E\left(
 * \varphi | \theta \right)\f$ are implemented in the GNU Scientific Library
 * (GSL) \cite Galassi2009 only for \f$ 0 \leq m < 1 \f$ with a parameter
 * \f$k\f$,
 *
 * \f[
 *      k^2 = m.
 * \f]
 *
 * However, for any real parameter \f$m\f$ that is negative or has an absolute
 * value larger than 1, the evaluation of the integral can be traced back to
 * an evaluation in the interval \f$ 0 \leq m < 1 \f$ according to Eqs.
 * (17.4.16) and (17.4.18) in \cite AbramowitzStegun1974. These
 * transformations assume a formulation of the incomplete elliptic integral in
 * terms of the parameter \f$u\f$ instead of \f$\varphi\f$, i.e. \f$E \left( u
 * | m \right)\f$. The quantity \f$u\f$ is the generalized angle \f$\varphi\f$
 * for the case of an ellipse , such that, for example, an elliptic equivalent
 * \f$\mathrm{sn}\f$ {a 'Jacobi Elliptic Function', see, e.g., Sec. 16 in Ref.
 * \cite abromowitzStegun1974} of the trigonometric function \f$\sin\f$ can be
 * defined {Eq. (17.2.2) in Ref. \cite AbramowitzStegun1974}:
 *
 * \f[
 *      \mathrm{sn} \left( u \right) = \sin \left( \varphi \right).
 * \f]
 *
 * The equation above provides a straightforward conversion from the Jacobi
 * elliptic function \f$\mathrm{sn} \left( u \right)\f$ to \f$\varphi\f$ via
 * the inverse sine. The Jacobi elliptic functions are also implemented in GSL
 * for the argument \f$u\f$, but it should be noted that this implementation
 * utilizes the parameter \f$m\f$ instead of \f$k\f$. The conversion from
 * \f$\varphi\f$ to \f$u\f$ is possible via the definition of the incomplete
 * elliptic integral of the first kind {Eq. (17.2.7) in Ref. \cite
 * AbramowitzStegun1974}:
 *
 * \f[
 *      F \left( \varphi | m \right) = u.
 * \f]
 *
 * For an arbitrary \f$m\f$, the function
 * elliptic_integral_2nd_kind_arbitrary_m will be called
 * recursively with transformed parameters until the parameter can be handled
 * by the GSL implementation.
 *
 * In fact, the implementation does not only use the formalism of Abramowitz
 * and Stegun, but the transformation to the interval \f$ 0 \leq m < 1 \f$ of
 * the incomplete elliptic integral of the second kind is performed using Eq.
 * (19.7.5) ('Imaginary Modulus Transformation') of the NIST Digital Library
 * of Mathematical Functions (DLMF) \cite DLMF2020. This equation has the
 * advantage that it does not require the back-and-forth conversion between
 * \f$\varphi\f$ and \f$u\f$, and that it handles \f$m < 0\f$ and \f$|m| >
 * 1\f$ at the same time.
 *
 * \param phi \f$\varphi\f$
 * \param m \f$m\f$
 *
 * \return \f$E\left( \varphi | m \right)\f$
 */
inline double ellint_2_arbitrary_m(const double phi, const double m) {

  // Use Eq. (19.6.9) in Ref. \cite DLMF2020.
  // The enoki implementation cannot handle this case.
  if (m == 1.) {
    return sin(phi);
  }

  const double abs_m = enoki::abs(m);
  const double abs_k = enoki::sqrt(abs_m);

  // Direct call of the library function possible.
  if (m >= 0. && abs_m < 1.) {
    return enoki::ellint_2(phi, abs_k);
  }

  // Use Eq. (19.7.5) in Ref. \cite DLMF2020.
  const double abs_k_squared = enoki::sqr(abs_k);
  const double kappa_prime = enoki::rsqrt(1. + abs_k_squared);
  const double kappa = abs_k * kappa_prime;
  const double kappa_squared = enoki::sqr(kappa);
  // phi == pi/2 corresponds to the complete elliptic integral
  // This is both a shortcut and an insurance against numerical instabilities.
  if (phi == 0.5 * pi) {
    return enoki::comp_ellint_2(kappa) / kappa_prime;
  }

  const double theta = enoki::asin(
      enoki::sin(phi) /
      (kappa_prime *
       enoki::sqrt(1. + abs_k * abs_k * enoki::sqr(enoki::sin(phi)))));

  return 1. / kappa_prime *
         (ellint_2_arbitrary_m(theta, kappa_squared) -
          kappa_squared * enoki::sin(theta) * enoki::cos(theta) *
              enoki::rsqrt(1. - kappa_squared * enoki::sqr(enoki::sin(theta))));

  // Alternative implementation of the cases \f$m < 0\f$ and \f$ |m| > 1 \f$
  // based on Eqs. (17.4.16) and (17.4.18) in \cite AbramowitzStegun1974, which
  // did not give the correct results. Possibly something went wrong in the
  // transformation between \f$\varphi\f$ and \f$u\f$.
  //
  // const double u = elliptic_integral_1st_kind_arbitrary_m(phi, m);
  //
  // if(m < 0.){
  //     const double sqrt_1_plus_abs_m = sqrt(1. + abs_m);
  //     const double abs_m_over_1_plus_abs_m = abs_m/(1. + abs_m);

  //     double sn_u_times_sqrt_1_plus_abs_m,
  //     cn_u_times_sqrt_1_plus_abs_m, dn_u_times_sqrt_1_plus_abs_m;

  //     gsl_sf_elljac_e(u*sqrt_1_plus_abs_m, abs_m_over_1_plus_abs_m,
  //     &sn_u_times_sqrt_1_plus_abs_m, &cn_u_times_sqrt_1_plus_abs_m,
  //     &dn_u_times_sqrt_1_plus_abs_m);

  //     return sqrt_1_plus_abs_m*(
  //         elliptic_integral_2nd_kind_arbitrary_m(
  //             asin(sn_u_times_sqrt_1_plus_abs_m),
  //             abs_m_over_1_plus_abs_m)
  //         -abs_m/sqrt_1_plus_abs_m*sn_u_times_sqrt_1_plus_abs_m*cn_u_times_sqrt_1_plus_abs_m/dn_u_times_sqrt_1_plus_abs_m
  //     );
  // }

  // const double inverse_abs_m = 1./abs_m;

  // double sn_u_times_sqrt_abs_m, cn_u_times_sqrt_abs_m, dn_u_times_sqrt_abs_m;

  // gsl_sf_elljac_e(u*sqrt_abs_m, inverse_abs_m, &sn_u_times_sqrt_abs_m,
  // &cn_u_times_sqrt_abs_m, &dn_u_times_sqrt_abs_m);

  // return sqrt_abs_m*elliptic_integral_2nd_kind_arbitrary_m(
  //     asin(sn_u_times_sqrt_abs_m),
  //     inverse_abs_m)
  // -(abs_m-1.)*u;
}

/**
 * \brief Elliptic integral of the first kind \f$F\left( \varphi | m
 * \right)\f$ for arbitrary real parameters
 *
 * Needed for the determination of the optimum value for \f$c\f$ according to
 * Eq. (14) in Ref. \cite Koay2011. The transformation from arbitrary real
 * \f$m\f$ to the range \f$0 \leq m < 1\f$ is given by Eqs. (19.7.5) in Ref.
 * \cite DLMF2020. See also the definition of
 * elliptic_integral_2nd_kind_arbitrary_m for more
 * information.
 *
 * \param phi \f$\varphi\f$
 * \param m \f$m\f$
 *
 * \return \f$F\left( \varphi | m \right)\f$
 */
inline double ellint_1_arbitrary_m(const double phi, const double m) {
  const double abs_m = enoki::abs(m);
  const double abs_k = enoki::sqrt(abs_m);

  if (m >= 0. && abs_m < 1.) {
    return enoki::ellint_1(phi, abs_k);
  }

  // Use Eq. (19.7.5) in Ref. \cite DLMF2020.
  const double abs_k_squared = enoki::sqr(abs_k);
  const double kappa_prime = enoki::rsqrt(1. + abs_k_squared);
  const double kappa = abs_k * kappa_prime;

  // phi == pi/2 corresponds to the complete elliptic integral
  // This is both a shortcut and an insurance against numerical instabilities.
  if (phi == 0.5 * pi) {
    return kappa_prime * enoki::comp_ellint_1(kappa);
  }

  const double theta = enoki::asin(
      enoki::sin(phi) /
      (kappa_prime *
       enoki::sqrt(1. + abs_k_squared * enoki::sqr(enoki::sin(phi)))));

  return kappa_prime * ellint_1_arbitrary_m(theta, kappa * kappa);
}
