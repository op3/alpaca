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

#include <gsl/gsl_sf.h>

#include "FCoefficient.hh"
#include "Transition.hh"

/**
 * \brief Class for an \f$E_\nu\f$ coefficient.
 * 
 * The coefficients \f$E_\nu\f$ are a generalization of the \f$A_\nu\f$ coefficients that are 
 * encountered in the definition of the dir-dir correlation {Eq. (I-2) in \cite FaggHanna1959}.
 * They are introduced in a textbook chapter by Biedenharn \cite AjzenbergSelove1960 to take into 
 * account a linear polarization of a photon and given by the expression 
 * {Eq. (85a) in Ref. \cite AjzenbergSelove1960}:
 * 
 * \f[
 * 		E_\nu \left( n \right) = \Bigg \lbrace \left( -1 \right)^{\sigma_{L_n}} F_\nu \left( L_n, L_n, j_n, j \right) \frac{2 \nu \left( \nu + 1 \right) L_n \left( L_n + 1 \right)}{\nu \left( \nu + 1 \right) - 2 L_n \left( L_n + 1 \right)}
 * \f]
 * \f[
 * 		+ 2 \delta_n \left( -1 \right)^{\sigma_{L_n^\prime}} F_\nu \left( L_n, L_n^\prime, j_n, j \right) \left( L_1^\prime - L_1 \right) \left( L_1^\prime + L_1 +1 \right)
 * \f]
 * \f[
 * 		+ \delta_n^2 \left( -1 \right)^{\sigma_{L_n^\prime}} F_\nu \left( L_n^\prime, L_n^\prime, j_n, j \right) \frac{2\nu \left( \nu + 1 \right) L_n^\prime \left( L_n^\prime + 1 \right)}{\nu \left( \nu + 1 \right) - 2 L_n^\prime \left( L_n^\prime + 1 \right)} \Bigg \rbrace
 * \f]
 * \f[
 * 		\times \frac{\left( \nu - 2 \right)!}{\left( \nu + 2 \right)!}
 * \f]
 * 
 * In the equation above, the symbol \f$\sigma_L\f$ is equal to \f$1\f$ if the transition with 
 * multipolarity \f$L\f$ has magnetic character, and \f$0\f$ for electric character.
 * Note that Eq. (85) in Ref. \cite AjzenbergSelove1960, which first introduces the coefficients 
 * \f$E_\nu\f$, includes an additional factor \f$\sin \left( \beta \right)\f$, where the angle
 * \f$\beta\f$ denotes the polar angle between the photon's direction of motion and its 
 * polarization. Here, \f$\beta = 90^\circ\f$, or \f$\sin \left( \beta \right) = 1\f$ was assumed.
 * 
 * See also the definition of the AvCoefficient class for more information.
 */

class EvCoefficient{
public:
	/**
	 * \brief Constructor
	 * 
	 * \param two_nu \f$2 \nu\f$
	 * \param em Primary electromagnetic character \f$\sigma_{L_n}\f$
	 * \param two_L Primary multipolarity \f$2 L\f$
	 * \param emp Secondary electromagnetic character \f$\sigma_{L_n^\prime}\f$
	 * \param two_Lp Secondary multipolarity \f$2 L^\prime\f$
	 * \param two_jn Angular momentum quantum number \f$2 j_n\f$ of the initial or final state
	 * 	of a transition 
	 * \param two_j Angular momentum quantum number \f$2 j\f$ of the intermediate state
	 * 	of a transition 
	 */
	EvCoefficient(const int two_nu, const EMCharacter em, const int two_L, const EMCharacter emp, const int two_Lp, const int two_jn, const int two_j) :
		two_nu(two_nu),
        em(em),
        two_L(two_L),
        emp(emp),
        two_Lp(two_Lp),
        two_jn(two_jn),
        two_j(two_j),
        sign_sigma_L_n((em == magnetic) ? -1 : 1),
        sign_sigma_Lp_n((emp == magnetic) ? -1 : 1),
        constant_f_coefficient(two_nu, two_L, two_L, two_jn, two_j),
        linear_f_coefficient(two_nu, two_L, two_Lp, two_jn, two_j),
		quadratic_f_coefficient(two_nu, two_Lp, two_Lp, two_jn, two_j)
	{}

	/**
	 * \brief Return value of a specific \f$E_\nu\f$ coefficient.
	 *
	 * \param delta Multipole mixing ratio \f$\delta\f$
	 *
	 * \return \f$E_\nu \left( L, L^\prime, j_n, j, \delta_n \right)\f$
	 */
	double operator()(const double delta) const {
		const int nu = two_nu/2;
		const double nu_times_nu_plus_one = nu*(nu+1);
		const int L = two_L/2;
		const double two_L_times_L_plus_one = 2*L*(L+1);
		const int Lp = two_Lp/2;
		const double two_Lp_times_Lp_plus_one = 2*Lp*(Lp+1);

		return (
			sign_sigma_L_n*constant_f_coefficient.get_value()*(nu_times_nu_plus_one*two_L_times_L_plus_one)/(nu_times_nu_plus_one - two_L_times_L_plus_one)
			+ 2.*delta*sign_sigma_Lp_n*linear_f_coefficient.get_value()*(Lp - L)*(Lp + L + 1)
			+ delta*delta*sign_sigma_Lp_n*quadratic_f_coefficient.get_value()*(nu_times_nu_plus_one*two_Lp_times_Lp_plus_one)/(nu_times_nu_plus_one - two_Lp_times_Lp_plus_one)
		)*gsl_sf_fact(nu-2)/gsl_sf_fact(nu+2);
	}

protected:
	int two_nu;
	EMCharacter em;
	int two_L;
	EMCharacter emp;
	int two_Lp;
	int two_jn;
	int two_j;
	int sign_sigma_L_n;
	int sign_sigma_Lp_n;

	FCoefficient constant_f_coefficient, linear_f_coefficient, quadratic_f_coefficient;
};
