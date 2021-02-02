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
#include <gsl/gsl_sf.h>

#include "UvCoefficient.hh"

double UvCoefficient::operator()(const unsigned int two_nu, const int two_j, const int two_L, const int two_jp) const {

    // Definition of Fagg and Hanna \cite FaggHanna1959 [Eq. (I-1') and the expression below that one].
    // Causes some tests to fail.
    //
    // const int phase_factor = (((two_jp - two_j - two_L)/2) % 2) == 0 ? 1 : -1;

    // return phase_factor
    // *sqrt(
    //     (two_jp + 1)*(two_j + 1)
    // )
    // *gsl_sf_coupling_6j(
    //     two_j, two_j, two_nu,
    //     two_jp, two_jp, two_L);

    // Definition of Biedenharn \cite AjzenbergSelove1960 (Sec. 1.a.1.iii)
    const int phase_factor = (((two_j + two_jp + two_L)/2) % 2) == 0 ? 1 : -1;

    return phase_factor
    *sqrt(
        (two_jp + 1)*(two_j + 1)
    )
    *gsl_sf_coupling_6j(
        two_j, two_nu, two_j,
        two_jp, two_L, two_jp);
}

double UvCoefficient::operator()(const unsigned int two_nu, const int two_j, const int two_L, const int two_Lp, const double delta, const int two_jp) const {

    if(delta != 0.){
        return (*this)(two_nu, two_j, two_L, two_jp) 
            + delta*delta*(*this)(two_nu, two_j, two_Lp, two_jp);
    }

    return (*this)(two_nu, two_j, two_L, two_jp);

}