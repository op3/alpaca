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

#include <stdexcept>
#include <string>

using std::invalid_argument;
using std::to_string;

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/W_dir_dir.hh"
#include "alpaca/W_pol_dir.hh"

namespace alpaca {

extern "C" {
double angular_correlation(const double theta, const double phi,
                           const size_t n_cas_ste, int *two_J, short *par,
                           short *em_char, int *two_L, short *em_charp,
                           int *two_Lp, double *delta, double *PhiThetaPsi) {
  State initial_state{two_J[0], static_cast<Parity>(par[0])};
  vector<CascadeStep> cascade_steps;

  for (size_t i = 0; i < n_cas_ste; ++i) {
    cascade_steps.push_back(
        {Transition{static_cast<EMCharacter>(em_char[i]), two_L[i],
                    static_cast<EMCharacter>(em_charp[i]), two_Lp[i], delta[i]},
         State{two_J[i + 1], static_cast<Parity>(par[i + 1])}});
  }

  return AngularCorrelation<double>(initial_state, cascade_steps)
      .
      operator()(theta, phi, {PhiThetaPsi[0], PhiThetaPsi[1], PhiThetaPsi[2]});
}

void *create_angular_correlation(const size_t n_cas_ste, int *two_J, short *par,
                                 short *em_char, int *two_L, short *em_charp,
                                 int *two_Lp, double *delta) {

  State initial_state{two_J[0], static_cast<Parity>(par[0])};
  vector<CascadeStep> cascade_steps;

  for (size_t i = 0; i < n_cas_ste; ++i) {
    cascade_steps.push_back(
        {Transition{static_cast<EMCharacter>(em_char[i]), two_L[i],
                    static_cast<EMCharacter>(em_charp[i]), two_Lp[i], delta[i]},
         State{two_J[i + 1], static_cast<Parity>(par[i + 1])}});
  }

  return new AngularCorrelation<double>(initial_state, cascade_steps);
}

void *
create_angular_correlation_with_transition_inference(const size_t n_cas_ste,
                                                     int *two_J, short *par) {

  State initial_state{two_J[0], static_cast<Parity>(par[0])};
  vector<State> cascade_states;

  for (size_t i = 0; i < n_cas_ste; ++i) {
    cascade_states.push_back(
        {State{two_J[i + 1], static_cast<Parity>(par[i + 1])}});
  }

  return new AngularCorrelation<double>(initial_state, cascade_states);
}

void evaluate_angular_correlation(
    AngularCorrelation<double> *angular_correlation, const size_t n_angles,
    double *theta, double *phi, double *result) {

  for (size_t i = 0; i < n_angles; ++i) {
    result[i] = angular_correlation->operator()(theta[i], phi[i]);
  }
}

void evaluate_angular_correlation_rotated(
    AngularCorrelation<double> *angular_correlation, const size_t n_angles,
    double *theta, double *phi, double *PhiThetaPsi, double *result) {

  for (size_t i = 0; i < n_angles; ++i) {
    result[i] = angular_correlation->operator()(
        theta[i], phi[i], {PhiThetaPsi[0], PhiThetaPsi[1], PhiThetaPsi[2]});
  }
}

void free_angular_correlation(AngularCorrelation<double> *angular_correlation) {
  delete angular_correlation;
}

void get_em_char(AngularCorrelation<double> *angular_correlation,
                 short *em_char) {

  vector<CascadeStep> cascade_steps = angular_correlation->get_cascade_steps();
  for (size_t i = 0; i < cascade_steps.size(); ++i) {
    em_char[i] = static_cast<short>(cascade_steps[i].first.em_char);
  }
}

void get_two_L(AngularCorrelation<double> *angular_correlation, int *two_L) {

  vector<CascadeStep> cascade_steps = angular_correlation->get_cascade_steps();
  for (size_t i = 0; i < cascade_steps.size(); ++i) {
    two_L[i] = cascade_steps[i].first.two_L;
  }
}
}

} // namespace alpaca
