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

#include <cassert>

#include "alpaca/State.hh"
#include "alpaca/Transition.hh"

using namespace alpaca;
using std::invalid_argument;

int main() {
  // Test IO of the Transition class.

  bool error_thrown = false;

  // Error: Both multipolarities are the same.
  // Test for both possible constructors.
  try {
    Transition transition(2, 2, 0.);
  } catch (const invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  try {
    Transition transition(EMCharacter::electric, 2, EMCharacter::magnetic, 2,
                          0.);
  } catch (const invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Check multipolarity IO

  // Error: Multipolarity smaller than zero.
  try {
    Transition transition(EMCharacter::electric, -2, EMCharacter::magnetic, 2,
                          0.);
  } catch (const invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Multipolarity zero.
  try {
    Transition transition(EMCharacter::electric, 0, EMCharacter::magnetic, 2,
                          0.);
  } catch (const invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Check string representation

  Transition transition(2, 4, 0.);
  try {
    transition.em_str_rep(EMCharacter::unknown);
  } catch (const runtime_error &e) {
    error_thrown = true;
  }

  assert(error_thrown);

  auto s0 = State(0, Parity::unknown);
  auto s0p = State(0, Parity::positive);
  auto s1 = State(2, Parity::unknown);
  auto s1p = State(2, Parity::positive);
  auto s1m = State(2, Parity::negative);
  auto s2 = State(4, Parity::unknown);
  auto s2p = State(4, Parity::positive);
  auto s2m = State(4, Parity::negative);

  auto trans_E1 = Transition(s0p, s1m);
  assert(trans_E1 == Transition::E1());

  auto trans_M1 = Transition(s0p, s1p);
  assert(trans_M1 == Transition::M1());

  auto trans_E2 = Transition(s0p, s2p);
  assert(trans_E2 == Transition::E2());

  auto trans_M2 = Transition(s0p, s2m);
  assert(trans_M2 == Transition::M2());

  auto trans_D1_0 = Transition(s0, s0);
  assert(trans_D1_0 == Transition::Dipole());

  auto trans_D1 = Transition(s2, s1);
  assert(trans_D1 == Transition::Dipole());

  auto trans_D2 = Transition(s0, s2);
  assert(trans_D2 == Transition::Quadrupole());

  auto trans_M1_L0 = Transition(s1p, s1p);
  assert(trans_M1_L0 == Transition::M1());

  auto trans_E1_L0 = Transition(s1p, s1m);
  assert(trans_E1_L0 == Transition::E1());

  auto trans_M1_0 = Transition(s0p, s0p);
  assert(trans_M1_0 == Transition::M1());
}
