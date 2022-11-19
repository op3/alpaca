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

#include "alpaca/State.hh"

namespace alpaca {

string State::parity_str_rep() const {

  if (parity == Parity::positive) {
    return "+";
  }
  if (parity == Parity::negative) {
    return "-";
  }

  throw std::runtime_error("No string representation for unknown parity.");
}

string State::spin_str_rep() const {

  if (two_J % 2 == 0) {
    return std::to_string(two_J / 2);
  }

  return std::to_string(two_J) + "/2";
}

string State::str_rep() const {

  if (parity != Parity::unknown) {
    return spin_str_rep() + "^" + parity_str_rep();
  }

  return spin_str_rep();
}

int State::check_two_J(const int two_J) {

  if (two_J < 0) {
    throw std::invalid_argument("two_J must be a nonnegative integer.");
  }

  return two_J;
}

double State::check_excitation_energy(const double e_x) {

  if (e_x < 0.) {
    throw std::invalid_argument("Excitation energy must not be negative.");
  }

  return e_x;
}

} // namespace alpaca
