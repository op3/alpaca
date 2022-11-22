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

    Copyright (C) 2022 Udo Friman-Gayer, Oliver Papst
*/

#include <cassert>
#include <cmath>
#include <iostream>

#include "alpaca/Special.hh"
#include "alpaca/TestUtilities.hh"

using namespace alpaca;

void test_legendre() {
  for (unsigned int l = 0; l < 10; ++l) {
    test_numerical_equality(legendre(l, 0.4), std::legendre(l, 0.4), 1e-10);
  }
}

void test_assoc_legendre() {
  for (unsigned int l = 0; l <= 4; ++l) {
    for (unsigned int m = 0; m <= l; ++m) {
      test_numerical_equality(assoc_legendre(l, m, 0.4),
                              std::assoc_legendre(l, m, 0.4), 1e-10);
    }
  }
  test_numerical_equality(assoc_legendre(5, 0, 0.4),
                          std::assoc_legendre(5, 0, 0.4), 1e-10);
  test_numerical_equality(assoc_legendre(5, 2, 0.4),
                          std::assoc_legendre(5, 2, 0.4), 1e-10);
  test_numerical_equality(assoc_legendre(6, 0, 0.4),
                          std::assoc_legendre(6, 0, 0.4), 1e-10);
  test_numerical_equality(assoc_legendre(6, 2, 0.4),
                          std::assoc_legendre(6, 2, 0.4), 1e-10);
}

int main() {
  test_legendre();
  test_assoc_legendre();
}
