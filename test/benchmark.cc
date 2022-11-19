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

#include <chrono>
#include <iostream>

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/CascadeRejectionSampler.hh"
#include "alpaca/State.hh"
#include "alpaca/Transition.hh"

using namespace alpaca;

/**
 * This script tests the correctness of SphereRejectionSampler by sampling from
 * a distribution that is 1 for \f$\varphi < \pi\f$ and 0 otherwise. For an
 * optimum choice \f$W_\mathrm{max} = 1\f$, this should result in an efficiency
 * of \f$\epsilon = 0.5\f$. For the choice \f$W_\mathrm{max} = 2\f$, this should
 * result in an efficiency of \f$\epsilon = 0.25\f$.
 *
 * It also tests the Euler-angle rotation functionality of
 * SphereRejectionSampler by sampling two vectors, one from an unrotated
 * distribution and another one from a rotated distribution. By using the same
 * random number seed and two different SphereRejectionSampler objects, it is
 * ensured that the random vector before the rotation was the same as the
 * unrotated one. By manually applying the rotation using the EulerAngleRotation
 * class, the unrotated vector is transformed into the rotated vector's
 * reference frame and both are tested for equality.
 */
int main() {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  vector<AngularCorrelation> cascade{
      AngularCorrelation(
          State(0, Parity::positive),
          {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
            State(2, Parity::positive)},
           {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
            State(0, Parity::unknown)}}),
  };

  static constexpr int random_number_seed = 1;

  CascadeRejectionSampler cas_rej_sam =
      CascadeRejectionSampler(cascade, random_number_seed, {0., 0., 0.}, false);

  double result = 1.;
  size_t count = 10000000;

  auto start = high_resolution_clock::now();
  for (size_t i = 0; i < count; i++) {
    vector<array<double, 2>> transitions = cas_rej_sam();
    result += transitions[0][0];
    result += transitions[0][1];
  }
  auto stop = high_resolution_clock::now();

  std::cout << "Result: " << result << "\n";
  std::cout << "Time for " << count << " iterations: "
            << duration_cast<milliseconds>(stop - start).count() << " ms\n";
}
