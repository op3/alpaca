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

#include <array>
#include <numeric>
#include <utility>
#include <vector>

using std::accumulate;
using std::array;
using std::vector;

#include "alpaca/SphereRejectionSampler.hh"

namespace alpaca {

double SphereRejectionSampler::estimate_efficiency(const unsigned int n_tries) {
  vector<unsigned int> required_tries(n_tries);

  pair<unsigned int, array<double, 2>> sampled_theta_phi;

  for (unsigned int i = 0; i < n_tries; ++i) {
    sampled_theta_phi = sample();
    required_tries[i] = sampled_theta_phi.first;
  }

  return static_cast<double>(n_tries) /
         static_cast<double>(
             accumulate(required_tries.begin(), required_tries.end(), 0));
}

} // namespace alpaca
