#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <alpaca/AngularCorrelation.hh>
#include <alpaca/CascadeRejectionSampler.hh>
#include <alpaca/State.hh>
#include <alpaca/Transition.hh>

namespace py = pybind11;
using namespace alpaca;

PYBIND11_MODULE(alpyca, m) {
  py::enum_<EMCharacter>(m, "EMCharacter")
      .value("electric", EMCharacter::electric)
      .value("magnetic", EMCharacter::magnetic)
      .value("unknown", EMCharacter::unknown);
  py::class_<Transition>(m, "Transition")
      .def(py::init<int, int, double>())
      .def(py::init<EMCharacter, int, EMCharacter, int, double>())
      .def_readwrite("em_char", &Transition::em_char)
      .def_readwrite("two_L", &Transition::two_L)
      .def_readwrite("em_charp", &Transition::em_charp)
      .def_readwrite("two_Lp", &Transition::two_Lp)
      .def_readwrite("delta", &Transition::delta);

  py::enum_<Parity>(m, "Parity")
      .value("negative", Parity::negative)
      .value("positive", Parity::positive)
      .value("parity_unknown", Parity::unknown);
  py::class_<State>(m, "State")
      .def(py::init<int>())
      .def(py::init<int, Parity>())
      .def(py::init<int, Parity, double>())
      .def(py::init<int, double>())
      .def_readwrite("two_J", &State::two_J)
      .def_readwrite("parity", &State::parity)
      .def_readwrite("excitation_energy", &State::excitation_energy)
      .def("__repr__", &State::str_rep);

  py::class_<AngularCorrelation<double>>(m, "AngularCorrelation")
      .def(py::init<State, std::vector<std::pair<Transition, State>>>())
      .def(py::init<State, std::vector<State>>())
      .def("calc", py::vectorize(py::overload_cast<const double, const double>(
                       &AngularCorrelation<double>::operator(), py::const_)))
      .def("calc", py::vectorize(py::overload_cast<const double, const double,
                                                   const EulerAngles<double>>(
                       &AngularCorrelation<double>::operator(), py::const_)));

  py::class_<CascadeRejectionSampler<double>>(m, "CascadeRejectionSampler")
      .def(py::init<std::vector<AngularCorrelation<double>> &, int,
                    unsigned int>())
      .def(py::init<std::vector<AngularCorrelation<double>> &, int,
                    EulerAngles<double>, bool>())
      .def("calc", &CascadeRejectionSampler<double>::operator())
      .def("calc", [](CascadeRejectionSampler<double> &c, py::ssize_t len) {
        auto dim = static_cast<py::ssize_t>(c.size());
        auto res = py::array_t<double>(std::vector<ptrdiff_t>{
            static_cast<ptrdiff_t>(len), static_cast<ptrdiff_t>(dim), 2});

        auto data = res.mutable_unchecked<3>();

        for (py::ssize_t i = 0; i < len; ++i) {
          auto sample = c();
          for (py::ssize_t j = 0; j < dim; ++j) {
            data(i, j, 0) = sample[j][0];
            data(i, j, 1) = sample[j][1];
          }
        }
        return res;
      });
}
