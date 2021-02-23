#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <boost/format.hpp>
#include <string>
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "eos_thermal.h"
#include "eos_thermal_file.h"


namespace etk = EOS_Toolkit;
namespace py = pybind11;

using real_t = etk::real_t;
using range = etk::interval<real_t>;



std::string range_to_str(const range& rg) 
{
  boost::format fmt("min=%.15e, max=%.15e");
  fmt % rg.min() % rg.max();
  return fmt.str();
}


PYBIND11_MODULE(pyreprimand, m) {
    m.doc() = "Python bindings for RePrimAnd library";
    
    py::class_<range>(m, "range")
        .def(py::init<real_t, real_t>(),
             py::arg("min")=0, py::arg("max")=0)
        .def("limit_to", py::vectorize(&range::limit_to),
             "Limit values to range")
        .def("contains", py::vectorize(&range::contains),
             "Test if values are within interval [min,max]")
        .def_property_readonly("min", &range::min)
        .def_property_readonly("max", &range::max)
        .def_property_readonly("length", &range::length)
        .def("__str__", &range_to_str);

    py::class_<etk::units>(m, "units")
        .def(py::init<const double, const double, const double>(),
             py::arg("length"), py::arg("time"), py::arg("mass"))
        .def(py::self / py::self)
        .def_property_readonly("length", &etk::units::length)
        .def_property_readonly("time", &etk::units::time)
        .def_property_readonly("freq", &etk::units::freq)
        .def_property_readonly("mass", &etk::units::mass)
        .def_property_readonly("velocity", &etk::units::velocity)
        .def_property_readonly("accel", &etk::units::accel)
        .def_property_readonly("force", &etk::units::force)
        .def_property_readonly("area", &etk::units::area)
        .def_property_readonly("volume", &etk::units::volume)
        .def_property_readonly("density", &etk::units::density)
        .def_property_readonly("pressure", &etk::units::pressure)
        .def_property_readonly("mom_inertia", &etk::units::mom_inertia)
        .def_readonly_static("c_SI", &etk::units::c_SI)
        .def_readonly_static("G_SI", &etk::units::G_SI)
        .def_readonly_static("M_sun_SI", &etk::units::M_sun_SI)
        .def_property_readonly_static("geom_solar", 
                 [](py::object) { return etk::units::geom_solar(); })
        .def_property_readonly_static("geom_meter", 
                 [](py::object) { return etk::units::geom_meter(); })
        .def("__str__", &etk::units::to_str);
    
    
    py::class_<etk::eos_thermal>(m, "eos_thermal")
        .def("press_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::press_at_rho_eps_ye),
             "Compute pressure from density, specific energy, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("csnd_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::csnd_at_rho_eps_ye),
             "Compute soundspeed from density, specific energy, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("temp_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::temp_at_rho_eps_ye),
             "Compute temperature from density, specific energy, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("sentr_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::sentr_at_rho_eps_ye),
             "Compute specific entropy from density, specific energy, "
             "and electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("dpress_drho_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::dpress_drho_at_rho_eps_ye),
             "Compute partial derivative of pressure with respect to "
             "density, given density, specific energy, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("dpress_deps_at_rho_eps_ye", 
             py::vectorize(&etk::eos_thermal::dpress_deps_at_rho_eps_ye),
             "Compute partial derivative of pressure with respect to "
             "specific energy, given density, specific energy, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("eps"),py::arg("ye"))
        .def("press_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::press_at_rho_temp_ye),
             "Compute pressure from density, temperature, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("csnd_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::csnd_at_rho_temp_ye),
             "Compute soundspeed from density, temperature, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("eps_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::eps_at_rho_temp_ye),
             "Compute temperature from density, temperature, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("sentr_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::sentr_at_rho_temp_ye),
             "Compute specific entropy from density, temperature, "
             "and electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("dpress_drho_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::dpress_drho_at_rho_temp_ye),
             "Compute partial derivative of pressure with respect to "
             "density, given density, temperature, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("dpress_deps_at_rho_temp_ye", 
             py::vectorize(&etk::eos_thermal::dpress_deps_at_rho_temp_ye),
             "Compute partial derivative of pressure with respect to "
             "specific energy, given density, temperature, and "
             "electron fraction.\n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"),py::arg("temp"),py::arg("ye"))
        .def("is_rho_valid", 
             py::vectorize(&etk::eos_thermal::is_rho_valid),
             "Check if given densities are within valid range.",
             py::arg("rho"))
        .def("is_ye_valid", 
             py::vectorize(&etk::eos_thermal::is_ye_valid),
             "Check if given electron fractions are within valid range.",
             py::arg("ye"))
        .def("is_rho_ye_valid", 
             py::vectorize(&etk::eos_thermal::is_rho_ye_valid),
             "Check if given densities and electron fractions are "
             "within valid range.",
             py::arg("rho"), py::arg("ye"))
        .def("is_rho_eps_ye_valid", 
             py::vectorize(&etk::eos_thermal::is_rho_eps_ye_valid),
             "Check if given densities, specific energies, and "
             "electron fractions are within valid range.",
             py::arg("rho"), py::arg("eps"), py::arg("ye"))
        .def("is_rho_temp_ye_valid", 
             py::vectorize(&etk::eos_thermal::is_rho_temp_ye_valid),
             "Check if given densities, temperatures, and "
             "electron fractions are within valid range.",
             py::arg("rho"), py::arg("temp"), py::arg("ye"))
        .def("range_eps", 
             &etk::eos_thermal::range_eps,
             "Get valid range of specific energy, given density and "
             "electron fraction.\n"
             "Note: this method does not accept arrays and requires "
             "valid density and electron fraction\n",
             py::arg("rho"),py::arg("ye"))
        .def("range_temp", 
             &etk::eos_thermal::range_temp,
             "Get valid range of temperature, given density and "
             "electron fraction.\n"
             "Note: this method does not accept arrays and requires "
             "valid density and electron fraction\n",
             py::arg("rho"),py::arg("ye"))
        .def_property_readonly("range_rho", 
             &etk::eos_thermal::range_rho,
             "Validity range for density")
        .def_property_readonly("range_ye", 
             &etk::eos_thermal::range_ye,
             "Validity range for electron fraction")
        .def_property_readonly("minimal_h", 
              &etk::eos_thermal::minimal_h,
             "Lower bound for enthalpy");

    m.def("load_eos_thermal", &etk::load_eos_thermal, 
          "Load thermal EOS from file",
          py::arg("path"),
          py::arg("units")=etk::units::geom_solar());


    py::class_<etk::eos_barotr>(m, "eos_barotr")
        .def("gm1_at_rho", 
             py::vectorize(&etk::eos_barotr::gm1_at_rho),
             "Compute pseudo enthalpy g-1 from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("press_at_rho", 
             py::vectorize(&etk::eos_barotr::press_at_rho),
             "Compute pressure from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("eps_at_rho", 
             py::vectorize(&etk::eos_barotr::eps_at_rho),
             "Compute specific internal energy from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("hm1_at_rho", 
             py::vectorize(&etk::eos_barotr::hm1_at_rho),
             "Compute enthalpy h - 1 from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("csdn_at_rho", 
             py::vectorize(&etk::eos_barotr::csnd_at_rho),
             "Compute soundspeed from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("temp_at_rho", 
             py::vectorize(&etk::eos_barotr::temp_at_rho),
             "Compute temperature from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("ye_at_rho", 
             py::vectorize(&etk::eos_barotr::ye_at_rho),
             "Compute electron fraction from density \n"
             "returns NAN outside EOS validity region.",
             py::arg("rho"))
        .def("rho_at_gm1", 
             py::vectorize(&etk::eos_barotr::rho_at_gm1),
             "Compute mass density from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("press_at_gm1", 
             py::vectorize(&etk::eos_barotr::press_at_gm1),
             "Compute pressure from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("eps_at_gm1", 
             py::vectorize(&etk::eos_barotr::eps_at_gm1),
             "Compute specific internal energy from pseudo "
             "enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("hm1_at_gm1", 
             py::vectorize(&etk::eos_barotr::hm1_at_gm1),
             "Compute enthalpy h - 1 from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("csdn_at_gm1", 
             py::vectorize(&etk::eos_barotr::csnd_at_gm1),
             "Compute soundspeed from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("temp_at_gm1", 
             py::vectorize(&etk::eos_barotr::temp_at_gm1),
             "Compute temperature from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def("ye_at_gm1", 
             py::vectorize(&etk::eos_barotr::ye_at_gm1),
             "Compute electron fraction from pseudo enthalpy g-1 \n"
             "returns NAN outside EOS validity region.",
             py::arg("gm1"))
        .def_property_readonly("is_isentropic", 
             &etk::eos_barotr::is_isentropic,
             "Whether EOS is isentropic")
        .def_property_readonly("is_zero_temp", 
             &etk::eos_barotr::is_zero_temp,
             "Whether EOS is zero-temperature")
        .def_property_readonly("has_temp", 
             &etk::eos_barotr::has_temp,
             "Whether EOS provides temperatures")
        .def_property_readonly("has_efrac", 
             &etk::eos_barotr::has_efrac,
             "Whether EOS provides electron fraction")
        .def_property_readonly("range_rho", 
             &etk::eos_barotr::range_rho,
             "Validity range for density")
        .def_property_readonly("range_gm1", 
             &etk::eos_barotr::range_gm1,
             "Validity range for pseudo enthalpy g - 1")
        .def("is_rho_valid", 
             py::vectorize(&etk::eos_barotr::is_rho_valid),
             "Check if given densities are within valid range.",
             py::arg("rho"))
        .def("is_gm1_valid", 
             py::vectorize(&etk::eos_barotr::is_gm1_valid),
             "Check if given pseudo enthalpies g - 1 are within "
             "valid range.", 
             py::arg("rho"))
        .def_property_readonly("minimal_h", 
              &etk::eos_barotr::minimal_h,
             "Lower bound for enthalpy");
             
    m.def("load_eos_barotr", &etk::load_eos_barotr, 
          "Load barotropic EOS from file",
          py::arg("path"),
          py::arg("units")=etk::units::geom_solar());

}
