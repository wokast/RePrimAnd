#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <boost/format.hpp>
#include <string>
#include <vector>
#include "reprimand/interpol.h"
#include "reprimand/eos_barotropic.h"
#include "reprimand/eos_barotr_file.h"
#include "reprimand/eos_barotr_table.h"
#include "reprimand/eos_barotr_spline.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/eos_barotr_pwpoly.h"
#include "reprimand/eos_thermal.h"
#include "reprimand/eos_thermal_file.h"
#include "reprimand/spherical_stars.h"
#include "reprimand/star_sequence.h"
#include "reprimand/star_seq_file.h"


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

std::string units_to_repr(const etk::units& u) 
{
  boost::format fmt("pyreprimand.units(%.20e, %.20e, %.20e)");
  fmt % u.length() % u.time() % u.mass();
  return fmt.str();
}


PYBIND11_MODULE(pyreprimand, m) {
    m.doc() = "Python bindings for RePrimAnd library";
    
    py::class_<range>(m, "range", "Represents a numeric range")
        .def(py::init<real_t, real_t>(), 
             "Construct from minimum and maximum",
             py::arg("min")=0, py::arg("max")=0)
        .def("limit_to", py::vectorize(&range::limit_to),
             "Limit values to range")
        .def("contains", 
             py::vectorize([](const range* that, real_t x) { 
               return that->contains(x); 
             }),
             "Test if values are within interval [min,max]")
        .def_property_readonly("min", &range::min)
        .def_property_readonly("max", &range::max)
        .def_property_readonly("length", &range::length)
        .def("__str__", &range_to_str);

    py::class_<etk::units>(m, "units", "Represents a unit system")
        .def(py::init<const double, const double, const double>(),
R"(Construct arbitrary unit system from base units length, 
time, and mass)",
             py::arg("length"), 
             py::arg("time"), 
             py::arg("mass"))
        .def(py::self / py::self, 
"Dividing unit objects A/B expresses units of A in unit system B")
        .def_property_readonly("length", &etk::units::length, 
                               "Unit of Length")
        .def_property_readonly("time", &etk::units::time, 
                               "Unit of time")
        .def_property_readonly("freq", &etk::units::freq, 
                               "Unit of frequency")
        .def_property_readonly("mass", &etk::units::mass, 
                               "Unit of mass")
        .def_property_readonly("velocity", &etk::units::velocity, 
                               "Unit of velocity")
        .def_property_readonly("accel", &etk::units::accel, 
                               "Unit of acceleration")
        .def_property_readonly("force", &etk::units::force, 
                               "Unit of force")
        .def_property_readonly("area", &etk::units::area, 
                               "Unit of area")
        .def_property_readonly("volume", &etk::units::volume, 
                               "Unit of volume")
        .def_property_readonly("density", &etk::units::density, 
                               "Unit of mass density")
        .def_property_readonly("pressure", &etk::units::pressure, 
                               "Unit of pressure")
        .def_property_readonly("mom_inertia", &etk::units::mom_inertia, 
                               "Unit of moment of inertia")
        .def_readonly_static("c_SI", &etk::units::c_SI,
                             "Speed of light in SI units")
        .def_readonly_static("G_SI", &etk::units::G_SI,
R"(Default value for gravitational constant in SI units.
Used in some methods unless other value is specified)")
        .def_readonly_static("M_sun_SI", &etk::units::M_sun_SI,
R"(Default value for solar mass in SI units.
Used in some methods unless other value is specified)")
        .def_static("geom_udensity", &etk::units::geom_udensity,
R"(Create geometric unit system with given mass density unit and 
"G=c=1. The constant G in SI units can be overridden.)",
             py::arg("umass"),
             py::arg("g_si")=etk::units::G_SI)
        .def_static("geom_umass", &etk::units::geom_umass, 
R"(Create geometric unit system (i.e. G=c=1) specified by 
mass unit. The constant G [SI units] can be overridden.)",
             py::arg("umass"),
             py::arg("g_si")=etk::units::G_SI)
        .def_static("geom_ulength", &etk::units::geom_ulength,
R"(Create geometric unit system (i.e. G=c=1) specified by
length unit. The constant G [SI units] can be overridden.)",
             py::arg("ulength"),
             py::arg("g_si")=etk::units::G_SI)
        .def_static("geom_solar", &etk::units::geom_solar,
R"(Create geometric unit system with G=c=M_sun=1.
Constants for G and M_sun in SI units can be overridden.)",
             py::arg("msun_si")=etk::units::M_sun_SI,
             py::arg("g_si")=etk::units::G_SI)
        .def_static("geom_meter", &etk::units::geom_meter,
R"(Create geometric unit system with G=c=1 and a length unit
of 1 m. The constant G [SI units] can be overridden.)",
             py::arg("g_si")=etk::units::G_SI)
        .def("__repr__", &units_to_repr)
        .def("__str__", &etk::units::to_str);
    
    
    py::class_<etk::eos_thermal>(m, "eos_thermal", 
    "Represents an EOS with thermal and composition effects")
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
             "Lower bound for enthalpy")
        .def_property_readonly("units_to_SI", 
              &etk::eos_thermal::units_to_SI,
             "Unit system (geometric) used by the EOS")
        .def("__str__", &etk::eos_thermal::descr_str);


    m.def("load_eos_thermal", &etk::load_eos_thermal, 
          "Load thermal EOS from file",
          py::arg("path"),
          py::arg("units")=etk::units::geom_solar());


    py::class_<etk::eos_barotr>(m, "eos_barotr",
        "Represents a barotropic EOS")
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
        .def("csnd_at_rho", 
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
        .def("csnd_at_gm1", 
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
             "Lower bound for enthalpy")
        .def_property_readonly("units_to_SI", 
              &etk::eos_barotr::units_to_SI,
             "Unit system (geometric) used by the EOS")
        .def("__str__", &etk::eos_barotr::descr_str);
             
    m.def("load_eos_barotr", &etk::load_eos_barotr, 
R"(Load barotropic EOS from file

Args:
    path (str): Path of the EOS file
    units (pyreprimand.units): Unit system, with respect to SI, to be 
        used by the returned EOS object (this has nothing to do with 
        the units inside the file, which are fixed)

Returns:
    pyreprimand.eos_barotr object
    
)",
          py::arg("path"),
          py::arg("units")=etk::units::geom_solar());

    m.def("save_eos_barotr", &etk::save_eos_barotr, 
R"(Save barotropic EOS to file

Args:
    path (str): Path of the EOS file (file should not exist)
    eos (pyreprimand.eos_barotr): The EOS object to be saved
    info (str): Optional arbitrary free-form string to be embedded 
        into the file (this is not used by the library)

)",
          py::arg("path"),
          py::arg("eos"),
          py::arg("info")="");


    m.def("make_eos_barotr_table",
          &etk::make_eos_barotr_table,
R"(DEPRECATED: create EOS based on linear interpolation from irregularly
spaced sample points. Use make_eos_barotr_spline instead.
)",
          py::arg("gm1"),
          py::arg("rho"),
          py::arg("eps"),
          py::arg("p_by_rho"),
          py::arg("csndsqr"),
          py::arg("temp"),
          py::arg("efrac"), 
          py::arg("isentropic"),            
          py::arg("n_poly"),
          py::arg("units")=etk::units::geom_solar());   



    m.def("make_eos_barotr_spline",
          [] (const std::vector<real_t>& gm1,
              const std::vector<real_t>& rho,
              const std::vector<real_t>& eps,
              const std::vector<real_t>& press,
              const std::vector<real_t>& csnd,
              const std::vector<real_t>& temp,
              const std::vector<real_t>& efrac, 
              bool isentr, range rg_rho,
              real_t n_poly, etk::units uc,
              std::size_t pts_per_mag) {
                return etk::make_eos_barotr_spline(gm1, rho, eps,
                press, csnd, temp, efrac, isentr, rg_rho, n_poly, 
                uc, pts_per_mag);
              },
R"(Create an EOS based on interpolation splines from irregularly spaced
sample points. This EOS type is internally using monotonic splines that
are regular in log-space. The given arbitrarily-spaced samples are first
resampled via irregular-spaced monotonic splines. Note that the cost of
evaluating the EOS is (almost) independent of the internal resolution.

Also note that the pseudo-enthalpy  g-1 is redundant, it is the 
responsability of the user to provide consistent values. Other versions 
exist that compute g-1 (or even epsilon and g-1) from the other 
quantities. Use those unless accurate values for g-1 and epsilon are 
readily available, e.g. for EOS given by analytic expressions. 

Args:
    gm1 (array_like): Pseudo enthalpy g-1
    rho (array_like): Baryonic mass density
    eps (array_like): Specific energy
    press (array_like): Pressure
    csnd (array_like): Adiabatic sound speed
    temperature (array_like): Temperature (or empty array if unavailable)
    efrac (array_like): Electron fraction (or empty array if unavailable)
    isentropic (bool): Whether EOS is supposed to be isentropic
    rg_rho (pyreprimand.range): Range of mass density where EOS is 
        computed via spline interpolation. 
    n_poly (float): Adiabatic index of the generalized polytropic EOS
        that will be used below the interpolated range
    units (pyreprimand.units): Unit system of the EOS and the given 
        sample points.
    pts_per_mag (int): How many sample points per magnitude should be
        used by the interpolation splines employed internally by the 
        EOS. Note this has nothing to do with the resolution of the 
        provided sample points. The default should be reasonable for
        most applications.

)",
              
          py::arg("gm1"),
          py::arg("rho"),
          py::arg("eps"),
          py::arg("press"),
          py::arg("csnd"),
          py::arg("temp"),
          py::arg("efrac"), 
          py::arg("isentropic"),
          py::arg("rg_rho"),            
          py::arg("n_poly"),
          py::arg("units")=etk::units::geom_solar(),
          py::arg("pts_per_mag")=200);   

    m.def("make_eos_barotr_spline",
          [] (const std::vector<real_t>& rho,
              const std::vector<real_t>& eps,
              const std::vector<real_t>& press,
              const std::vector<real_t>& csnd,
              const std::vector<real_t>& temp,
              const std::vector<real_t>& efrac, 
              bool isentr, range rg_rho,
              real_t n_poly, etk::units uc,
              std::size_t pts_per_mag) {
                return etk::make_eos_barotr_spline(rho, eps,
                press, csnd, temp, efrac, isentr, rg_rho, n_poly, 
                uc, pts_per_mag);
              },
R"(Create an EOS based on interpolation splines from irregularly spaced
sample points. This EOS type is internally using monotonic splines that
are regular in log-space. The given arbitrarily-spaced samples are first
resampled via irregular-spaced monotonic splines. Note that the cost of
evaluating the EOS is (almost) independent of the internal resolution.

This version of the function computes the required pseudo-enthalpy g-1 
from given density, pressure, and specific energy via numerical 
integration of the corresponding interpolation splines (thus ensuring 
consistency). 

Args:
    rho (array_like): Baryonic mass density
    eps (array_like): Specific energy
    press (array_like): Pressure
    csnd (array_like): Adiabatic sound speed
    temperature (array_like): Temperature (or empty array if unavailable)
    efrac (array_like): Electron fraction (or empty array if unavailable)
    isentropic (bool): Whether EOS is supposed to be isentropic
    rg_rho (pyreprimand.range): Range of mass density where EOS is 
        computed via spline interpolation. 
    n_poly (float): Adiabatic index of the generalized polytropic EOS
        that will be used below the interpolated range
    units (pyreprimand.units): Unit system of the EOS and the given 
        sample points.
    pts_per_mag (int): How many sample points per magnitude should be
        used by the interpolation splines employed internally by the 
        EOS. Note this has nothing to do with the resolution of the 
        provided sample points. The default should be reasonable for
        most applications.

)",              
          py::arg("rho"),
          py::arg("eps"),
          py::arg("press"),
          py::arg("csnd"),
          py::arg("temp"),
          py::arg("efrac"), 
          py::arg("isentropic"),
          py::arg("rg_rho"),            
          py::arg("n_poly"),
          py::arg("units")=etk::units::geom_solar(),
          py::arg("pts_per_mag")=200);   

    m.def("make_eos_barotr_spline",
          [] (const std::vector<real_t>& rho,
              const std::vector<real_t>& press,
              const std::vector<real_t>& csnd,
              const std::vector<real_t>& temp,
              const std::vector<real_t>& efrac, 
              range rg_rho,
              real_t n_poly, 
              real_t eps_0, 
              etk::units uc,
              std::size_t pts_per_mag) {
                return etk::make_eos_barotr_spline(rho, press, csnd, 
                temp, efrac, rg_rho, n_poly, eps_0,
                uc, pts_per_mag);
              },
R"(Create an EOS based on interpolation splines from irregularly spaced
sample points. This EOS type is internally using monotonic splines that
are regular in log-space. The given arbitrarily-spaced samples are first
resampled via irregular-spaced monotonic splines. Note that the cost of
evaluating the EOS is (almost) independent of the internal resolution.

This version of the function computes the required specific energy from
pressure and density under the assumption that the EOS is isentropic.
This is done via numerical integration of the corresponding interpolation 
splines (thus ensuring consistency). The required pseudo-enthalpy g-1 is 
similarly computed from specific energy and given density and pressure.

Args:
    rho (array_like): Baryonic mass density
    press (array_like): Pressure
    csnd (array_like): Adiabatic sound speed
    temperature (array_like): Temperature (or empty array if unavailable)
    efrac (array_like): Electron fraction (or empty array if unavailable)
    rg_rho (pyreprimand.range): Range of mass density where EOS is 
        computed via spline interpolation. 
    n_poly (float): Adiabatic index of the generalized polytropic EOS
        that will be used below the interpolated range
    eps_0 (float): Specific energy density at zero density
    units (pyreprimand.units): Unit system of the EOS and the given 
        sample points.
    pts_per_mag (int): How many sample points per magnitude should be
        used by the interpolation splines employed internally by the 
        EOS. Note this has nothing to do with the resolution of the 
        provided sample points. The default should be reasonable for
        most applications.

)",              
          py::arg("rho"),
          py::arg("press"),
          py::arg("csnd"),
          py::arg("temp"),
          py::arg("efrac"), 
          py::arg("rg_rho"),            
          py::arg("n_poly"),
          py::arg("eps_0"),
          py::arg("units")=etk::units::geom_solar(),
          py::arg("pts_per_mag")=200);   

    m.def("make_eos_barotr_spline",
          [] (etk::eos_barotr eos, 
              range rg_rho, real_t n_poly, 
              std::size_t pts_per_mag) {
                return etk::make_eos_barotr_spline(eos, rg_rho, 
                n_poly, pts_per_mag);
              },
R"(Create an EOS based on interpolation splines from an existing 
barotropic EOS of (any type) via sampling. 

Args:
    eos (pyreprimand.eos_barotr): The EOS to be approximated
    rg_rho (pyreprimand.range): Range of mass density where EOS is 
        computed via spline interpolation. 
    n_poly (float): Adiabatic index of the generalized polytropic EOS
        that will be used below the interpolated range
    pts_per_mag (int): How many sample points per magnitude should be
        used by the interpolation splines employed internally by the 
        EOS. 
        
Returns:
    eos_barotr_spline object (using same units as the source EOS)

)",              
          py::arg("eos"),
          py::arg("rg_rho"),            
          py::arg("n_poly"),
          py::arg("pts_per_mag")=200);   

    m.def("make_eos_barotr_pwpoly", 
          &etk::make_eos_barotr_pwpoly,
R"(Create a picewise polytropic EOS.

Args:
    rho_poly_0 (float): The polytropic density scale of the first segment.
    rho_segm_bounds (array_like): Mass densities of lower boundary 
        for each segment. The first segment has to start at zero.
    segm_gammas (array_like): Adiabatic exponent for each segment.
    rho_max (float): Maximum valid mass density.
    units (pyreprimand.units): Unit system of the EOS and the parameters
        above.

Returns:
    pyreprimand.eos_barotr object.

)",
          py::arg("rho_poly_0"),
          py::arg("rho_segm_bounds"),
          py::arg("segm_gammas"),
          py::arg("rho_max"),
          py::arg("units")=etk::units::geom_solar());

    m.def("make_eos_barotr_poly",
          &etk::make_eos_barotr_poly, 
R"(Create an polytropic EOS.

Args:
    n_poly (float): Polytropic index
    rho_poly (float): The polytropic density scale.
    rho_max (float): Maximum valid mass density.
    units (pyreprimand.units): Unit system of the EOS and the parameters
        above.

Returns:
    pyreprimand.eos_barotr object.

)", 
          py::arg("n_poly"), 
          py::arg("rho_poly"),
          py::arg("rho_max"),
          py::arg("units")=etk::units::geom_solar());

    py::class_<etk::interpolator>(m, "interpolator")
        .def("__call__", py::vectorize(&etk::interpolator::operator()))
        .def_property_readonly("range_x", 
             &etk::interpolator::range_x,
             "Valid range")
        .def_property_readonly("range_y", 
             &etk::interpolator::range_y,
             "Value range");
        
    m.def("make_interpol_regspl",
          [] (std::vector<real_t> y, etk::interval<real_t> rg) {
            return etk::make_interpol_regspl(y,rg);
          },
          py::arg("y"),
          py::arg("rg"));


    m.def("make_interpol_regspl",
          [] (const etk::interpolator f, etk::interval<real_t> rg, 
              std::size_t ns) {
            return etk::make_interpol_regspl(f,rg,ns);
          },
          py::arg("f"),
          py::arg("rg"),
          py::arg("nsamp"));

    m.def("make_interpol_logspl",
          [] (std::vector<real_t> y, etk::interval<real_t> rg) {
            return etk::make_interpol_logspl(y,rg);
          },
          py::arg("y"),
          py::arg("rg"));


    m.def("make_interpol_logspl",
          [] (const etk::interpolator f, etk::interval<real_t> rg, 
              std::size_t ns) {
            return etk::make_interpol_logspl(f,rg,ns);
          },
          py::arg("f"),
          py::arg("rg"),
          py::arg("nsamp"));


    m.def("make_interpol_llogspl",
          [] (std::vector<real_t> y, etk::interval<real_t> rg) {
            return etk::make_interpol_llogspl(y,rg);
          },
          py::arg("y"),
          py::arg("rg"));


    m.def("make_interpol_llogspl",
          [] (const etk::interpolator f, etk::interval<real_t> rg, 
              std::size_t ns) {
            return etk::make_interpol_llogspl(f,rg,ns);
          },
          py::arg("f"),
          py::arg("rg"),
          py::arg("nsamp"));


    m.def("make_interpol_pchip_spline",
          [] (const std::vector<real_t>& x, 
              const std::vector<real_t>& y) 
          {
            return etk::make_interpol_pchip_spline(x,y);
          },
          py::arg("x"),
          py::arg("y"));

    py::class_<etk::spherical_star_tidal>(m, "spherical_star_tidal",
    "Describes tidal deformability" )
        .def_readonly("k2", &etk::spherical_star_tidal::k2,
                      "Love number k_2")
        .def_readonly("lambda_tidal", 
                      &etk::spherical_star_tidal::lambda,
                      "Dimensionless tidal deformability Lambda");

    
    py::class_<etk::spherical_star_profile>(m, "spherical_star_profile",
R"(Represents the radial profiles of a spherical star

The profiles are provided as functions of circumferential
radius that can be evaluated both inside and outside the star.
All quantities are in the same geometric units used by the EOS.
)" )
        .def_property_readonly("eos", 
             &etk::spherical_star_profile::eos,
             "Star EOS")
        .def_property_readonly("center_gm1", 
             &etk::spherical_star_profile::center_gm1,
             "Central pseudo enthalpy g - 1")
        .def_property_readonly("surf_circ_radius", 
             &etk::spherical_star_profile::surf_circ_radius,
             "Surface circumferential radius")
        .def("nu_from_rc", 
             py::vectorize(&etk::spherical_star_profile::nu_from_rc),
             "Metric function nu at given circumf. radius",
             py::arg("rc"))
        .def("lambda_from_rc", 
             py::vectorize(&etk::spherical_star_profile::lambda_from_rc),
             "Metric function lambda at circumf. radius",
             py::arg("rc"))
        .def("gm1_from_rc", 
             py::vectorize(&etk::spherical_star_profile::gm1_from_rc),
             "Pseudo enthalpy g - 1 at circumf. radius",
             py::arg("rc"))
        .def("mbary_from_rc", 
             py::vectorize(&etk::spherical_star_profile::mbary_from_rc),
             "Baryonic mass within circumf. radius",
             py::arg("rc"))
        .def("pvol_from_rc", 
             py::vectorize(&etk::spherical_star_profile::pvol_from_rc),
             "Proper volume within circumf. radius",
             py::arg("rc"));
  
        
    py::class_<etk::spherical_star_properties>(m, 
            "spherical_star_properties",
R"(Represents properties of a spherical neutron star

This collects scalar measures of a spherical neutron star
solution, but not the radial profile. It also provides the EOS.
All quantities are in the same geometric units used by the EOS.

)" )
        .def_property_readonly("eos", 
             &etk::spherical_star_properties::eos,
             "Star EOS")
        .def_property_readonly("center_rho", 
             &etk::spherical_star_properties::center_rho,
             "Central baryonic mass density")
        .def_property_readonly("center_eps", 
             &etk::spherical_star_properties::center_eps,
             "Central specific internal energy")
        .def_property_readonly("center_press", 
             &etk::spherical_star_properties::center_press,
             "Central pressure")
        .def_property_readonly("center_csnd", 
             &etk::spherical_star_properties::center_csnd,
             "Central sound speed")
        .def_property_readonly("center_ye", 
             &etk::spherical_star_properties::center_ye,
             "Central electron fraction")
        .def_property_readonly("grav_mass", 
             &etk::spherical_star_properties::grav_mass,
             "Star gravitational mass")
        .def_property_readonly("bary_mass", 
             &etk::spherical_star_properties::bary_mass,
             "Star baryonic mass")
        .def_property_readonly("binding_energy", 
             &etk::spherical_star_properties::binding_energy,
             "Star binding energy")
        .def_property_readonly("circ_radius", 
             &etk::spherical_star_properties::circ_radius,
             "Star circumferential radius")
        .def_property_readonly("proper_volume", 
             &etk::spherical_star_properties::proper_volume,
             "Star proper volume")
        .def_property_readonly("moment_inertia", 
             &etk::spherical_star_properties::moment_inertia,
             "Star moment of inertia")
        .def_property_readonly("has_bulk", 
             &etk::spherical_star_properties::has_bulk,
             "Whether bulk properties are available")
        .def_property_readonly("bulk", 
             &etk::spherical_star_properties::bulk,
             "Star bulk properties")
        .def_property_readonly("has_deform", 
             &etk::spherical_star_properties::has_deform,
             "Whether tidal deformability is available")
        .def_property_readonly("deformability", 
             &etk::spherical_star_properties::deformability,
             "Star tidal deformability and love number");
        
        
    py::class_<etk::spherical_star>(m, "spherical_star",
R"(Represents a spherical neutron star model

This provides everything spherical_star_properties does,
and in addition the radial profiles for metric and matter 
quantities. The profiles are provided as functions of circumferential
radius that can be evaluated both inside and outside the star.
All quantities are in the same geometric units used by the EOS.

)" )
        .def_property_readonly("eos", 
             &etk::spherical_star::eos,
             "Star EOS")
        .def_property_readonly("center_rho", 
             &etk::spherical_star::center_rho,
             "Central baryonic mass density")
        .def_property_readonly("center_eps", 
             &etk::spherical_star::center_eps,
             "Central specific internal energy")
        .def_property_readonly("center_press", 
             &etk::spherical_star::center_press,
             "Central pressure")
        .def_property_readonly("center_csnd", 
             &etk::spherical_star::center_csnd,
             "Central sound speed")
        .def_property_readonly("center_ye", 
             &etk::spherical_star::center_ye,
             "Central electron fraction")
        .def_property_readonly("grav_mass", 
             &etk::spherical_star::grav_mass,
             "Star gravitational mass")
        .def_property_readonly("bary_mass", 
             &etk::spherical_star::bary_mass,
             "Star baryonic mass")
        .def_property_readonly("binding_energy", 
             &etk::spherical_star::binding_energy,
             "Star binding energy")
        .def_property_readonly("circ_radius", 
             &etk::spherical_star::circ_radius,
             "Star circumferential radius")
        .def_property_readonly("proper_volume", 
             &etk::spherical_star::proper_volume,
             "Star proper volume")
        .def_property_readonly("moment_inertia", 
             &etk::spherical_star::moment_inertia,
             "Star moment of inertia")
        .def_property_readonly("has_bulk", 
             &etk::spherical_star::has_bulk,
             "Whether bulk properties are available")
        .def_property_readonly("bulk", 
             &etk::spherical_star::bulk,
             "Star bulk properties")
        .def_property_readonly("has_deform", 
             &etk::spherical_star::has_deform,
             "Whether tidal deformability is available")
        .def_property_readonly("deformability", 
             &etk::spherical_star::deformability,
             "Star tidal deformability and love number")
        .def("nu_from_rc", 
             py::vectorize(&etk::spherical_star::nu_from_rc),
             "Metric function nu at given circumf. radius",
             py::arg("rc"))
        .def("lambda_from_rc", 
             py::vectorize(&etk::spherical_star::lambda_from_rc),
             "Metric function lambda at circumf. radius",
             py::arg("rc"))     
        .def("gm1_from_rc", 
             py::vectorize(&etk::spherical_star::gm1_from_rc),
             "Pseudo enthalpy g - 1 at circumf. radius",
             py::arg("rc"))
        .def("mbary_from_rc", 
             py::vectorize(&etk::spherical_star::mbary_from_rc),
             "Baryonic mass within circumf. radius",
             py::arg("rc"))
        .def("pvol_from_rc", 
             py::vectorize(&etk::spherical_star::pvol_from_rc),
             "Proper volume within circumf. radius",
             py::arg("rc"))
        .def("rho_from_rc", 
             py::vectorize(&etk::spherical_star::rho_from_rc),
             "Baryonic mass density at circumf. radius",
             py::arg("rc"))
        .def("press_from_rc", 
             py::vectorize(&etk::spherical_star::press_from_rc),
             "Pressure at circumf. radius",
             py::arg("rc"))
        .def("eps_from_rc", 
             py::vectorize(&etk::spherical_star::eps_from_rc),
             "Specific internal energy at circumf. radius",
             py::arg("rc"))
        .def("csnd_from_rc", 
             py::vectorize(&etk::spherical_star::csnd_from_rc),
             "Soundspeed at circumf. radius",
             py::arg("rc"))
        .def("ye_from_rc", 
             py::vectorize(&etk::spherical_star::ye_from_rc),
             "Electron fraction at circumf. radius ",
             py::arg("rc"))
        .def("temp_from_rc", 
             py::vectorize(&etk::spherical_star::temp_from_rc),
             "Temperature at circumf. radius ",
             py::arg("rc"));

    py::class_<etk::tov_acc_simple>(m, "tov_acc_simple",
R"(Accuracy parameters for solving TOV and tidal deformability ODEs.


This contains heuristic parameters used in the adaptive ODE solvers.
Accuracy for TOV and tidal deformability ODEs are specified seperately.
There are no guarantees on the absolute errors of the various 
quantities, but lowering the values should lead to more accurate results.
    
)" )
        .def(py::init<real_t, real_t, std::size_t>(),
R"(
Args:
    tov (float): accuracy parameter for solving TOV ODE
    deform (float): accuracy parameter for solving deformability ODE
    minsteps (int): Minimum number of steps (controls maximum stepsize)
    
)",
        
             py::arg("tov")=1e-8, 
             py::arg("deform")=1e-6, 
             py::arg("minsteps")=500)
        .def_readonly("tov", &etk::tov_acc_simple::tov)
        .def_readonly("deform", &etk::tov_acc_simple::deform)
        .def_readonly("minsteps", &etk::tov_acc_simple::minsteps);

    py::class_<etk::tov_acc_precise>(m, "tov_acc_precise",
R"(Accuracy parameters for solving TOV and tidal deformability ODEs.

This allows to seperately specify the required tolerance for each 
quantity. This will repeatedly solve the ODEs with increasing accuracy
until the estimated residuals fall below the given tolerances.

)" )
        .def(py::init<real_t, real_t, real_t, real_t, 
                       std::size_t, real_t>(),
R"(Args:
    mass (float): Relative tolerance for baryonic and gravitational mass
    radius (float): Relative tolerance for radius and volume^(1/3)
    minertia (float): Relative tolerance for moment of inertia
    deform (float): Relative tolerance for tidal deformability (Lambda
        and k_2)
    minsteps (int): Minimum number of steps (controls maximum stepsize)
    acc_min (float): Give up if adaptive ODE solver tolerances get 
        any smaller.

)",
             py::arg("mass")=1e-8, 
             py::arg("radius")=1e-8, 
             py::arg("minertia")=1e-8, 
             py::arg("deform")=1e-6, 
             py::arg("minsteps")=500,
             py::arg("acc_min")=1e-14)
        .def_readonly("mass", &etk::tov_acc_precise::mass)
        .def_readonly("radius", &etk::tov_acc_precise::radius)
        .def_readonly("minertia", &etk::tov_acc_precise::minertia)
        .def_readonly("deform", &etk::tov_acc_precise::deform)
        .def_readonly("minsteps", &etk::tov_acc_precise::minsteps)
        .def_readonly("acc_min", &etk::tov_acc_precise::acc_min);

    py::class_<etk::star_seq>(m, "star_seq",
R"(Represents sequences of neutron stars or similar.

This class allows to store precomputed properties for a sequence 
of spherical stars and provide those as functions of the central 
pseudo-enthalpy. For this, regular spaced monotonic spline 
interpolation is used. The unit system can be chosen when creating
sequences, but is assumed to be geometric. It is stored for 
bookkeeping.

)")
        .def("grav_mass_from_center_gm1",  
             py::vectorize(&etk::star_seq::grav_mass_from_center_gm1),
             "Gravitational mass from central pseudo enthalpy",
             py::arg("gm1"))
        .def("bary_mass_from_center_gm1",  
             py::vectorize(&etk::star_seq::bary_mass_from_center_gm1),
             "Baryonic mass from central pseudo enthalpy",
             py::arg("gm1"))
        .def("circ_radius_from_center_gm1",  
             py::vectorize(&etk::star_seq::circ_radius_from_center_gm1),
             "Circumferential radius from central pseudo enthalpy",
             py::arg("gm1"))
        .def("moment_inertia_from_center_gm1",  
             py::vectorize(&etk::star_seq::moment_inertia_from_center_gm1),
             "Moment of inertia from central pseudo enthalpy",
             py::arg("gm1"))
        .def("lambda_tidal_from_center_gm1",  
             py::vectorize(&etk::star_seq::lambda_tidal_from_center_gm1),
             "Tidal deformability from central pseudo enthalpy",
             py::arg("gm1"))
        .def("contains_gm1",  
             py::vectorize(&etk::star_seq::contains_gm1),
             "If sequence contains a given central pseudo enthalpy",
             py::arg("gm1"))
        .def_property_readonly("range_center_gm1", 
             &etk::star_seq::range_center_gm1,
             "Range of central pseudo enthalpy g - 1")
        .def_property_readonly("units_to_SI", 
             &etk::star_seq::units_to_SI,
             "Unit system");

    py::class_<etk::star_branch, etk::star_seq>(m, "star_branch",
R"(Representing a stable branch of a star sequence.

This provides star properties as function of gravitational mass, in  
addition to the methods available from pyreprimand.star_seq. 

)")
        .def("grav_mass_from_center_gm1",  
             py::vectorize(
               &etk::star_branch::grav_mass_from_center_gm1),
             "Gravitational mass from central pseudo enthalpy",
             py::arg("gm1"))
        .def("bary_mass_from_center_gm1",  
             py::vectorize(
               &etk::star_branch::bary_mass_from_center_gm1),
             "Baryonic mass from central pseudo enthalpy",
             py::arg("gm1"))
        .def("circ_radius_from_center_gm1",  
             py::vectorize(
               &etk::star_branch::circ_radius_from_center_gm1),
             "Circumferential radius from central pseudo enthalpy",
             py::arg("gm1"))
        .def("moment_inertia_from_center_gm1",  
             py::vectorize(
               &etk::star_branch::moment_inertia_from_center_gm1),
             "Moment of inertia from central pseudo enthalpy",
             py::arg("gm1"))
        .def("lambda_tidal_from_center_gm1",  
             py::vectorize(
               &etk::star_branch::lambda_tidal_from_center_gm1),
             "Tidal deformability from central pseudo enthalpy",
             py::arg("gm1"))
        .def("center_gm1_from_grav_mass",  
             py::vectorize(
               &etk::star_branch::center_gm1_from_grav_mass),
             "Central pseudo enthalpy from grav. mass",
             py::arg("mg"))
        .def("bary_mass_from_grav_mass",  
             py::vectorize(
               &etk::star_branch::bary_mass_from_grav_mass),
             "Baryonic mass from grav. mass",
             py::arg("mg"))
        .def("circ_radius_from_grav_mass",  
             py::vectorize(
               &etk::star_branch::circ_radius_from_grav_mass),
             "Circumferential radius from grav. mass",
             py::arg("mg"))
        .def("moment_inertia_from_grav_mass",  
             py::vectorize(
               &etk::star_branch::moment_inertia_from_grav_mass),
             "Moment of inertia from grav. mass",
             py::arg("mg"))
        .def("lambda_tidal_from_grav_mass",  
             py::vectorize(
               &etk::star_branch::lambda_tidal_from_grav_mass),
             "Tidal deformability  from grav. mass",
             py::arg("mg"))
        .def_property_readonly("range_grav_mass", 
             &etk::star_branch::range_grav_mass,
             "Range of gravitational masses")
        .def("contains_grav_mass",  
             py::vectorize(
               &etk::star_branch::contains_grav_mass),
             "If sequence contains a given grav. mass.",
             py::arg("gm1"))
        .def_property_readonly("range_center_gm1", 
             &etk::star_branch::range_center_gm1,
             "Range of central pseudo enthalpy g - 1")
        .def("contains_gm1",  
             py::vectorize(&etk::star_branch::contains_gm1),
             "If sequence contains a given central pseudo enthalpy",
             py::arg("gm1"))
        .def_property_readonly("includes_maximum", 
             &etk::star_branch::includes_maximum,
             "If sequence includes maximum mass model")
        .def_property_readonly("grav_mass_maximum", 
             &etk::star_branch::grav_mass_maximum,
             "Maximum gravitational mass")
        .def_property_readonly("bary_mass_maximum", 
             &etk::star_branch::bary_mass_maximum,
             "Baryonic mass of model with maximum gravitational mass")
        .def_property_readonly("center_gm1_maximum", 
             &etk::star_branch::center_gm1_maximum,
             "Central pseudo enthalpy of model with maximum "
             "gravitational mass")
        .def_property_readonly("units_to_SI", 
             &etk::star_seq::units_to_SI,
             "Unit system");

            
    m.def("make_tov_star", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const etk::tov_acc_simple acc, const bool find_bulk,
             const bool find_tidal) 
          { 
             return etk::make_tov_star(eos, rho_center, acc, 
                                       find_bulk, find_tidal);
          }, 
R"(Compute a TOV solution 

Use this if you also need the radial profile, otherwise use
get_tov_star_properties. This function provides only basic control
on the accuracy via pyreprimand.tov_acc_simple. There is also a
version using pyreprimand.tov_acc_precise.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    acc (pyreprimand.tov_acc_simple): Tolerances for adaptive solver
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability

Returns:
    pyreprimand.spherical_star
)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("acc"),
          py::arg("find_bulk")=false,
          py::arg("find_tidal")=true);

    m.def("get_tov_star_properties", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const etk::tov_acc_simple acc, const bool find_bulk,
             const bool find_tidal) 
          { 
             return etk::get_tov_star_properties(eos, rho_center, 
                                          acc, find_bulk, find_tidal);
          }, 
R"(Compute a TOV solution 

Use this if you do not need the radial profile, otherwise use
make_tov_star. This function provides only basic control
on the accuracy via pyreprimand.tov_acc_simple. There is also a
version using pyreprimand.tov_acc_precise.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    acc (pyreprimand.tov_acc_simple): Tolerances for adaptive solver
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability

Returns:
    pyreprimand.spherical_star_properties

)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("acc"),
          py::arg("find_bulk")=false,
          py::arg("find_tidal")=true);

    m.def("make_tov_star", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const etk::tov_acc_precise acc, const bool find_bulk,
             const bool find_tidal) 
          { 
             return etk::make_tov_star(eos, rho_center, acc, 
                                       find_bulk, find_tidal);
          }, 
R"(Compute a TOV solution 

Use this if you also need the radial profile, otherwise use
get_tov_star_properties. This function provides strict control
on the accuracy via pyreprimand.tov_acc_precise. There is also a
faster version using pyreprimand.tov_acc_simple.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    acc (pyreprimand.tov_acc_simple): Tolerances for adaptive solver
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability

Returns:
    pyreprimand.spherical_star
)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("acc"),
          py::arg("find_bulk")=false,
          py::arg("find_tidal")=true);

    m.def("get_tov_star_properties", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const etk::tov_acc_precise acc, const bool find_bulk,
             const bool find_tidal) 
          { 
             return etk::get_tov_star_properties(eos, rho_center, 
                                      acc, find_bulk, find_tidal);
          }, 
R"(Compute a TOV solution 

Use this if you do not need the radial profile, otherwise use
make_tov_star. This function provides strict control
on the accuracy via pyreprimand.tov_acc_precise. There is also a
faster version using pyreprimand.tov_acc_simple.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    acc (pyreprimand.tov_acc_simple): Tolerances for adaptive solver
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability

Returns:
    pyreprimand.spherical_star_properties

)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("acc"),
          py::arg("find_bulk")=false,
          py::arg("find_tidal")=true);

    m.def("find_rhoc_tov_max_mass", 
          &etk::find_rhoc_tov_max_mass, 
          "Find maximum gravitational mass TOV model",
          py::arg("eos"),
          py::arg("rhobr0"),
          py::arg("rhobr1"),
          py::arg("nbits"),
          py::arg("acc"),
          py::arg("max_steps"));

    m.def("find_rhoc_tov_of_mass", 
          &etk::find_rhoc_tov_of_mass, 
          "Find TOV model with given gravitational mass",
          py::arg("eos"),
          py::arg("mg"),
          py::arg("rhobr0"),
          py::arg("rhobr1"),
          py::arg("acc"),
          py::arg("max_steps"));

    m.def("make_tov_seq", 
          &etk::make_tov_seq, 
R"(Compute sequence of TOV solutions

This computes the sequence of solutions within a given range of
central pseudo-enthalpy. If you only want to obtain the stable branch,
use make_tov_branch_stable instead.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the NSs. 
    acc (pyreprimand.tov_acc_simple): Tolerances for adaptive ODE solver.
    rg_gm1 (pyreprimand.range): Range of central pseudo enthalpy \f$ g-1 \f$
    num_samp (int): Number of sample points along the sequence.

Returns:
    pyreprimand.star_seq

)",
          py::arg("eos"),
          py::arg("acc"),
          py::arg("rg_gm1"),
          py::arg("num_samp")=500);

    m.def("make_star_seq", 
          &etk::make_star_seq, 
R"(Create star sequence from given arrays with star properties

This creates a star_seq object from given stellar properties for a 
seqquence. The provided sample points must be uniformly spaced in 
the central pseudo enthalpy g-1.

Args:
    mg (array_like): Samples for gravitational mass
    mb (array_like): Samples for baryonic mass
    rc (array_like): Samples for circumferential proper radius
    mi (array_like): Samples for moment of inertia
    lt (array_like): Samples for tidal deformability
    range_gm1 (pyreprimand.range): Range of central 
        pseudo-enthalpy g-1
    seq_units (pyreprimand.units): Unit system of the samples and
        the returned sequence object
    
Returns:
    pyreprimand.star_seq
)",
          py::arg("mg"),
          py::arg("mb"),
          py::arg("rc"),
          py::arg("mi"),
          py::arg("lt"),
          py::arg("range_gm1"),
          py::arg("seq_units"));

    m.def("load_star_seq", 
          &etk::load_star_seq, 
R"(Load sequence from file

Args:
    path (str): path of the sequence file.
    units (pyreprimand.units): Unit system to use for the sequence 
        object that is returned. This does *not* refer to the units
        inside the file, which are fixed.

Returns:
    pyreprimand.star_seq

)",
          py::arg("path"),
          py::arg("units"));

    m.def("save_star_seq", 
          &etk::save_star_seq, 
          "Save sequence to file",
          py::arg("path"),
          py::arg("seq"));

    m.def("load_star_branch", 
          &etk::load_star_branch, 
R"(Load sequence branch from file

Args:
    path (str): path of the branch file.
    units (pyreprimand.units): Unit system to use for the branch
        object that is returned. This does *not* refer to the units
        inside the file, which are fixed.

Returns:
    pyreprimand.star_branch

)",
          py::arg("path"),
          py::arg("units"));

    m.def("save_star_branch", 
          &etk::save_star_branch, 
          "Save sequence branch to file",
          py::arg("path"),
          py::arg("seq"));



    m.def("make_tov_branch_stable", 
          &etk::make_tov_branch_stable, 
R"(Compute stable branch of TOV solutions

This function employs heuristic algorithm to find the stable branch
of TOV solutions. Since there may be more than one such branch, one 
has to provide a central pseudo-enthalpy to indicate the correct one.
By increasing/decreasing this value successively by some factor,
a search interval is expanded until it brackets the maximum mass
or until the upper bound exceeds the EOS validity range. The initial 
guess is not required to be within a stable branch, and the default
value should work for any remotely realistic NS EOS.
In a similar fashion, a central pseudo-enthalpy for which the NS mass 
falls below the parameter mgrav_min is determined. Next, the maximum
is determined using a maximum search. However, the maximum might be 
located at the EOS validity bound. The maxmimum is considered physical 
based on a simple heuristics: (g_max-1)*(1+max_margin) < g_eos, where 
g_max and g_eos are the the central pseudo-enthalpy of the maximum mass 
model and the EOS upper validity bound. This criterion can later be 
queried using the includes_maximum() method of the branch object.
Finally, the branch is sampled with resolution given by num_samp.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the NSs. 
    acc (pyreprimand.tov_acc_simple) Tolerances for adaptive ODE solver.
    mgrav_min (float) Minimum gravitational mass that should be covered.
    num_samp (int) Sample resolution of the sequence.
    gm1_initial (float) Central enthalpy to indicate desired branch.
    max_margin Defines when maximum is considered physical.

Returns:
    pyreprimand.star_branch

)",
          py::arg("eos"),
          py::arg("acc"),
          py::arg("mgrav_min")=0.5,
          py::arg("num_samp")=500,
          py::arg("gm1_initial")=1.2,
          py::arg("max_margin")=1e-4);


}
