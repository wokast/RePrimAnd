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
        .def_readonly_static("FORMAL_BARYON_MASS_SI", 
                             &etk::units::FORMAL_BARYON_MASS_SI,
R"(Convention for the mass constant defining baryonic mass density
in terms of baryon number density, in SI units.)")
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

    py::class_<etk::spherical_star_bulk>(m, "spherical_star_bulk",
    "Describes bulk properties" )
        .def_readonly("circ_radius", 
                      &etk::spherical_star_bulk::circ_radius,
                      "Bulk circumferential radius")
        .def_readonly("rho", 
                      &etk::spherical_star_bulk::rho,
                      "Bulk baryonic mass density")
        .def_readonly("proper_volume", 
                      &etk::spherical_star_bulk::proper_volume,
                      "Bulk enclosed proper volume")
        .def_readonly("bary_mass", 
                      &etk::spherical_star_bulk::bary_mass,
                      "Bulk enclosed baryonic mass");

    
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
        .def_property_readonly("center_gm1", 
             &etk::spherical_star_properties::center_gm1,
             "Central pseudo enthalpy")
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



    py::class_<etk::star_accuracy_spec>(m, "star_accuracy_spec",
R"(Accuracy requirements for NS model properties.

This allows to seperately specify the desired tolerance for each 
quantity, and whether the optional variables tidal deformability
and bulk radius are required.

)" )
        .def_readonly("acc_mass", &etk::star_accuracy_spec::acc_mass)
        .def_readonly("acc_radius", &etk::star_accuracy_spec::acc_radius)
        .def_readonly("acc_minertia", &etk::star_accuracy_spec::acc_minertia)
        .def_readonly("minsteps", &etk::star_accuracy_spec::minsteps)
        .def_readonly("need_deform", &etk::star_accuracy_spec::need_deform)
        .def_readonly("acc_deform", &etk::star_accuracy_spec::acc_deform)
        .def_readonly("need_bulk", &etk::star_accuracy_spec::need_bulk);

    m.def("star_acc_simple", &etk::star_acc_simple,
R"(Simplified accuracy specification for NS models.

Args:
    need_deform (boolean): If tidel deformability is needed
    need_bulk (bool): If bulk radius is needed
    acc_tov (float): Relative error tolerance for masses, radii, 
                     and  moment of inertia.
    acc_deform (float): Relative error tolerance for tidal deformability
    minsteps (int): Minimum resolution inside NS.
    
Returns:
    pyreprimand.star_accuracy_spec
)",
        py::kw_only(),
        py::arg("need_deform")=true, 
        py::arg("need_bulk")=false, 
        py::arg("acc_tov")=1e-6,
        py::arg("acc_deform")=1e-3, 
        py::arg("minsteps")=20);



    m.def("star_acc_detailed", &etk::star_acc_detailed,
R"(Detailed accuracy specification for NS models.

Args:
    need_deform (boolean): If tidel deformability is needed
    need_bulk (bool): If bulk radius is needed
    acc_mass (float): Relative error tolerance for mass
    acc_radius (float): Relative error tolerance for radius
    acc_minertia (float): Relative error tolerance for moment of inertia
    acc_deform (float): Relative error tolerance for tidal deformability
    minsteps (int): Minimum resolution inside NS.
    
Returns:
    pyreprimand.star_accuracy_spec
)",
          py::kw_only(),
          py::arg("need_deform")=true,
          py::arg("need_bulk")=false,
          py::arg("acc_mass")=1e-6, 
          py::arg("acc_radius")=1e-6, 
          py::arg("acc_minertia")=1e-5,
          py::arg("acc_deform")=1e-3, 
          py::arg("minsteps")=20);


    m.def("get_tov_properties", &etk::get_tov_properties,
R"(Compute properties of TOV solution for given EOS and central density.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density
    acc (pyreprimand.star_accuracy_spec): Specifies desired accuracy and
      which optional quantities are needed. To set up the accuracy
      spec, use star_acc_detailed() or star_acc_simple(). 
      
Returns:
    pyreprimand.spherical_star_properties
)",    
    
    py::arg("eos"), 
    py::arg("rho_center"), 
    py::arg("acc"));


    m.def("get_tov_properties", 
          [] (const etk::eos_barotr eos, 
              const etk::real_t rho_center) 
          { 
            return etk::get_tov_properties(eos, rho_center, 
                                           etk::star_acc_simple());
          },
R"(Compute properties of TOV solution for given EOS and central density.

This variant uses a default accuracy specification equivalent to
get_tov_properties(eos, rho_center, star_acc_simple())
 
Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density
      
Returns:
    pyreprimand.spherical_star_properties
)",    
    
    py::arg("eos"), 
    py::arg("rho_center"));


 

    m.def("get_tov_star", &etk::get_tov_star,
R"(Compute TOV solution for given EOS and central density.

Use this if you also need the radial profile, otherwise use
get_tov_properties().

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density
    acc (pyreprimand.star_accuracy_spec): Specifies desired accuracy and
      which optional quantities are needed. To set up the accuracy
      spec, use star_acc_detailed() or star_acc_simple(). 
      
Returns:
    pyreprimand.spherical_star
)",    
    
    py::arg("eos"), 
    py::arg("rho_center"), 
    py::arg("acc"));


    m.def("get_tov_star", 
          [] (const etk::eos_barotr eos, 
              const etk::real_t rho_center) 
          { 
            return etk::get_tov_star(eos, rho_center, 
                                           etk::star_acc_simple());
          },
R"(Compute TOV solution for given EOS and central density.

Use this if you also need the radial profile, otherwise use
get_tov_properties().
This variant uses a default accuracy specification equivalent to
get_tov_properties(eos, rho_center, star_acc_simple())
 
Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density
      
Returns:
    pyreprimand.spherical_star
)",    
    
    py::arg("eos"), 
    py::arg("rho_center"));


    m.def("get_tov_properties_fixstep", 
          [] (const etk::eos_barotr eos, const real_t rho_center, 
              const etk::star_accuracy_spec acc) 
          {
            return etk::get_tov_properties_fixstep(eos, 
                                                   rho_center, acc); 
          },  
R"(Compute properties of TOV solution for given EOS and central density.

This function is intended for development and debugging.
It uses an implementation based on fixed-stepsize ODE solvers. 
For use in production, use get_tov_properties() instead. 

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density
    acc (pyreprimand.star_accuracy_spec): Specifies desired accuracy and
      which optional quantities are needed. To set up the accuracy
      spec, use star_acc_detailed() or star_acc_simple(). 
      
Returns:
    pyreprimand.spherical_star_properties
)",    
    
    py::arg("eos"), 
    py::arg("rho_center"), 
    py::arg("acc"));
 

    m.def("get_tov_properties_fixstep", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const bool find_bulk,
             const bool find_tidal,
             const std::size_t nsamp_tov, 
             const std::size_t nsub_tidal,
             const real_t wdiv_tidal, 
             const real_t bulk_acc) 
          { 
             return etk::get_tov_properties_fixstep(eos, rho_center, 
                                find_bulk, find_tidal,
                                nsamp_tov, nsub_tidal, wdiv_tidal, 
                                bulk_acc);
          }, 
R"(Compute a TOV solution with fixed ode step size

This low-level function is only meant for testing and development.
It might change or disappear at any new version.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    nsamp_tov (int): In how many steps to divide TOV ODE integration interval.
    nsub_tidal (int): Factor to divide TOV ODE step size for tidal ODE.
    wdiv_tidal (float): where to switch between the two tidal ODE variants
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability
    bulk_acc (float): Root finding accuracy for bulk radius

Returns:
    pyreprimand.spherical_star_properties

)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("find_bulk"),
          py::arg("find_tidal"),
          py::arg("nsamp_tov"),
          py::arg("nsub_tidal")=2,
          py::arg("wdiv_tidal")=0.91,
          py::arg("bulk_acc")=1e-8);

    m.def("get_tov_properties_adaptive", 
          [](etk::eos_barotr eos, const real_t rho_center, 
             const std::size_t nsamp_tov, 
             const real_t acc_tov,
             const std::size_t nsamp_tidal, 
             const real_t acc_tidal,
             const real_t wdiv_tidal, 
             const bool find_bulk,
             const bool find_tidal, const real_t bulk_acc) 
          { 
             return etk::get_tov_properties_adaptive(eos, rho_center, 
                                nsamp_tov, acc_tov, nsamp_tidal, 
                                acc_tidal, wdiv_tidal, 
                                find_bulk, find_tidal, bulk_acc);
          }, 
R"(Compute a TOV solution with adaptive ode step size

This low-level function is only meant for testing and development.
It might change or disappear at any new version.

Args:
    eos (pyreprimand.eos_barotr): The EOS of the star
    rho_center (float): The central baryonic mass density [EOS Units].
    nsamp_tov (int): Sample resolution for TOV profile
    acc_tov (int): Adaptive accuracy parameter for TOV ODE
    nsamp_tidal (int): Minimum number of steps for tidal ODEs
    acc_tidal (float): Adaptive accuracy parameter for tidal ODE 
    wdiv_tidal (float): where to switch between the two tidal ODE variants
    find_bulk (bool): Whether to also compute the "bulk" properties
    find_tidal (bool): Whether to compute the tidal deformability
    bulk_acc (float): Root finding accuracy for bulk radius

Returns:
    pyreprimand.spherical_star_properties

)",
          py::arg("eos"),
          py::arg("rho_center"),
          py::arg("nsamp_tov"),
          py::arg("acc_tov"),
          py::arg("nsamp_tidal"),
          py::arg("acc_tidal"),
          py::arg("wdiv_tidal"),
          py::arg("find_bulk"),
          py::arg("find_tidal"),
          py::arg("bulk_acc"));


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
    acc (pyreprimand.star_accuracy_spec): Specifies desired accuracy and
      which optional quantities are needed. To set up the accuracy
      spec, use star_acc_detailed() or star_acc_simple(). 
    rg_gm1 (pyreprimand.range): Range of central pseudo enthalpy \f$ g-1 \f$
    num_samp (int): Number of sample points along the sequence.

Returns:
    pyreprimand.star_seq

)",
          py::arg("eos"),
          py::arg("rg_gm1"),
          py::arg("acc"),
          py::arg("num_samp")=500);

    m.def("make_star_seq", 
          [] (std::vector<real_t> mg, 
              std::vector<real_t> mb, 
              std::vector<real_t> rc, 
              std::vector<real_t> mi, 
              std::vector<real_t> lt, 
              etk::star_seq::range_t rg_gm1, 
              etk::units u) 
          {
            return etk::make_star_seq(mg,mb,rc,mi,lt,rg_gm1, u);
          }, 
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
If the given value is on a stable branch, that branch will be returned.
If the value lies on an unstable branch, the adjacent stable branch
at lower densities will be returned. The default value should work for 
any remotely realistic NS EOS.

For most astrophysical applications, the low-mass part of a branch 
is not relevant. One can specify a low-mass cutoff factor relative 
to the maximum mass to indicate that lower masses are not needed.
The default cutoff is 1/5 of the maximum mass. One can also specify 
an additional absolute mass cutoff, but beware that an exception 
is thrown if it happens to exceed the maximum mass. If the mass cutoff
is below the actual minimum mass of the stable branch, the branch
will extend down to the minimum mass. Set the cutoffs to zero if the
goal is to find the minimum mass.

The parameter for the accuracy refers to the accuracy of the TOV
solutions, but does not include the interpolation error of the 
resulting star sequence. The algorithm samples the TOV sequence 
using regular steps in log(g-1), g being the central pseudo-enthalpy. 
The sample point of maximum mass will be moved closer to the true
maximum using a quadratic approximation around the maximum. 

The stepsize can be controled by a parameter. The default stepsize is 
chosen such that the total error is dominated by the TOV solution 
error for the default TOV solver accurracy. Changing the step size 
should only be required when extreme accuracy is required, or when 
accuracy is less important than speed. The computational costs are 
roughly antiproportional to the stepsize parameter.

For some EOS, the maximum NS mass is limited by the 
EOS validity range and not by the physical maximum of TOV solutions. 
This case can be queried using the includes_maximum() method of the 
returned branch object. For this, the maxmimum is considered physical
based on a simple heuristics: (g_max-1)*(1+max_margin) < g_eos, where 
g_max and g_eos are the the central pseudo-enthalpy of the maximum mass 
model and the EOS upper validity bound. 

Args:
    eos (pyreprimand.eos_barotr): The EOS of the NSs. 
    acc (pyreprimand.star_accuracy_spec): Specifies desired accuracy and
      which optional quantities are needed. To set up the accuracy
      spec, use star_acc_detailed() or star_acc_simple().
    mg_cut_low_rel (float) Low-mass cutoff in terms of maximum mass.
    mg_cut_low_abs (float) Low-mass cutoff in terms of absolute mass.
    gm1_initial (float) Central enthalpy g-1 to indicate desired branch.
    gm1_step (int) Sample resolution in terms of (delta g) / (g-1)
    max_margin Distance to EOS validity bound needed to consider maximum as physical.

Returns:
    pyreprimand.star_branch

)",
          py::arg("eos"),
          py::arg("acc"),
          py::kw_only(),
          py::arg("mg_cut_low_rel")=0.2,
          py::arg("mg_cut_low_abs")=0.0,
          py::arg("gm1_initial")=1.2,
          py::arg("gm1_step")=0.004,
          py::arg("max_margin")=1e-2);


    m.def("k2_from_ym2_mbr_stable", 
          py::vectorize(&etk::k2_from_ym2_mbr_stable),
R"(Numerically stable implementation of formula for computing 
   love number k2
  
   This computes k2 from y-2 and beta, where y is the surface value
   of the tidal ODE, and beta=M/R is the NS compactness. For small
   compactness beta, a Taylor expansion has to be used because the 
   exact formula suffers from extreme cancellation errors. The 
   coefficient of the highest order term beta^6 is not computed 
   analytically but from the difference between the terms of lower 
   order and the exact formula, evaluated at a value were the latter 
   is still more accurate. This is reducing the overall error 
   further. 
   
   This function is intended for development and debugging. 
)",
          py::arg("ym2"),
          py::arg("beta"),
          py::arg("b_thresh")=5e-2);

}
