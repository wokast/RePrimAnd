#include "datastore.h"
#include "hdf5store.h"
#include "eos_barotr_file.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_gpoly_impl.h"
#include "eos_barotr_pwpoly_impl.h"
#include "eos_barotr_poly_impl.h"
#include "eos_barotr_spline_impl.h"
#include "eos_barotr_table_impl.h"


namespace EOS_Toolkit {

void ugly_hack_to_trick_stupid_linker()
{ // This is only here to force the linker to include
  // _seemingly_ unused handler registration code 
  volatile bool builtin_handlers_registered {
    implementations::eos_barotr_gpoly::file_handler_registered &&
    implementations::eos_barotr_pwpoly::file_handler_registered &&
    implementations::eos_barotr_poly::file_handler_registered &&
    implementations::eos_barotr_spline::file_handler_registered &&
    implementations::eos_barotr_table::file_handler_registered
  };
  assert(builtin_handlers_registered); 
}

eos_barotr detail::load_eos_barotr(const datasource s, const units& u)
{
  ugly_hack_to_trick_stupid_linker();

  if (s.has_data("eos_type"))  //old format
  {
    std::string spec = s["eos_type"];
    auto& r = implementations::registry_reader_eos_barotr::get(spec);
    return r.load(s / (std::string("eos_")+spec), u);  
  }
  auto g = s / "eos_barotropic";
  std::string spec = g["eos_type"];
  auto& r = implementations::registry_reader_eos_barotr::get(spec);
  return r.load(g, u);  
}

void detail::save_eos_barotr(datasink g, eos_barotr eos)
{
  eos.save(g / "eos_barotropic");
}

eos_barotr load_eos_barotr(std::string fname, const units& u)
{
  auto g = make_hdf5_file_source(fname);
  return detail::load_eos_barotr(g,u);
}

void save_eos_barotr(std::string fname, eos_barotr eos, 
                     std::string info)
{
  auto g = make_hdf5_file_sink(fname);
  g["eos_info"] = info;
  save_eos_barotr(g, eos);
}

} // namespace EOS_Toolkit
