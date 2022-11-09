#include "datastore.h"
#include "hdf5store.h"
#include "eos_thermal_file.h"
#include "eos_thermal_file_impl.h"
#include "eos_idealgas_impl.h"
#include "eos_hybrid_impl.h"

namespace EOS_Toolkit {


void ugly_hack_to_trick_stupid_linker2()
{ // This is only here to force the linker to include
  // _seemingly_ unused handler registration code 
  volatile bool builtin_handlers_registered {
    implementations::eos_idealgas::file_handler_registered &&
    implementations::eos_hybrid::file_handler_registered
  };
  assert(builtin_handlers_registered); 
}


eos_thermal detail::load_eos_thermal(const datasource s, 
                                    const units& u)
{
  ugly_hack_to_trick_stupid_linker2();
  
  if (s.has_data("eos_type"))  //old format
  {
    std::string spec = s["eos_type"];
    auto& r = implementations::registry_reader_eos_thermal::get(spec);
    return r.load(s / (std::string("eos_")+spec), u);  
  }
  auto g = s / "eos_thermal";
  std::string spec = g["eos_type"];
  auto& r = implementations::registry_reader_eos_thermal::get(spec);
  return r.load(g, u);  
}

eos_thermal load_eos_thermal(std::string fname, const units& u)
{
  auto g = make_hdf5_file_source(fname);
  return detail::load_eos_thermal(g,u);
}


void save_eos_thermal(std::string fname, eos_thermal eos, 
                     std::string info)
{
  auto g = make_hdf5_file_sink(fname);
  g["eos_info"] = info;
  eos.save(g / "eos_thermal");
}

} // namespace EOS_Toolkit
