#include "hdf5imple.h"
#include "eos_thermal_file.h"
#include "eos_thermal_file_impl.h"

namespace EOS_Toolkit {


eos_thermal detail::load_eos_thermal(const h5grp& g, const units& u)
{
  std::string spec = get_attribute<std::string>(g, "eos_type");
  
  h5grp g2(g, std::string("eos_")+spec);

  auto& r = implementations::registry_reader_eos_thermal::get(spec);
  return r.load(g2,u);  
}

eos_thermal load_eos_thermal(std::string fname, const units& u)
{
  h5file f(fname);
  h5grp g(f, "/");
  return detail::load_eos_thermal(g,u);
}


} // namespace EOS_Toolkit
