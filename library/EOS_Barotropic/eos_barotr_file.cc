#include "hdf5imple.h"
#include "eos_barotr_file.h"
#include "eos_barotr_file_impl.h"

namespace EOS_Toolkit {


eos_barotr detail::load_eos_barotr(const h5grp& g, const units& u)
{
  std::string spec = get_attribute<std::string>(g, "eos_type");
  
  h5grp g2(g, std::string("eos_")+spec);

  auto& r = implementations::registry_reader_eos_barotr::get(spec);
  return r.load(g2,u);  
}

eos_barotr load_eos_barotr(std::string fname, const units& u)
{
  h5file f(fname);
  h5grp g(f, "/");
  return detail::load_eos_barotr(g,u);
}


} // namespace EOS_Toolkit
