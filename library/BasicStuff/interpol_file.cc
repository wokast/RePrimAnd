#include "datastore.h"
#include "interpol.h"
#include "interpol_linear.h"
#include "interpol_pchip_spline.h"
#include "interpol_regspl.h"
#include "interpol_logspl.h"


namespace EOS_Toolkit {

auto make_interpolator(datasource s) 
->interpolator
{
  std::string t = s["interpolator_type"];

  if (t == detail::interpol_reglin_impl::datastore_id) 
  {
    return make_interpol_reglin(s);
  }

  if (t == detail::interpol_loglin_impl::datastore_id) 
  {
    return make_interpol_loglin(s);
  }

  if (t == detail::interpol_pchip_impl::datastore_id) 
  {
    return make_interpol_pchip_spline(s);
  }

  if (t == detail::interpol_regspl_impl::datastore_id) 
  {
    return make_interpol_regspl(s);
  }

  if (t == detail::interpol_logspl_impl::datastore_id) 
  {
    return make_interpol_logspl(s);
  }

  if (t == detail::interpol_llogspl_impl::datastore_id) 
  {
    return make_interpol_llogspl(s);
  }
  

  throw std::runtime_error("interpolator: encountered invalid " 
                             "datastore_id while reading");
}


} //namespace EOS_Toolkit 
