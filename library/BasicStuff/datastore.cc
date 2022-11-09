#include "datastore.h"
#include <stdexcept>

namespace EOS_Toolkit {


bool datasource::has_group(std::string n) const 
{
  return pimpl->has_group(n);
}

datasource datasource::operator/(std::string n) const 
{
  return datasource(pimpl->group(n), pimpl);
}

bool datasource::has_data(std::string n) const 
{
  return pimpl->has_data(n);
}

detail::source_proxy datasource::operator[](std::string n) const 
{
  return detail::source_proxy(*this, n);
}

datasink datasink::operator/(std::string n)
{
  return datasink(pimpl->group(n), pimpl);
}

detail::sink_proxy datasink::operator[](std::string n) 
{
  return detail::sink_proxy(*this, n);
}


} //namespace EOS_Toolkit 
