#ifndef EOS_BAROTROPIC_INTERNALS_H
#define EOS_BAROTROPIC_INTERNALS_H

#include "config.h"
#include <memory>
#include "eos_barotropic_impl.h"


namespace EOS_Toolkit {

namespace detail {

class eos_barotr_base {
  public:

  using impl_t = implementations::eos_barotr_impl;
  using spimpl_t = std::shared_ptr<const impl_t>;

  protected:

  explicit eos_barotr_base(spimpl_t eos) : pimpl(std::move(eos)) 
  {
    assert(pimpl); 
  }

  eos_barotr_base()                                  = default;
  eos_barotr_base(const eos_barotr_base &)           = default;
  eos_barotr_base(eos_barotr_base &&)                = default;
  eos_barotr_base& operator=(const eos_barotr_base&) = default;
  eos_barotr_base& operator=(eos_barotr_base&&)      = default;
  ~eos_barotr_base()                                 = default;

  auto impl() const -> const impl_t& {return *pimpl;}
  static auto ibad() -> const impl_t& {return *pbad;}

  class state_base {
    const impl_t& p { eos_barotr_base::ibad() };    
    const bool ok{false};
    
    protected:
      
    state_base() = default;

    explicit state_base(const impl_t& i) : p{i}, ok{true} {}
    
    state_base(const state_base&)             = default;
    state_base(state_base &&)                 = default;
    state_base& operator=(const state_base&)  = delete;
    state_base& operator=(state_base&&)       = delete;
    ~state_base()                             = default;

    auto am_ok() const -> bool {return ok;}
    auto impl() const -> const impl_t& {return p;}

  };


  private:

  const static spimpl_t pbad; ///< Dummy implementation that always throws
  spimpl_t pimpl{pbad};  ///< Pointer to implementation
};

} // namespace detail


} // namespace EOS_Toolkit 

#endif
