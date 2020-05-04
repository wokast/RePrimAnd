#ifndef GLOBAL_REGISTRY_H
#define GLOBAL_REGISTRY_H

#include <string>
#include <unordered_map>
#include <cassert> 

namespace EOS_Toolkit {

/**\brief Global static registry of named polymorphic objects

This allows defining a map of names to objects subclassing an abstract
interface before main(). Used for example to register reader methods
for different EOS types in an extensible way.

\tparam F The base class of the objects that can be registered 
**/
template<class F> 
class global_registry {
  using map_t = std::unordered_map<std::string, const F*>;
  
  map_t reg;  
  
  static global_registry& single() 
  {
    static global_registry s;
    return s;
  }
  
  global_registry() = default;
  
  ~global_registry() {
     for (auto e : reg) delete e.second;
  }
  
  bool add_(std::string name, const F* f) 
  {
    assert(f != nullptr);
    if (reg.find(name) != reg.end()) {
      return false;
    }
    reg[name] = f;
    return true;
  }
  
  const F& get_(std::string name) 
  {
    auto i = reg.find(name);
    if (i == reg.end()) {
      std::string msg = std::string("GlobalRegistry: entry ") + name 
                         + std:: string(" not found") ;
      throw std::runtime_error(msg);
    }
    return *i->second;
  }
  
  public:
  
  static bool add(std::string name, const F* f)
  {
    return single().add_(name, f);
  }

  static const F& get(std::string name) 
  {
    return single().get_(name);
  }
};


} // namespace EOS_Toolkit 


#endif
