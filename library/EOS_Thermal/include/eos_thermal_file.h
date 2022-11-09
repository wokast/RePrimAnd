#ifndef EOS_THERMAL_FILE_H
#define EOS_THERMAL_FILE_H

#include "eos_thermal.h"
#include "unitconv.h"
#include "datastore.h"
#include <string>


namespace EOS_Toolkit {


/**\brief Load thermal EOS from hdf5 file. 

@param fname Filename of EOS file
@param Unit system the returned EOS should use. The unit system needs
to be geometric, i.e. \f$ G=c=1 \f$. Default is to fix the mass unit to
\f$ 1 M_\odot \f$.

@return Generic interface to the EOS  
**/ 
eos_thermal load_eos_thermal(std::string fname, 
                             const units& u=units::geom_solar());

/**\brief Save thermal EOS to hdf5 file. 

@param fname Filename of EOS file
@param eos EOS object to save
@param info Free-form description text (optional)
**/ 
void save_eos_thermal(std::string fname,  eos_thermal eos, 
                     std::string info="");


namespace detail {

eos_thermal load_eos_thermal(
  const datasource g,                 ///< Source to load from
  const units& u=units::geom_solar()  ///< Unit system to be used by EOS
);

}

}

#endif

