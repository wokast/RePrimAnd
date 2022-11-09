#ifndef STAR_SEQ_FILE_H
#define STAR_SEQ_FILE_H

#include "star_sequence.h"
#include "unitconv.h"
#include "datastore.h"
#include <string>


namespace EOS_Toolkit {

/**\brief Load NS sequence from file. 

@param fname File name of sequence file
@param Unit system the returned sequence should use. The unit system needs
to be geometric, i.e. \f$ G=c=1 \f$. Default is to fix the mass unit to
\f$ 1 M_\odot \f$. Note this does not refer to the units in the file, which is always in SI units.

@return Star sequence
**/ 
auto load_star_seq(std::string fname, 
                   const units& u=units::geom_solar()) 
->star_seq;


/**\brief Save NS sequence to file. 

@param fname File name
@param seq Sequence to save

This will save the sequence to file. Internally, the file
format uses SI units. Conversion is preformed using the units
stored in the sequence object.
**/ 
void save_star_seq(std::string fname, const star_seq& seq);

/**\brief Load NS branch sequence from file. 

@param fname Filename of branch sequence file
@param Unit system the returned sequence should use. 
@return Star sequence branch
**/ 
auto load_star_branch(std::string fname, 
                   const units& u=units::geom_solar()) 
->star_branch;

/**\brief Save NS branch to file. 

@param fname File name
@param seq Branch to save

This will save the branch to file. Internally, the file
format uses SI units. Conversion is preformed using the units
stored in the branch object.
**/ 
void save_star_branch(std::string fname, const star_branch& seq);

namespace detail {
auto load_star_seq(
  const datasource g,                 ///< Source to load from
  const units& u=units::geom_solar()  ///< Unit system to be used 
)
-> star_seq;


auto load_star_branch(
  const datasource g,                 ///< Source to load from
  const units& u=units::geom_solar()  ///< Unit system to be used 
) 
-> star_branch;


}

}

#endif

