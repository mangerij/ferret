/**
 * @file   ElectrostrictiveCoefficientR4.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Nov 25 10:32:53 2013
 * @brief  ElectrostrictiveCoefficient is a rank four tensor; the entries are indexed by m,n,k,l equal to 0,1,2
 *         It holds 81 entries, however, only 54 entries are independent due to the constrait symmetry constraint Q(m,n,k,l)=Q(m,n,l,k)
 */
#ifndef ELECTROSTRICTIVECOEFFICIENTR4_H
#define ELECTROSTRICTIVECOEFFICIENTR4_H

#include "MooseError.h"
#include "RankFourTensor.h"
class ElectrostrictiveCoefficientR4: public RankFourTensor
{
public:
  virtual ~ElectrostrictiveCoefficientR4(){}
  virtual void fillFromInputVector(const std::vector<Real>& input, std::string type);
};

#endif //ELECTROSTRICTIVECOEFFICIENTR4_H
