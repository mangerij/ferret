/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "LandauMaterialBulkFourthDerivative.h"
#include "LandauTensorTools.h"
#include "ComputeRankTwoLandauTensor.h"
#include "ComputeRankFourLandauTensor.h"

registerMooseObject("FerretApp", LandauMaterialBulkFourthDerivative);

InputParameters LandauMaterialBulkFourthDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy which is an fourth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

LandauMaterialBulkFourthDerivative::LandauMaterialBulkFourthDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _rank_two_landau_tensor(getMaterialProperty<RankTwoTensor>("rank_two_landau_tensor")),
   _rank_four_landau_tensor(getMaterialProperty<RankFourTensor>("rank_four_landau_tensor")),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
LandauMaterialBulkFourthDerivative::computeQpResidual()
{
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (_component == 0)
  {
    return LandauTensorTools::landauTwoProductDerivative(_rank_two_landau_tensor[_qp], 0, p) + LandauTensorTools::landauFourProductDerivative(_rank_four_landau_tensor[_qp], 0, p);
  }
  else if (_component == 1)
  {
    return LandauTensorTools::landauTwoProductDerivative(_rank_two_landau_tensor[_qp], 1, p) + LandauTensorTools::landauFourProductDerivative(_rank_four_landau_tensor[_qp], 1, p);
  }
  else if (_component == 2)
  {
    return LandauTensorTools::landauTwoProductDerivative(_rank_two_landau_tensor[_qp], 2, p) + LandauTensorTools::landauFourProductDerivative(_rank_four_landau_tensor[_qp], 2, p);
  }
  else 
  {
    return 0.0;
  }
}

Real
LandauMaterialBulkFourthDerivative::computeQpJacobian()
{  
  if (_component == 0)
  {
    return 0.0;
  }
  else if (_component == 1)
  {
    return 0.0;
  }
  else if (_component == 2)
  {
    return 0.0;
  }
  else 
  {
    return 0.0;
  }
}


Real
LandauMaterialBulkFourthDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if(_component == 0) {
    if (jvar == _polar_y_var)
    {
        return 0.0;
    }
    else if (jvar == _polar_z_var)
    {
        return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else if(_component == 1) {
    if(jvar == _polar_x_var) {
        return 0.0;
    } 
    else if(jvar == _polar_z_var) {
        return 0.0;
    }
    else 
    {
        return 0.0;
    }
  }
  else if (_component == 2) {
    if(jvar == _polar_x_var) {
        return 0.0;
    } 
    else if(jvar == _polar_y_var) {
        return 0.0;
    }
    else
    {
        return 0.0;
    }
  }
  else 
  {
    return 0.0;
  }
}
