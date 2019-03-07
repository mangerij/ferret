/**
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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef CORRELATEDRANDOMFIELDIC_H
#define CORRELATEDRANDOMFIELDIC_H

#include "InitialCondition.h"

// Forward Declarations
class CorrelatedRandomFieldIC;
namespace libMesh { class Point; }

template<>
InputParameters validParams<CorrelatedRandomFieldIC>();

class CorrelatedRandomFieldIC : public InitialCondition
{
public:
  CorrelatedRandomFieldIC(const InputParameters & parameters);

  virtual std::vector<std::vector<std::vector<std::vector<Real>>>> fourierCoeffs();

  virtual Real evaluate(const Point & p, std::vector<std::vector<std::vector<std::vector<Real>>>> output);
  virtual Real evaluateNoArray(const Point & p);
  virtual Real value(const Point & p);

protected:

private:
  //correlation length scale
  const Real _Lcorr;
  /// The dimension of the mesh
  MooseEnum _dim;
  /// The min/max values for x,y,z component
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
  //mesh object
  const MeshBase & _mesh;
};

#endif //CORRELATEDRANDOMFIELDIC_H
