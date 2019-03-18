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
#include "RandomData.h"
#include "MooseRandom.h"

// Forward Declarations
class CorrelatedRandomFieldIC;
namespace libMesh 
{ 
class Point; 
}

template<>
InputParameters validParams<CorrelatedRandomFieldIC>();

class CorrelatedRandomFieldIC : public InitialCondition
{
public:
  CorrelatedRandomFieldIC(const InputParameters & parameters);

  void initialSetup() override;
  virtual std::vector<std::vector<std::vector<std::vector<Real>>>> fourierCoeffs();

  virtual Real evaluate(const Point & p, std::vector<std::vector<std::vector<std::vector<Real>>>> output);
  virtual Real evaluateNoArray(const Point & p);
  virtual Real value(const Point & p);

protected:
  /**
   * Generate a uniformly distributed random number on the interval from 0 to 1
   * @return random number
   */
  Real generateRandom();

  /// Determines whether a variable basis is elemental or nodal
  const bool _is_nodal;

  /// Boolean to indicate whether we want to use the old (deprecated) generation pattern
  const bool _use_legacy;
private:
  /// RandomData element object, we cannot inherit from RandomInterface in an InitialCondition
  std::unique_ptr<RandomData> _elem_random_data;

  /// RandomData node object, we cannot inherit from RandomInterface in an InitialCondition
  std::unique_ptr<RandomData> _node_random_data;

  /// Elemental random number generator
  MooseRandom * _elem_random_generator;

  /// Nodal random number generator
  MooseRandom * _node_random_generator;

  /// Random numbers per element (currently limited to a single value at a time)
  std::map<dof_id_type, Real> _elem_numbers;

  /// Random numbers per node (currently limited to a single value at a time)
  std::map<dof_id_type, Real> _node_numbers;

////////
// User defined things:
//
////////
  //correlation length scale
  const Real _Lcorr;
  /// The dimension of the mesh
  MooseEnum _dim;
  /// The min/max values for x,y,z component
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
  //
  Real _Nnodes;
  //mesh object
  const MeshBase & _mesh;
};

#endif //CORRELATEDRANDOMFIELDIC_H
