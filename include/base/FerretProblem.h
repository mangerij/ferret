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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef FERRETPROBLEM_H
#define FERRETPROBLEM_H

#include "FEProblem.h"

class FerretProblem : public FEProblem
{
public:
  FerretProblem(const InputParameters & params);
  static InputParameters validParams();

  virtual ~FerretProblem();
  virtual void initialSetup();
  virtual bool shouldUpdateSolution();

  /**
   * Does the bounding by modifying vec_solution, and then ghosted_solution
   * @param vec_solution is the solution that Petsc says we should use.
   * @param ghosted_solution is a ghosted version of vec_solution.
   * @return true if vec_solution was changed at a node in order to respect the bounds
   */
  virtual bool updateSolution(NumericVector<Number> & vec_solution,
                              NumericVector<Number> & ghosted_solution);

protected:
  NonlinearVariableName _polar_var_name;
  NonlinearVariableName _azimuth_phi_var_name;
  unsigned int _polar_var_num;
  unsigned int _azimuth_phi_var_num;
};

#endif /* FERRETPROBLEM_H */
