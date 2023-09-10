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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "FerretProblem.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "NonlinearSystem.h"

registerMooseObject("FerretApp", FerretProblem);

InputParameters
FerretProblem::validParams()
{
  InputParameters params = FEProblemBase::validParams();
  params.addRequiredParam<NonlinearVariableName>("polar_var", "Polar angle variable to be prescribed unique on the 0, pi interval.");
  params.addRequiredParam<NonlinearVariableName>("azimuth_phi_var","Azimuthal angle variable to be prescribed unique on the 0, 2 pi interval");
  return params;
}

FerretProblem::FerretProblem(const InputParameters & params)
  : FEProblem(params),
    // in the following have to get the names of the variables, and then find their numbers in
    // initialSetup,
    // as their numbers won't be defined at the moment of instantiation of this class
    _polar_var_name(params.get<NonlinearVariableName>("polar_var")),
    _azimuth_phi_var_name(params.get<NonlinearVariableName>("azimuth_phi_var")),
    _polar_var_num(0),
    _azimuth_phi_var_num(0)
{
}

FerretProblem::~FerretProblem() {}

void
FerretProblem::initialSetup()
{
  // the first argument to getVariable is threadID - i hope the following always works
  unsigned int tid = 0;

  // We are going to do more specific checks below, hence allowing
  // more specific error messages to be printed in case something goes
  // wrong. Therefore we just pass VAR_ANY here.
  MooseVariableFEBase & polar = getVariable(
      tid, _polar_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
  MooseVariableFEBase & azimuth = getVariable(
      tid, _azimuth_phi_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

  // some checks
  if (!polar.isNodal() || !azimuth.isNodal())
    mooseError("Both the polar_var and azimuth_phi_var variables must be nodal variables in FerretProblem");
  if (polar.feType().family != azimuth.feType().family)
    mooseError("Both the polar_var and azimuth_phi_var variables must belong to the same element family, eg LAGRANGE, in FerretProblem");
  if (polar.feType().order != azimuth.feType().order)
    mooseError("Both the polar_var and azimuth_phi_var variables must have the same order, eg FIRST, in FerretProblem");

  // extract the required info
  _polar_var_num = polar.number();
  _azimuth_phi_var_num = azimuth.number();

  FEProblemBase::initialSetup();
}

bool
FerretProblem::shouldUpdateSolution()
{
  return true;
}

bool
FerretProblem::updateSolution(NumericVector<Number> & vec_solution,
                                          NumericVector<Number> & ghosted_solution)
{
  bool updatedSolution =
      false; // this gets set to true if we needed to enforce the bound at any node

  unsigned int sys_num = getNonlinearSystemBase().number();

  // For parallel procs i believe that i have to use local_nodes_begin, rather than just nodes_begin
  // _mesh comes from SystemBase (_mesh = getNonlinearSystemBase().subproblem().mesh(), and
  // subproblem is this object)
  for (const auto & node : _mesh.getMesh().local_node_ptr_range())
  {
    // dofs[0] is the dof number of the polar variable at this node
    // dofs[1] is the dof number of the azimuth variable at this node
    std::vector<dof_id_type> dofs(2);
    dofs[0] = node->dof_number(sys_num, _polar_var_num, 0);
    dofs[1] = node->dof_number(sys_num, _azimuth_phi_var_num, 0);

    // soln[0] is the value of the polar variable at this node
    // soln[1] is the value of the azimuth variable at this node
    std::vector<Number> soln(2);
    vec_solution.get(dofs, soln);

    // do the bounding
    if (soln[0] > libMesh::pi)
    {
     // Moose::out << " dof: ";
     // std::cout << dofs[0];
      Moose::out << " ";
      Moose::out << " soln ";
      std::cout << soln[0];
      Moose::out << " ";
      vec_solution.set(dofs[0], soln[0] - libMesh::pi);
      updatedSolution = true;
    }
    else if (soln[0] < 0.0)
    {
      vec_solution.set(dofs[0], soln[0] + libMesh::pi);
      updatedSolution = true;
    }

    if (soln[1] > 2*libMesh::pi)
    {
      vec_solution.set(dofs[1], soln[1] - 2*libMesh::pi);
      updatedSolution = true;
    }
    else if (soln[1] < 0.0)
    {
      vec_solution.set(dofs[1], soln[1] + 2*libMesh::pi);
      updatedSolution = true;
    }
    //NOTE segfaults exist from block restriction
  }


  // The above vec_solution.set calls potentially added "set" commands to a queue
  // The following actions the queue (doing MPI commands if necessary), so
  // vec_solution will actually be modified by this following command
  vec_solution.close();

  // if any proc updated the solution, all procs will know about it
  _communicator.max(updatedSolution);

  if (updatedSolution)
  {
    ghosted_solution = vec_solution;
    ghosted_solution.close();
  }

  return updatedSolution;
}
