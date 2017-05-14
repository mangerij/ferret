/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "TensorMechanicsActionScaled.h"

#include "Factory.h"
#include "FEProblem.h"
#include "Parser.h"

template<>
InputParameters validParams<TensorMechanicsActionScaled>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<NonlinearVariableName>("disp_x", "The x displacement");
  params.addParam<NonlinearVariableName>("disp_y", "The y displacement");
  params.addParam<NonlinearVariableName>("disp_z", "The z displacement");
  params.addParam<NonlinearVariableName>("temp", "The temperature");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<bool>("use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");

  params.addParam<std::vector<AuxVariableName> >("save_in_disp_x", "Auxiliary variables to save the x displacement residuals.");
  params.addParam<std::vector<AuxVariableName> >("save_in_disp_y", "Auxiliary variables to save the y displacement residuals.");
  params.addParam<std::vector<AuxVariableName> >("save_in_disp_z", "Auxiliary variables to save the z displacement residuals.");

  return params;
}

TensorMechanicsActionScaled::TensorMechanicsActionScaled(const InputParameters & parameters) :
    Action(parameters)
{
}

void
TensorMechanicsActionScaled::act()
{
  unsigned int dim = 1;
  std::vector<std::string> keys;
  std::vector<VariableName> vars;
  std::string type("StressDivergenceTensorsScaled");

  std::vector<std::vector<AuxVariableName> > save_in;
  ///Prepare displacements and set value for dim
  keys.push_back("disp_x");
  vars.push_back(getParam<NonlinearVariableName>("disp_x"));

  if (isParamValid("disp_y"))
  {
    ++dim;
    keys.push_back("disp_y");
    vars.push_back(getParam<NonlinearVariableName>("disp_y"));
    if (isParamValid("disp_z"))
    {
      ++dim;
      keys.push_back("disp_z");
      vars.push_back(getParam<NonlinearVariableName>("disp_z"));
    }
  }

  save_in.resize(dim);

  if (isParamValid("save_in_disp_x"))
    save_in[0] = getParam<std::vector<AuxVariableName> >("save_in_disp_x");

  if (isParamValid("save_in_disp_y"))
    save_in[1] = getParam<std::vector<AuxVariableName> >("save_in_disp_y");

  if (isParamValid("save_in_disp_z"))
    save_in[2] = getParam<std::vector<AuxVariableName> >("save_in_disp_z");

  ///Add in the temperature
  unsigned int num_coupled(dim);
  if (isParamValid("temp"))
  {
    ++num_coupled;
    keys.push_back("temp");
    vars.push_back(getParam<NonlinearVariableName>("temp"));
  }

  InputParameters params = _factory.getValidParams(type);
  for (unsigned int j = 0; j < num_coupled; ++j)
  {
    params.addCoupledVar(keys[j], "");
    params.set<std::vector<VariableName> >(keys[j]) = std::vector<VariableName>(1, vars[j]);
  }

  params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

  if (isParamValid("base_name"))
    params.set<std::string>("base_name") = getParam<std::string>("base_name");

  std::string short_name = "TensorMechanicsScaled";

  for (unsigned int i = 0; i < dim; ++i)
  {
    std::stringstream name;
    name << short_name;
    name << i;

    params.set<unsigned int>("component") = i;
    params.set<NonlinearVariableName>("variable") = vars[i];
    params.set<std::vector<AuxVariableName> >("save_in") = save_in[i];

    _problem->addKernel(type, name.str(), params);
  }
}
