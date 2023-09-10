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

#include "ABO3CoupledPhaseFieldAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("FerretApp", ABO3CoupledPhaseFieldAction, "add_kernel");

//registerMooseAction("FerretApp", ABO3CoupledPhaseFieldAction, "add_variable");
//registerMooseAction("FerretApp", ABO3CoupledPhaseFieldAction, "add_aux_variable");

registerMooseAction("FerretApp", ABO3CoupledPhaseFieldAction, "add_material");

registerMooseAction("FerretApp", ABO3CoupledPhaseFieldAction, "add_postprocessor");

InputParameters
ABO3CoupledPhaseFieldAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up homogeneous or inhomogeneous ferroelectric materials problem. This can be steady-state or time-dependent. The modes of coupling can be to the inhomogeneous strain fields or renormalized in the potential.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "variables", "Polarization vector components, electrostatic potential, and elastic displacements components");
  /*MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  params.addParam<MooseEnum>("family",
                             families,
                             "Specifies the family of FE "
                             "shape function to use for the order parameters");
  params.addParam<MooseEnum>("order",
                             orders,
                             "Specifies the order of the FE "
                             "shape function to use for the order parameters");*/
  params.addRequiredParam<std::vector<std::string>>(
      "alpha_ijkl", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("alpha_ijkl_val", "Landau expansion coefficients");
  params.addRequiredParam<std::vector<std::string>>(
      "G_ij", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("G_ij_val", "gradient coefficients");
  params.addRequiredParam<std::vector<std::string>>(
      "Q_ij", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("Q_ij_val", "Electrostrictive Q tensor coefficients");
  params.addRequiredParam<std::vector<std::string>>(
      "q_ij", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("q_ij_val", "Electrostrictive q tensor coefficients");
  params.addRequiredParam<std::vector<std::string>>(
      "C_ij", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("C_ij_val", "Elastic stiffness tensor coefficients");
  params.addRequiredParam<std::vector<std::string>>(
      "permittivity", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<Real>>("permittivity_val", "static (high freq) permitivitty");
  params.addRequiredParam<bool>("coupled_problem","Asks if the problem is coupled (adds the elasticity variables).");
  params.addRequiredParam<bool>("polar_time_dependence","Asks if the polarization is time dependent.");
  params.addParam<bool>("u_time_dependence", false, "Asks if the elastic displacement is time dependent.");
  params.addParam<bool>("phi_time_dependence", false, "Asks if the electrostatic potential is time dependent.");
  params.addParam<bool>("is_renormalized", false, "Asks if the thermodynamic potential is renormalized.");
  params.addParam<bool>("is_permittivity_anisotropic", false, "Asks if the permittivity is anisotropic.");
  return params;
}

ABO3CoupledPhaseFieldAction::ABO3CoupledPhaseFieldAction(const InputParameters & params)
    : Action(params),
      _coupled_problem(getParam<bool>("coupled_problem")),
      _polar_time_dependence(getParam<bool>("polar_time_dependence")),
      _u_time_dependence(getParam<bool>("u_time_dependence")),
      _phi_time_dependence(getParam<bool>("phi_time_dependence")),
      _is_renormalized(getParam<bool>("is_renormalized")),
      _is_permittivity_anisotropic(getParam<bool>("is_permittivity_anisotropic"))
{
  // Do some error checking
  if(_coupled_problem==true)
  {
    std::vector<NonlinearVariableName> variables =
      getParam<std::vector<NonlinearVariableName>>("variables");
    if (variables.size() != 7)
      mooseError("There should be 7 variables for the coupled problem.");
  }
  else if (_coupled_problem==false)
  {
    std::vector<NonlinearVariableName> variables =
      getParam<std::vector<NonlinearVariableName>>("variables");
    if (variables.size() != 4)
      mooseError("There should be 4 variables for the uncoupled problem.");
  }
}

void
ABO3CoupledPhaseFieldAction::act()
{
  std::vector<NonlinearVariableName> variables =
    getParam<std::vector<NonlinearVariableName>>("variables");

  //if (_current_task == "add_variable")
  //{
  //}

  if (_current_task == "add_kernel")
  {
    unsigned int _ord_num = 3; //only polar_x polar_y and polar_z will be added in this loop
    for (unsigned int kk = 0; kk < _ord_num; ++kk)
    {
      if(_polar_time_dependence==true)
      {
        InputParameters params = _factory.getValidParams("TimeDerivativeScaled");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.applyParameters(parameters());

        std::string kernel_name = "tdsp_" + Moose::stringify(kk);
        _problem->addKernel("TimeDerivativeScaled", kernel_name, params);
      }
      {
        InputParameters params = _factory.getValidParams("BulkEnergyDerivativeEighth");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("polar_x") = {variables[0]};
        params.set<std::vector<VariableName>>("polar_y") = {variables[1]};
        params.set<std::vector<VariableName>>("polar_z") = {variables[2]};
        params.applyParameters(parameters());

        std::string kernel_name = "bed_" + Moose::stringify(kk);
        _problem->addKernel("BulkEnergyDerivativeEighth", kernel_name, params);
      }
      {
        InputParameters params = _factory.getValidParams("WallEnergyDerivative");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("polar_x") = {variables[0]};
        params.set<std::vector<VariableName>>("polar_y") = {variables[1]};
        params.set<std::vector<VariableName>>("polar_z") = {variables[2]};
        params.applyParameters(parameters());

        std::string kernel_name = "wed_" + Moose::stringify(kk);
        _problem->addKernel("WallEnergyDerivative", kernel_name, params);
      }
      {
        InputParameters params = _factory.getValidParams("PolarElectricPStrong");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("potential_E_int") = {variables[3]};
        params.applyParameters(parameters());

        std::string kernel_name = "pep_" + Moose::stringify(kk);
        _problem->addKernel("PolarElectricPStrong", kernel_name, params);
      }
      if(_coupled_problem==true)
      {
        InputParameters params = _factory.getValidParams("ElectrostrictiveCouplingPolarDerivative");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("polar_x") = {variables[0]};
        params.set<std::vector<VariableName>>("polar_y") = {variables[1]};
        params.set<std::vector<VariableName>>("polar_z") = {variables[2]};
        params.set<std::vector<VariableName>>("u_x") = {variables[4]};
        params.set<std::vector<VariableName>>("u_y") = {variables[5]};
        params.set<std::vector<VariableName>>("u_z") = {variables[6]};
        params.applyParameters(parameters());

        std::string kernel_name = "ecpd_" + Moose::stringify(kk);
        _problem->addKernel("ElectrostrictiveCouplingPolarDerivative", kernel_name, params);
      }
    }
    if(_phi_time_dependence==true)
    {
      InputParameters params = _factory.getValidParams("TimeDerivativeScaled");
      params.set<NonlinearVariableName>("variable") = variables[3];
      params.applyParameters(parameters());

      std::string kernel_name = "tdsphi";
      _problem->addKernel("TimeDerivativeScaled", kernel_name, params);
    }
    {
      InputParameters params = _factory.getValidParams("PolarElectricEStrong");
      params.set<NonlinearVariableName>("variable") = variables[3];
      params.set<std::vector<VariableName>>("polar_x") = {variables[0]};
      params.set<std::vector<VariableName>>("polar_y") = {variables[1]};
      params.set<std::vector<VariableName>>("polar_z") = {variables[2]};
      params.applyParameters(parameters());

      std::string kernel_name = "pees";
      _problem->addKernel("PolarElectricEStrong", kernel_name, params);
    }
    if (_is_permittivity_anisotropic==false)
    {
      InputParameters params = _factory.getValidParams("Electrostatics");
      params.set<NonlinearVariableName>("variable") = variables[3];
      params.applyParameters(parameters());

      std::string kernel_name = "e";
      _problem->addKernel("Electrostatics", kernel_name, params);
    }
    if(_coupled_problem==true)
    {
      unsigned int _ord_num = 3;
      for (unsigned int kk = 0; kk < _ord_num; ++kk)
      {
        if(_u_time_dependence==true)
        {
          InputParameters params = _factory.getValidParams("TimeDerivativeScaled");
          params.set<NonlinearVariableName>("variable") = variables[kk+4];
          params.applyParameters(parameters());

          std::string kernel_name = "tdsu_" + Moose::stringify(kk);
          _problem->addKernel("TimeDerivativeScaled", kernel_name, params);
        }
        InputParameters params = _factory.getValidParams("ElectrostrictiveCouplingDispDerivative");
        params.set<NonlinearVariableName>("variable") = variables[kk+4];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("polar_x") = {variables[0]};
        params.set<std::vector<VariableName>>("polar_y") = {variables[1]};
        params.set<std::vector<VariableName>>("polar_z") = {variables[2]};
        params.applyParameters(parameters());

        std::string kernel_name = "ecdd_" + Moose::stringify(kk);
        _problem->addKernel("ElectrostrictiveCouplingDispDerivative", kernel_name, params);
      }
    }
  }
  if (_current_task == "add_material")
  //note that these will change in the future because we want consistent Euler angle rotations across every material property
  {
    {
      InputParameters params = _factory.getValidParams("GenericConstantMaterial");

      params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("alpha_ijkl");
      params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("alpha_ijkl_val");

      _problem->addMaterial("GenericConstantMaterial", "Landau_FE_bulk_material", params);
    }
    {
      InputParameters params = _factory.getValidParams("GenericConstantMaterial");
      params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("G_ij");
      params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("G_ij_val");

      _problem->addMaterial("GenericConstantMaterial", "Landau_FE_gradient_material", params);
    }
    if (_is_permittivity_anisotropic==false)
    {
      InputParameters params = _factory.getValidParams("GenericConstantMaterial");
      params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("permittivity");
      params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("permittivity_val");

      _problem->addMaterial("GenericConstantMaterial", "FE_Poisson_material", params);
    }
    if(_coupled_problem==true)
    {
      {
        InputParameters params = _factory.getValidParams("GenericConstantMaterial");
        params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("Q_ij");
        params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("Q_ij_val");

        _problem->addMaterial("GenericConstantMaterial", "FE_Q_ij_material", params);
      }
      {
        InputParameters params = _factory.getValidParams("GenericConstantMaterial");
        params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("q_ij");
        params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("q_ij_val");

        _problem->addMaterial("GenericConstantMaterial", "FE_q_ij_material", params);
      }
      {
        InputParameters params = _factory.getValidParams("GenericConstantMaterial");
        params.set<std::vector<std::string>>("prop_names") = getParam<std::vector<std::string>>("C_ij");
        params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>("C_ij_val");

        _problem->addMaterial("GenericConstantMaterial", "FE_C_ij_material", params);
      }
    }
    //if (_current_task == "add_postprocessor")
    //{
    //}
  }
}
