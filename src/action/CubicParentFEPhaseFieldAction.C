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

#include "CubicParentFEPhaseFieldAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_variable");
registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_aux_variable");
registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_aux_kernel");
registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_kernel");
registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_material");
registerMooseAction("FerretApp", CubicParentFEPhaseFieldAction, "add_postprocessor");

InputParameters
CubicParentFEPhaseFieldAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Set up a cubic-parent-phase ferroelectric problem: the polarization variables, the Landau "
      "bulk and gradient kernels, the optional electrostatic and electrostrictive couplings, the "
      "constant material properties and the free energy postprocessors. The mechanics themselves "
      "are left to the SolidMechanics QuasiStatic physics, which this action is designed to be "
      "used alongside; couple the two by listing 'eigenstrain_name' in its 'eigenstrain_names'.");

  params.addParam<std::vector<NonlinearVariableName>>(
      "polar_vars",
      std::vector<NonlinearVariableName>{"polar_x", "polar_y", "polar_z"},
      "The three polarization vector components.");
  params.addParam<NonlinearVariableName>(
      "potential_var", "potential_E_int", "The internal electrostatic potential variable.");
  params.addCoupledVar("displacements",
                       "The displacement variables, as given to the SolidMechanics action. When "
                       "supplied the P<-u Jacobian block is assembled.");

  params.addParam<bool>("electrostatics", true, "Solve Poisson's equation for the depolarization "
                                                "field and couple it to the polarization.");
  params.addParam<bool>("elastic", true, "Couple the polarization to the strain through the "
                                         "cubic-parent electrostrictive energy.");
  params.addParam<bool>("polar_time_dependence", true, "Add the TimeDerivative kernels that make "
                                                       "this a time-dependent (TDLGD) problem.");
  params.addParam<bool>("add_variables", true, "Add the polarization variables, and the potential "
                                               "when 'electrostatics' is true.");
  params.addParam<bool>("add_elastic_materials", true,
                        "Also add ComputeElasticityTensor and ComputeLinearElasticStress from "
                        "'C_ij_val'. Turn off to supply your own stress model.");
  params.addParam<bool>(
      "add_elastic_energy_aux", true,
      "Add the ElasticEnergyAux route to the elastic energy: an f_el aux variable holding "
      "1/2 sigma:elastic_strain, integrated into the 'Felastic_true' postprocessor. Because "
      "the strain calculator has already removed EVERY eigenstrain in 'eigenstrain_names', "
      "this is the true elastic energy whatever else contributes one, whereas "
      "CubicParentElasticEnergy subtracts only its own. When on, 'Ftotal' uses it.");
  params.addParam<bool>("add_postprocessors", true, "Add the free energy postprocessors.");

  params.addRequiredParam<std::vector<std::string>>("alpha_ijkl",
                                                    "Names of the Landau expansion coefficients.");
  params.addRequiredParam<std::vector<Real>>("alpha_ijkl_val",
                                             "Values of the Landau expansion coefficients.");
  params.addRequiredParam<std::vector<std::string>>("G_ij", "Names of the gradient coefficients.");
  params.addRequiredParam<std::vector<Real>>("G_ij_val", "Values of the gradient coefficients.");
  params.addParam<std::vector<std::string>>(
      "Q_ij", std::vector<std::string>{"Q11", "Q12", "Q44"},
      "Names of the electrostrictive Q tensor coefficients.");
  params.addParam<std::vector<Real>>("Q_ij_val",
                                     "Values of the electrostrictive Q tensor coefficients.");
  params.addParam<std::vector<std::string>>(
      "C_ij", std::vector<std::string>{"C11", "C12", "C44"},
      "Names of the elastic stiffness coefficients. Must be 'C11 C12 C44' (cubic parent).");
  params.addParam<std::vector<Real>>("C_ij_val", "Values of the elastic stiffness coefficients.");
  params.addParam<std::vector<std::string>>(
      "permittivity", std::vector<std::string>{"permittivity"},
      "Name of the permittivity material property.");
  params.addParam<std::vector<Real>>("permittivity_val", "Value of the permittivity.");

  params.addParam<std::string>("eigenstrain_name", "ferro",
                               "Name of the spontaneous polar eigenstrain. List this same name in "
                               "the SolidMechanics action's 'eigenstrain_names'.");
  params.addParam<std::string>("base_name", "Material property base name.");
  params.addParam<Real>("len_scale", 1.0, "The length scale of the unit.");
  params.addParam<Real>("energy_scale", 1.0, "The energy scale of the postprocessors.");
  params.addParam<std::vector<SubdomainName>>(
      "block", "The blocks this action operates on. Defaults to the whole mesh.");

  return params;
}

CubicParentFEPhaseFieldAction::CubicParentFEPhaseFieldAction(const InputParameters & params)
  : Action(params),
    _polar_vars(getParam<std::vector<NonlinearVariableName>>("polar_vars")),
    _potential_var(getParam<NonlinearVariableName>("potential_var")),
    _electrostatics(getParam<bool>("electrostatics")),
    _elastic(getParam<bool>("elastic")),
    _polar_time_dependence(getParam<bool>("polar_time_dependence")),
    _add_variables(getParam<bool>("add_variables")),
    _add_elastic_materials(getParam<bool>("add_elastic_materials")),
    _add_elastic_energy_aux(getParam<bool>("add_elastic_energy_aux")),
    _add_postprocessors(getParam<bool>("add_postprocessors")),
    _eigenstrain_name(getParam<std::string>("eigenstrain_name"))
{
  if (_polar_vars.size() != 3)
    paramError("polar_vars",
               "The cubic parent functional is written in all three polarization components, so "
               "exactly 3 names are needed (got ",
               _polar_vars.size(),
               "). In 2D still pass three and constrain the third through its initial condition.");

  if (_elastic)
  {
    if (!isParamValid("Q_ij_val"))
      paramError("Q_ij_val", "Required when 'elastic' is true.");
    if (!isParamValid("C_ij_val"))
      paramError("C_ij_val", "Required when 'elastic' is true.");
    if (getParam<std::vector<std::string>>("Q_ij").size() !=
        getParam<std::vector<Real>>("Q_ij_val").size())
      paramError("Q_ij_val", "Must have the same number of entries as 'Q_ij'.");
    if (getParam<std::vector<std::string>>("C_ij").size() !=
        getParam<std::vector<Real>>("C_ij_val").size())
      paramError("C_ij_val", "Must have the same number of entries as 'C_ij'.");
    if (_add_elastic_materials && getParam<std::vector<Real>>("C_ij_val").size() != 3)
      paramError("C_ij_val",
                 "With 'add_elastic_materials' the elasticity tensor is built from the cubic "
                 "triple 'C11 C12 C44', so exactly 3 values are needed.");
  }

  if (_electrostatics && !isParamValid("permittivity_val"))
    paramError("permittivity_val", "Required when 'electrostatics' is true.");

  if (getParam<std::vector<std::string>>("alpha_ijkl").size() !=
      getParam<std::vector<Real>>("alpha_ijkl_val").size())
    paramError("alpha_ijkl_val", "Must have the same number of entries as 'alpha_ijkl'.");
  if (getParam<std::vector<std::string>>("G_ij").size() !=
      getParam<std::vector<Real>>("G_ij_val").size())
    paramError("G_ij_val", "Must have the same number of entries as 'G_ij'.");
}

void
CubicParentFEPhaseFieldAction::setPolarCoupling(InputParameters & params) const
{
  params.set<std::vector<VariableName>>("polar_x") = {_polar_vars[0]};
  params.set<std::vector<VariableName>>("polar_y") = {_polar_vars[1]};
  params.set<std::vector<VariableName>>("polar_z") = {_polar_vars[2]};
}

void
CubicParentFEPhaseFieldAction::addConstantMaterial(const std::string & prop_param,
                                                   const std::string & val_param,
                                                   const std::string & object_name)
{
  InputParameters params = _factory.getValidParams("GenericConstantMaterial");
  params.set<std::vector<std::string>>("prop_names") =
      getParam<std::vector<std::string>>(prop_param);
  params.set<std::vector<Real>>("prop_values") = getParam<std::vector<Real>>(val_param);
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") =
        getParam<std::vector<SubdomainName>>("block");
  _problem->addMaterial("GenericConstantMaterial", object_name, params);
}

void
CubicParentFEPhaseFieldAction::act()
{
  const bool have_disp = isParamValid("displacements");

  if (_current_task == "add_variable" && _add_variables)
  {
    InputParameters params = _factory.getValidParams("MooseVariable");
    params.set<MooseEnum>("order") = "FIRST";
    params.set<MooseEnum>("family") = "LAGRANGE";
    if (isParamValid("block"))
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");

    for (const auto & pvar : _polar_vars)
      _problem->addVariable("MooseVariable", pvar, params);

    if (_electrostatics)
      _problem->addVariable("MooseVariable", _potential_var, params);
  }

  if (_current_task == "add_aux_variable" && _elastic && _add_elastic_energy_aux)
  {
    InputParameters params = _factory.getValidParams("MooseVariableConstMonomial");
    params.set<MooseEnum>("order") = "CONSTANT";
    params.set<MooseEnum>("family") = "MONOMIAL";
    if (isParamValid("block"))
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    _problem->addAuxVariable("MooseVariableConstMonomial", "f_el", params);
  }

  if (_current_task == "add_aux_kernel" && _elastic && _add_elastic_energy_aux)
  {
    InputParameters params = _factory.getValidParams("ElasticEnergyAux");
    params.set<AuxVariableName>("variable") = "f_el";
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
    if (isParamValid("base_name"))
      params.set<std::string>("base_name") = getParam<std::string>("base_name");
    if (isParamValid("block"))
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    _problem->addAuxKernel("ElasticEnergyAux", "cpfe_f_el", params);
  }

  if (_current_task == "add_kernel")
  {
    for (unsigned int kk = 0; kk < 3; ++kk)
    {
      if (_polar_time_dependence)
      {
        InputParameters params = _factory.getValidParams("TimeDerivative");
        params.set<NonlinearVariableName>("variable") = _polar_vars[kk];
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel("TimeDerivative", "cpfe_time_" + Moose::stringify(kk), params);
      }
      {
        InputParameters params = _factory.getValidParams("BulkEnergyDerivativeEighth");
        params.set<NonlinearVariableName>("variable") = _polar_vars[kk];
        params.set<unsigned int>("component") = kk;
        setPolarCoupling(params);
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel(
            "BulkEnergyDerivativeEighth", "cpfe_bed_" + Moose::stringify(kk), params);
      }
      {
        InputParameters params = _factory.getValidParams("WallEnergyDerivative");
        params.set<NonlinearVariableName>("variable") = _polar_vars[kk];
        params.set<unsigned int>("component") = kk;
        setPolarCoupling(params);
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel("WallEnergyDerivative", "cpfe_wed_" + Moose::stringify(kk), params);
      }
      if (_elastic)
      {
        InputParameters params = _factory.getValidParams("CubicParentElasticPDerivative");
        params.set<NonlinearVariableName>("variable") = _polar_vars[kk];
        params.set<unsigned int>("component") = kk;
        setPolarCoupling(params);
        if (have_disp)
          params.set<std::vector<VariableName>>("displacements") =
              getParam<std::vector<VariableName>>("displacements");
        if (isParamValid("base_name"))
          params.set<std::string>("base_name") = getParam<std::string>("base_name");
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel(
            "CubicParentElasticPDerivative", "cpfe_estr_" + Moose::stringify(kk), params);
      }
      if (_electrostatics)
      {
        InputParameters params = _factory.getValidParams("PolarElectricPStrong");
        params.set<NonlinearVariableName>("variable") = _polar_vars[kk];
        params.set<unsigned int>("component") = kk;
        params.set<std::vector<VariableName>>("potential_E_int") = {_potential_var};
        params.set<Real>("len_scale") = getParam<Real>("len_scale");
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel("PolarElectricPStrong", "cpfe_pep_" + Moose::stringify(kk), params);
      }
    }

    if (_electrostatics)
    {
      {
        InputParameters params = _factory.getValidParams("Electrostatics");
        params.set<NonlinearVariableName>("variable") = _potential_var;
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel("Electrostatics", "cpfe_poisson", params);
      }
      {
        InputParameters params = _factory.getValidParams("PolarElectricEStrong");
        params.set<NonlinearVariableName>("variable") = _potential_var;
        setPolarCoupling(params);
        params.set<Real>("len_scale") = getParam<Real>("len_scale");
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addKernel("PolarElectricEStrong", "cpfe_divP", params);
      }
    }
  }

  if (_current_task == "add_material")
  {
    addConstantMaterial("alpha_ijkl", "alpha_ijkl_val", "cpfe_landau_bulk");
    addConstantMaterial("G_ij", "G_ij_val", "cpfe_landau_gradient");

    if (_electrostatics)
      addConstantMaterial("permittivity", "permittivity_val", "cpfe_permittivity");

    if (_elastic)
    {
      addConstantMaterial("Q_ij", "Q_ij_val", "cpfe_electrostrictive_Q");
      addConstantMaterial("C_ij", "C_ij_val", "cpfe_elastic_C");

      {
        InputParameters params =
            _factory.getValidParams("ComputeCubicParentElectrostrictiveStrain");
        params.set<std::string>("eigenstrain_name") = _eigenstrain_name;
        setPolarCoupling(params);
        if (isParamValid("base_name"))
          params.set<std::string>("base_name") = getParam<std::string>("base_name");
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addMaterial(
            "ComputeCubicParentElectrostrictiveStrain", "cpfe_eigenstrain", params);
      }

      if (_add_elastic_materials)
      {
        const std::vector<Real> & c = getParam<std::vector<Real>>("C_ij_val");
        const Real c11 = c[0], c12 = c[1], c44 = c[2];
        {
          InputParameters params = _factory.getValidParams("ComputeElasticityTensor");
          params.set<MooseEnum>("fill_method") = "symmetric9";
          params.set<std::vector<Real>>("C_ijkl") = {
              c11, c12, c12, c11, c12, c11, c44, c44, c44};
          if (isParamValid("base_name"))
            params.set<std::string>("base_name") = getParam<std::string>("base_name");
          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          _problem->addMaterial("ComputeElasticityTensor", "cpfe_elasticity_tensor", params);
        }
        {
          InputParameters params = _factory.getValidParams("ComputeLinearElasticStress");
          if (isParamValid("base_name"))
            params.set<std::string>("base_name") = getParam<std::string>("base_name");
          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          _problem->addMaterial("ComputeLinearElasticStress", "cpfe_stress", params);
        }
      }
    }
  }

  if (_current_task == "add_postprocessor" && _add_postprocessors)
  {
    std::vector<PostprocessorName> pp_names;
    std::vector<Real> pp_coefs;

    {
      InputParameters params = _factory.getValidParams("BulkEnergyEighth");
      setPolarCoupling(params);
      params.set<Real>("energy_scale") = getParam<Real>("energy_scale");
      params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
      if (isParamValid("block"))
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      _problem->addPostprocessor("BulkEnergyEighth", "Fbulk", params);
      pp_names.push_back("Fbulk");
      pp_coefs.push_back(1.0);
    }
    {
      InputParameters params = _factory.getValidParams("WallEnergy");
      setPolarCoupling(params);
      params.set<Real>("energy_scale") = getParam<Real>("energy_scale");
      params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
      if (isParamValid("block"))
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      _problem->addPostprocessor("WallEnergy", "Fwall", params);
      pp_names.push_back("Fwall");
      pp_coefs.push_back(1.0);
    }
    if (_elastic)
    {
      {
        InputParameters params = _factory.getValidParams("CubicParentElasticEnergy");
        setPolarCoupling(params);
        params.set<Real>("energy_scale") = getParam<Real>("energy_scale");
        params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
        if (isParamValid("base_name"))
          params.set<std::string>("base_name") = getParam<std::string>("base_name");
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addPostprocessor("CubicParentElasticEnergy", "Felastic", params);
      }
      if (_add_elastic_energy_aux)
      {
        InputParameters params =
            _factory.getValidParams("ElementIntegralVariablePostprocessor");
        params.set<std::vector<VariableName>>("variable") = {"f_el"};
        params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
        if (isParamValid("block"))
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        _problem->addPostprocessor(
            "ElementIntegralVariablePostprocessor", "Felastic_true", params);
      }
      pp_names.push_back(_add_elastic_energy_aux ? "Felastic_true" : "Felastic");
      pp_coefs.push_back(1.0);
    }
    if (_electrostatics)
    {
      InputParameters params = _factory.getValidParams("ElectrostaticEnergy");
      setPolarCoupling(params);
      params.set<std::vector<VariableName>>("potential_E_int") = {_potential_var};
      params.set<Real>("len_scale") = getParam<Real>("len_scale");
      params.set<Real>("energy_scale") = getParam<Real>("energy_scale");
      params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
      if (isParamValid("block"))
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      _problem->addPostprocessor("ElectrostaticEnergy", "Felec", params);
      pp_names.push_back("Felec");
      pp_coefs.push_back(1.0);
    }
    {
      InputParameters params = _factory.getValidParams("LinearCombinationPostprocessor");
      params.set<std::vector<PostprocessorName>>("pp_names") = pp_names;
      params.set<std::vector<Real>>("pp_coefs") = pp_coefs;
      params.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_END};
      _problem->addPostprocessor("LinearCombinationPostprocessor", "Ftotal", params);
    }
  }
}
