#include "PolarizationVortexAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

template<>
InputParameters validParams<PolarizationVortexAction>()
{
  InputParameters params = validParams<Action>();
  params.addParam<Real>("a_x", 0, "Vortex location along x");
  params.addParam<Real>("a_y", 0, "Vortex location along y");
  params.addParam<Real>("c", 0, "Vortex size");
  params.addParam<std::string>("time", "a_x", "Vortex evolution parameter: a_x, a_y or c.");

  return params;
}

PolarizationVortexAction::PolarizationVortexAction(const std::string & name, InputParameters params) :
  Action(name, params),
  _P_x(getParam<NonlinearVariableName>("P_x")),
  _P_y(getParam<NonlinearVariableName>("P_y")),
  _P_z(getParam<NonlinearVariableName>("P_z"))
{
}

void
PolarizationVortexAction::act()
{
  Real a_x = getParam<Real>("a_x");
  Real a_y = getParam<Real>("a_y");
  Real c = getParam<Real>("c");
  std::string p = getParam<std::string>("time");


  /**
   * We need to manually set up the PolarizationVortex kernel parameters.
   * Much of the syntax below is usualy hidden by the parser system, but 
   * we have to set things up ourselves this time.
   */

  // Do some error checking
  {
    std::stringstream err;
    err << "Expected nonnegative c, received " << c;
    mooseAssert(c >= 0, err.str());
  }
  {
    std::stringstream err;
    err << "Expected evolution parameter to be a_x, a_y or c, received " <<  p;
    mooseAssert(p != "a_x" && p != "a_y" && p != "c", err.str());
  }

  // Set up PolarizationVortex kernel on the P_x variable
  InputParameters polarization_params = Factory::instance()->getValidParams("Polarization");
  polarization_params.set<Real>("a_x") = a_x;
  polarization_params.set<Real>("a_y") = a_y;
  polarization_params.set<Real>("c") = c;
  polarization_params.set<std::string>("time") = p;
  NonlinearVariableName var, P_1, P_2;
  for(unsigned int i = 0; i < 3; ++i) {
    switch(i) {
    case 0:
      var = _P_x;
      P_1 = _P_y;
      P_2 = _P_z;
      break;
    case 1:
      var = _P_y;
      P_1 = _P_z;
      P_2 = _P_x;
      break;
    case 2:
      var = _P_z;
      P_1 = _P_x;
      P_2 = _P_y;
      break;
    }
    polarization_params.set<NonlinearVariableName>("variable") = var;
    polarization_params.set<NonlinearVariableName>("P_1") = P_1;
    polarization_params.set<NonlinearVariableName>("P_2") = P_2;
    
    std::stringstream name;
    name << "polarization_" << i;    
    _problem->addKernel("Polarization", name.str(), polarization_params);
  }
}

