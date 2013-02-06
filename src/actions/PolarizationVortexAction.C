#include "PolarizationVortexAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

template<>
InputParameters validParams<PolarizationVortexAction>()
{
  InputParameters params = validParams<Action>();
  params.addParam<std::vector<SubdomainName> >("block", "The list of ids of the blocks (subdomain) that these kernels will be applied to");
  params.addParam<Real>("R", 1, "Magnetic dot radius");
  params.addParam<Real>("L", 1, "Magnetic dot thickness");
  params.addParam<NonlinearVariableName>("P_x", "", "The x component of polarization");
  params.addParam<NonlinearVariableName>("P_y", "", "The y component of polarization");
  params.addParam<NonlinearVariableName>("P_z", "", "The z component of polarization");
  params.addParam<Real>("a_x", 0, "Vortex location along x");
  params.addParam<Real>("a_y", 0, "Vortex location along y");
  params.addParam<Real>("c", 1, "Vortex size");
  params.addParam<std::string>("p", "a_x", "Vortex evolution parameter: a_x, a_y or c.");
  params.addParam<bool>("debug", false, "Debugging flag");
  params.addParam<bool>("kernel_debug", false, "Debugging flag passed to PolarizationVortex kernels");
  return params;
}

PolarizationVortexAction::PolarizationVortexAction(const std::string & name, InputParameters params) :
  Action(name, params),
  _P_x(getParam<NonlinearVariableName>("P_x")),
  _P_y(getParam<NonlinearVariableName>("P_y")),
  _P_z(getParam<NonlinearVariableName>("P_z")),
  _a_x(getParam<Real>("a_x")),
  _a_y(getParam<Real>("a_y")),
  _c(getParam<Real>("c")),
  _R(getParam<Real>("R")),
  _L(getParam<Real>("L")),
  _p(getParam<std::string>("p")),
  _debug(getParam<bool>("debug"))
{
  if(_debug){
    libMesh::out << "PolarizationVortexAction: debug = " << _debug << "\n";
    libMesh::out << "_P_x = " << _P_x << ", _P_y = " << _P_y << ", _P_z = " << _P_z << ", _a_x = " << _a_x << ", _a_y = " << _a_y << ", _c = " << _c << ", R = " << _R << ", L = " << _L << ", _p = " << _p << "\n";
  }
  // Do some error checking
  {
    std::stringstream err;
    err << "Expected nonnegative c, received " << _c;
    mooseAssert(_c >= 0, err.str());
  }
  {
    std::stringstream err;
    err << "Expected positive R, received " << _R;
    mooseAssert(_R >  0, err.str());
  }
  {
    std::stringstream err;
    err << "Expected positive L, received " << _L;
    mooseAssert(_L >  0, err.str());
  }
  {
    std::stringstream err;
    err << "Expected evolution parameter to be a_x, a_y or c, received " <<  _p;
    mooseAssert(_p == "a_x" || _p == "a_y" || _p == "c", err.str());
  }
}

void
PolarizationVortexAction::act()
{

  /*
   * We need to manually set up the PolarizationVortex kernel parameters.
   * Much of this is usualy hidden by the parser system, but 
   * we have to set things up ourselves this time.
   */
 
  std::set<SubdomainID> subdomains;
  std::vector<SubdomainName> subdomainnames;
  // Set up PolarizationVortex kernel, one per each component of polarization and one per mesh block.
  if(isParamValid("block")) // Should it be restricted to certain blocks?
    {
      std::cout<<"Restricting to blocks!"<<std::endl;
      std::vector<SubdomainName> block = getParam<std::vector<SubdomainName> >("block");
      for(unsigned int i=0; i < block.size(); i++)
	subdomains.insert(_problem->mesh().getSubdomainID(block[i]));
    }  
  else // Put it everywhere
    subdomains = _problem->mesh().meshSubdomains();

 for (std::set<SubdomainID>::const_iterator it = subdomains.begin(); it != subdomains.end(); ++it)
  {
    SubdomainID sid = *it;
      // Convert the SubdomainID into SubdomainName since kernel params take SubdomainNames (this is retarded...)
    std::stringstream ss;
    ss << sid;
    SubdomainName sname = ss.str();
    subdomainnames.push_back(sname);
  }

  InputParameters polarization_vortex_params = _factory.getValidParams("PolarizationVortex");
  polarization_vortex_params.set<Real>("a_x") = _a_x;
  polarization_vortex_params.set<Real>("a_y") = _a_y;
  polarization_vortex_params.set<Real>("c") = _c;
  polarization_vortex_params.set<std::string>("p") = _p;
  NonlinearVariableName P;
  std::vector<NonlinearVariableName> Pnormal;
  for(unsigned int i = 0; i < 3; ++i) {
    switch(i) {
    case 0:
      P = _P_x;
      Pnormal.push_back(_P_y);
      Pnormal.push_back(_P_z);
      break;
    case 1:
     P = _P_y;
      Pnormal.push_back(_P_z);
      Pnormal.push_back(_P_x);
      break;
    case 2:
     P = _P_z;
      Pnormal.push_back(_P_x);
      Pnormal.push_back(_P_y);
      break;
    }
    polarization_vortex_params.set<Real>("R") = _R;
    polarization_vortex_params.set<Real>("L") = _L;
    polarization_vortex_params.set<unsigned int>("i") = i;
    polarization_vortex_params.set<NonlinearVariableName>("variable") = P;
    polarization_vortex_params.set<std::vector<NonlinearVariableName> >("Pnormal") = Pnormal;
    polarization_vortex_params.set<std::vector<SubdomainName> >("block") = subdomainnames;
    polarization_vortex_params.set<bool>("debug") = getParam<bool>("kernel_debug");
    std::stringstream name;
    name << "polarization_" << i;    
    _problem->addKernel("PolarizationVortex", name.str(), polarization_vortex_params);
  }
}

