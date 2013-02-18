#include "FerretApp.h"
#include "Ferret.h"
#include "Moose.h"
#include "Elk.h"

template<>
InputParameters validParams<FerretApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

FerretApp::FerretApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(libMesh::processor_id());
  
  Moose::registerObjects(_factory);
  Elk::registerObjects(_factory);
  Elk::associateSyntax(_syntax, _action_factory);
  Ferret::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  Ferret::associateSyntax(_syntax, _action_factory);
}

FerretApp::~FerretApp()
{
}

