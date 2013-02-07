#include "FerretApp.h"
#include "Ferret.h"
#include "Moose.h"
#include "Elk.h"

FerretApp::FerretApp(int argc, char * argv[]) :
    MooseApp(argc, argv)
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

