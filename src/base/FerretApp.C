#include "FerretApp.h"
#include "Ferret.h"
#include "Moose.h"

FerretApp::FerretApp(int argc, char * argv[]) :
    MooseApp(argc, argv)
{
  srand(libMesh::processor_id());

  init();

  Ferret::registerObjects();
  Ferret::associateSyntax(_syntax);
}

FerretApp::~FerretApp()
{
}

