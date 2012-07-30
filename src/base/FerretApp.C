#include "FerretApp.h"
#include "Ferret.h"
#include "Moose.h"
#include "Elk.h"

FerretApp::FerretApp(int argc, char * argv[]) :
    MooseApp(argc, argv)
{
  srand(libMesh::processor_id());

  init();

  Elk::registerObjects();
  Elk::associateSyntax(_syntax);

  Ferret::registerObjects();
}

FerretApp::~FerretApp()
{
}

