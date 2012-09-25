#include "Moose.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "Parser.h"

#include "Elk.h"

#include "PolarizationVortex.h"
#include "PolarizationVortexAction.h"
//#include "VortexSurfaceEnergy.h"

namespace Ferret
{
  void registerObjects()
  {
    Elk::registerObjects();
    registerKernel(PolarizationVortex);
    //registerPostprocessor(VortexSurfaceEnergy);
  }
  
  void associateSyntax(Syntax& syntax)
  {
    Elk::associateSyntax(syntax);
    syntax.registerActionSyntax("PolarizationVortexAction","PolarizationVortex");
    registerAction(PolarizationVortexAction, "add_kernel");
  }
}
