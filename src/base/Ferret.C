#include "Moose.h"
#include "Factory.h"
#include "VortexSurfaceCharge.h"
#include "VortexSurfaceEnergy.h"

namespace Ferret
{
  void registerObjects()
  {
    registerBoundaryCondition(VortexSurfaceCharge);
    registerPostprocessor(VortexSurfaceEnergy);
  }
}
