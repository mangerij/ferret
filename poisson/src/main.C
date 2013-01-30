/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/**
 * Solve: -\Delta u=\nabla\cdot P
 * with DirichletBC and PeriodicBC
 * where P is a vector field 
 */

// Moose Includes
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "Factory.h"

#include "Magnetizing.h"
#include "PolarMaterial.h"
#include "ExampleMaterial.h"
// libMesh includes
#include "perf_log.h"

// Create a performance log
PerfLog Moose::perf_log("Poisson with DirichletBC and PeriodicBC");

// Begin the main program.
int main (int argc, char** argv)
{
  MooseInit init (argc, argv);
  MooseApp app(argc, argv);
  app.init();

  registerKernel(Magnetizing);
  registerMaterial(PolarMaterial);
  registerMaterial(ExampleMaterial);
  app.run();

  return 0;
}
