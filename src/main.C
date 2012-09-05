#include "FerretApp.h"
//Moose Includes
#include "Moose.h"
#include "MooseInit.h"

PerfLog Moose::perf_log("Ferret");

 // Begin the main program.
int main (int argc, char** argv)
{
  Moose::perf_log.push("main()","Ferret");

  MooseInit init (argc, argv);
  FerretApp app(argc, argv);
  app.setCheckUnusedFlag(true);
  app.run();

  Moose::perf_log.pop("main()","Ferret");

  return 0;
}
