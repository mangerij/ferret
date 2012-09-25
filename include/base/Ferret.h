#ifndef FERRET_H
#define FERRET_H
#include "Parser.h"
namespace Ferret
{
  /**
   * Registers all Kernels and BCs
   */
  void registerObjects();

  /**
   * Associate actions to sections of input file.
   */
  void associateSyntax(Syntax&);
}

#endif //FERRET_H
