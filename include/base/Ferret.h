#ifndef FERRET_H
#define FERRET_H
#include "Parser.h"
namespace Ferret
{
  /**
   * Registers all Kernels and BCs
   */
  void registerObjects(Factory & factory);

  /**
   * Associate actions to sections of input file.
   */
  void associateSyntax(Syntax&, ActionFactory & action_factory);
}

#endif //FERRET_H
