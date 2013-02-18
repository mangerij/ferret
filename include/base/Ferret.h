#ifndef FERRET_H
#define FERRET_H

class Factory;
class ActionFactory;
class Syntax;

namespace Ferret
{
  /**
   * Register this application and any it depends on.
   */
  void registerApps();

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
