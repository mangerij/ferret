#ifndef FERRETAPP_H
#define FERRETAPP_H

#include "MooseApp.h"

class FerretApp;

template<>
InputParameters validParams<FerretApp>();

class FerretApp : public MooseApp
{
public:
  FerretApp(const InputParameters & parameters);
  virtual ~FerretApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax& syntax, ActionFactory & action_factory);
};

#endif /* FERRETAPP_H */
