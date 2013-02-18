#ifndef FERRETAPP_H
#define FERRETAPP_H

#include "MooseApp.h"

class FerretApp;

template<>
InputParameters validParams<FerretApp>();

class FerretApp : public MooseApp
{
public:
  FerretApp(const std::string & name, InputParameters parameters);
  virtual ~FerretApp();
};

#endif /* FERRETAPP_H */
