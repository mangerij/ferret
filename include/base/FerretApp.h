#ifndef FERRETAPP_H
#define FERRETAPP_H

#include "MooseApp.h"

class FerretApp : public MooseApp
{
public:
  FerretApp(int argc, char * argv[]);
  virtual ~FerretApp();
};

#endif /* FERRETAPP_H */
