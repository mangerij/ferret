#ifndef FERRETBASE_H
#define FERRETBASE_H

#include <string>
#include <set>
#include "InputParameters.h"
#include "MooseObject.h"

class FerretBase 
{
 public:
  FerretBase(InputParameters parameters);
 protected:
  std::vector<std::string> _debug_vec;
  std::string              _class;
  std::set<std::string>    _debug_set;
  bool debug(const std::string&);
};

#endif //FERRETBASE_H
