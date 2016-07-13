#include "FerretBase.h"

template<>
InputParameters validParams<FerretBase>()
{
  InputParameters params = emptyInputParameters();
  params.addParam<std::vector<std::string> >("Debug", std::vector<std::string>(), "Debug strings");
  params.addParam<std::string>("class","","Class name");
  return params;
}


FerretBase::FerretBase(InputParameters parameters):
  _debug_vec(parameters.get<std::vector<std::string> >("Debug")),
  _class(parameters.get<std::string>("class"))
{
  /// std::copy(_debug_vec.begin(),_debug_vec.end(),std::insert_iterator<std::set<std::string> >(_debug_set,_debug_set.begin()));
  for (unsigned int i = 0; i < _debug_vec.size(); ++i) {
    _debug_set.insert(_debug_vec[i]);
  }
}

bool 
FerretBase::debug(const std::string& str)
{
  std::string s;
  if (_class.size()) s = _class+"::"+str; else s = str;
  unsigned int count = _debug_set.count(s);
  bool res = (count > 0);
  return res;
}
