#ifndef POLARIZATIONNWEMARKER_H
#define POLARIZATIONNWEMARKER_H

#include "QuadraturePointMarker.h"

class PolarizationNWEMarker;

template<>
InputParameters validParams<PolarizationNWEMarker>();

class PolarizationNWEMarker : public QuadraturePointMarker
{
public:
  PolarizationNWEMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeQpMarker() override;

  bool _coarsen_set;
  Real _coarsen;
  bool _refine_set;
  Real _refine;
  const PostprocessorValue & _PolarMag;
  bool _invert;
  bool _AMRoff;
  MarkerValue _third_state;

  const VariableValue & _u;
  const Real _Bulk_Polar;

};

#endif /* VALUETHRESHOLDMARKER_H */
