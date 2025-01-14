#ifndef BOUNDARYINTEGRALFMM_H
#define BOUNDARYINTEGRALFMM_H

#include "GeneralUserObject.h"
#include "MooseMesh.h"

// Forward declarations
class SystemBase;
class BoundaryIntegralFMM : public GeneralUserObject
{
public:
  BoundaryIntegralFMM(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void initialize() override;

  virtual void execute() override;

  virtual void finalize() override;

protected:


private:
  /// FMM parameters
  Real _cx, _cy, _cz, _boxWidth;
  unsigned int _TreeHeight;
};

#endif
