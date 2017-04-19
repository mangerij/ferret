/**
 * @file   ComputeRotatedElectrostrictiveTensorBase.C
 * @author L. Kuna <lukasz.kuna@uconn.edu>
 *
 * @brief  Marker class to refine at specific points in the
 *         Landau-Ginzburg-Devonshire evolution
 *
 */


#include "PolarizationNWEMarker.h"
#include "FEProblem.h"
#include "MooseEnum.h"
#include <math.h>

template<>
InputParameters validParams<PolarizationNWEMarker>()
{
  InputParameters params = validParams<Marker>();

  MooseEnum third_state("DONT_MARK=-1 COARSEN DO_NOTHING REFINE", "DONT_MARK");
  params.addParam<MooseEnum>("third_state", third_state, "The Marker state to apply to values falling in-between the coarsen and refine thresholds.");
  params.addParam<Real>("coarsen", "The threshold value for coarsening.  Elements with variable values beyond this will be marked for coarsening.");
  params.addParam<Real>("refine", "The threshold value for refinement.  Elements with variable values beyond this will be marked for refinement.");
  params.addParam<bool>("invert", false, "If this is true then values _below_ 'refine' will be refined and _above_ 'coarsen' will be coarsened.");
  params.addParam<bool>("AMRoff", false, "If this is true then values _below_ 'refine' will be refined and _above_ 'coarsen' will be coarsened.");
  params.addParam<PostprocessorName>("ExtremeValue", "The magnitude of the maximum polarization");

  params.addRequiredCoupledVar("variable", "The values of this variable will be compared to 'refine' and 'coarsen' to see what should be done with the element");
  params.addRequiredParam<Real>("Bulk_Polar", "The bulk polarization of the material");
  params.addClassDescription("The the refinement state based on a threshold value compared to the specified variable.");
  return params;
}


PolarizationNWEMarker::PolarizationNWEMarker(const InputParameters & parameters) :
    QuadraturePointMarker(parameters),
    _coarsen_set(parameters.isParamValid("coarsen")),
    _coarsen(parameters.get<Real>("coarsen")),
    _refine_set(parameters.isParamValid("refine")),
    _refine(parameters.get<Real>("refine")),
    _PolarMag(getPostprocessorValue("ExtremeValue")),

    _invert(parameters.get<bool>("invert")),
    _AMRoff(parameters.get<bool>("AMRoff")),
    _third_state(getParam<MooseEnum>("third_state").getEnum<MarkerValue>()),

    _u(coupledValue("variable")),
    _Bulk_Polar(getParam<Real>("Bulk_Polar"))
{
}

Marker::MarkerValue //TODO: is this line needed? -John
PolarizationNWEMarker::computeQpMarker()
{
  Real ten_per = 0.20 * _Bulk_Polar;
  if (!_AMRoff)
  {
    if (_PolarMag > ten_per)
    {
      if (!_invert)
      {
        if (_refine_set && _u[_qp] > _refine)
          return REFINE;

        if (_coarsen_set && _u[_qp] < _coarsen)
          return COARSEN;
      }
      else
      {
        if (_refine_set && _u[_qp] < _refine)
          return REFINE;

        if (_coarsen_set && _u[_qp] > _coarsen)
          return COARSEN;
      }

      return _third_state;
    }
    return REFINE;
  }
  return REFINE;
}
