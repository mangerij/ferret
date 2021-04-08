#include "ThermoelectricMaterial.h"
#include "Function.h"

#include "libmesh/quadrature.h"

registerMooseObject("FerretApp", ThermoelectricMaterial);
registerMooseObject("FerretApp", ADThermoelectricMaterial);

template <bool is_ad>
InputParameters
ThermoelectricMaterialTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();

  params.addCoupledVar("temp", "Coupled Temperature");

  params.addParam<Real>("thC", "The thermal conductivity value");
  params.addParam<FunctionName>(
      "thC_temperature_function", "", "Thermal conductivity as a function of temperature.");
  params.addParam<Real>("ecC", "The electrical conductivity value");
  params.addParam<FunctionName>(
      "ecC_temperature_function", "", "electrical conductivity as a function of temperature.");
  params.addParam<Real>("sbC", "The seebeck coefficient value");
  params.addParam<FunctionName>(
      "sbC_temperature_function", "", "seebeck coefficient as a function of temperature.");
  params.addClassDescription("General-purpose material model for thermoelectrics");

  return params;
}

template <bool is_ad>
ThermoelectricMaterialTempl<is_ad>::ThermoelectricMaterialTempl(const InputParameters & parameters)
  : Material(parameters),

    _has_temp(isCoupled("temp")),
    _T((_has_temp && !is_ad) ? coupledValue("temp") : _zero),
    _ad_T((_has_temp && is_ad) ? adCoupledValue("temp") : _ad_zero),
    _my_thC(isParamValid("thC") ? getParam<Real>("thC") : 0),
    _my_ecC(isParamValid("ecC") ? getParam<Real>("ecC") : 0),
    _my_sbC(isParamValid("sbC") ? getParam<Real>("sbC") : 0),

    _thC(declareGenericProperty<Real, is_ad>("thC")),
    _thC_dT(declareProperty<Real>("thC_dT")),
    _thC_temperature_function(getParam<FunctionName>("thC_temperature_function") != ""
                                  ? &getFunction("thC_temperature_function")
                                  : NULL),

    _ecC(declareGenericProperty<Real, is_ad>("ecC")),
    _ecC_dT(declareProperty<Real>("ecC_dT")),
    _ecC_temperature_function(getParam<FunctionName>("ecC_temperature_function") != ""
                                  ? &getFunction("ecC_temperature_function")
                                  : NULL),

    _sbC(declareGenericProperty<Real, is_ad>("sbC")),
    _sbC_dT(declareProperty<Real>("sbC_dT")),
    _sbC_temperature_function(getParam<FunctionName>("sbC_temperature_function") != ""
                                  ? &getFunction("sbC_temperature_function")
                                  : NULL)
{
  if (_thC_temperature_function && !_has_temp)
  {
    mooseError("Must couple with temperature if using thermal conductivity function");
  }
  if (isParamValid("thC") && _thC_temperature_function)
  {
    mooseError(
        "Cannot define both thermal conductivity and thermal conductivity temperature function");
  }
  if (_ecC_temperature_function && !_has_temp)
  {
    mooseError("Must couple with temperature if using electrical conductivity function");
  }
  if (isParamValid("ecC") && _ecC_temperature_function)
  {
    mooseError("Cannot define both electrical conductivity and electrical conductivity temperature "
               "function");
  }
  if (_sbC_temperature_function && !_has_temp)
  {
    mooseError("Must couple with temperature if using seebeck coefficient function");
  }
  if (isParamValid("sbC") && _sbC_temperature_function)
  {
    mooseError(
        "Cannot define both seebeck coefficient and seebeck coefficient temperature function");
  }
}

template <bool is_ad>
void
ThermoelectricMaterialTempl<is_ad>::setDerivatives(GenericReal<is_ad> & prop,
                                                   Real dprop_dT,
                                                   const ADReal & ad_T)
{
  if (ad_T < 0)
    prop.derivatives() = 0;
  else
    prop.derivatives() = dprop_dT * ad_T.derivatives();
}

template <>
void
ThermoelectricMaterialTempl<false>::setDerivatives(Real &, Real, const ADReal &)
{
  mooseError("Mistaken call of setDerivatives in a non-AD ThermoelectricMaterial version");
}

template <bool is_ad>
void
ThermoelectricMaterialTempl<is_ad>::computeProperties()
{
  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    Real qp_T = 0;
    if (_has_temp)
    {
      if (is_ad)
        qp_T = MetaPhysicL::raw_value(_ad_T[qp]);
      else
        qp_T = _T[qp];
      if (qp_T < 0)
      {
        std::stringstream msg;
        msg << "WARNING:  In ThermoelectricMaterial:  negative temperature!\n"
            << "\tResetting to zero.\n"
            << "\t_qp: " << qp << "\n"
            << "\ttemp: " << qp_T << "\n"
            << "\telem: " << _current_elem->id() << "\n"
            << "\tproc: " << processor_id() << "\n";
        mooseWarning(msg.str());
        qp_T = 0;
      }
    }
    if (_thC_temperature_function)
    {
      Point p;
      _thC[qp] = _thC_temperature_function->value(qp_T, p);
      // A terrible exploitation of the Function API to get a derivative with respect to a
      // non-linear variable
      _thC_dT[qp] = _thC_temperature_function->timeDerivative(qp_T, p);
      if (is_ad)
        setDerivatives(_thC[qp], _thC_dT[qp], _ad_T[qp]);
    }
    else
    {
      _thC[qp] = _my_thC;
      _thC_dT[qp] = 0;
    }

    if (_ecC_temperature_function)
    {
      Point p;
      _ecC[qp] = _ecC_temperature_function->value(qp_T, p);
      // A terrible exploitation of the Function API to get a derivative with respect to a
      // non-linear variable
      _ecC_dT[qp] = _ecC_temperature_function->timeDerivative(qp_T, p);
      if (is_ad)
        setDerivatives(_ecC[qp], _ecC_dT[qp], _ad_T[qp]);
    }
    else
    {
      _ecC[qp] = _my_ecC;
      _ecC_dT[qp] = 0;
    }

    if (_sbC_temperature_function)
    {
      Point p;
      _sbC[qp] = _sbC_temperature_function->value(qp_T, p);
      // A terrible exploitation of the Function API to get a derivative with respect to a
      // non-linear variable
      _sbC_dT[qp] = _sbC_temperature_function->timeDerivative(qp_T, p);
      if (is_ad)
        setDerivatives(_sbC[qp], _sbC_dT[qp], _ad_T[qp]);
    }
    else
    {
      _sbC[qp] = _my_sbC;
      _sbC_dT[qp] = 0;
    }
  }
}

template class ThermoelectricMaterialTempl<false>;
template class ThermoelectricMaterialTempl<true>;
