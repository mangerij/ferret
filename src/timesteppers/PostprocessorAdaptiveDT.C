/**
 * @file   PostprocessorAdaptiveDT.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Fri Nov 15 13:17:57 2013
 *
 * @brief  Use pps for adaptive time stepping
 *
 *
 */
#include "PostprocessorAdaptiveDT.h"
#include "FEProblem.h"
//libMesh includes
#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"
// C++ Includes
#include <iomanip>
#include<limits>

template<>
InputParameters validParams<PostprocessorAdaptiveDT>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor that computes the dt");
  params.addParam<Real>("dt", 1., "The initial time step size.");
  params.addParam<Real>("increase_rate",2.0, "the rate the step size increases");
  params.addParam<Real>("decrease_rate",0.5,"the rate the step size decrease");
  params.addParam<unsigned int>("update_step",4,"try to increase step size every k successful steps, 0 means never");
  params.addParam<size_t>("max_record_size",100,"maximum pps record size");
  return params;
}

PostprocessorAdaptiveDT::PostprocessorAdaptiveDT(const std::string & name, InputParameters parameters) :
  TimeStepper(name, parameters),
  PostprocessorInterface(parameters),
    _on_trial(false),
    _pps_good(true),
    _increase_rate(getParam<Real>("increase_rate")),
    _decrease_rate(getParam<Real>("decrease_rate")),
    _update_step(getParam<unsigned int>("update_step")),
    _pps_value(getPostprocessorValue("postprocessor")),
    _max_record_size(getParam<size_t>("max_record_size")),
    _u_saved(NULL),
    _u_older_saved(NULL),
    _aux_saved(NULL),
    _aux_older_saved(NULL),
    _counter(0)
{}

void
PostprocessorAdaptiveDT::preExecute()
{
  TimeStepper::preExecute();
  TransientNonlinearImplicitSystem & nl_sys = _fe_problem.getNonlinearSystem().sys();
  _u_saved = &nl_sys.add_vector("u_saved", false, GHOSTED);
  _u_older_saved = &nl_sys.add_vector("u_older_saved", false, GHOSTED);

  TransientExplicitSystem & aux_sys = _fe_problem.getAuxiliarySystem().sys();
  _aux_saved = &aux_sys.add_vector("aux_saved", false, GHOSTED);
  _aux_older_saved = &aux_sys.add_vector("aux_older_saved", false, GHOSTED);
}

void
PostprocessorAdaptiveDT::preSolve()
{
  TransientNonlinearImplicitSystem & nl_sys = _fe_problem.getNonlinearSystem().sys();
  TransientExplicitSystem & aux_sys = _fe_problem.getAuxiliarySystem().sys();

  // save solution vectors
  *_u_saved = *nl_sys.current_local_solution;
  *_u_older_saved = *nl_sys.older_local_solution;

  *_aux_saved = *aux_sys.current_local_solution;
  *_aux_older_saved = *aux_sys.older_local_solution;

  _u_saved->close();
  _u_older_saved->close();
  _aux_saved->close();
  _aux_older_saved->close();
}

Real
PostprocessorAdaptiveDT::computeInitialDT()
{
  _pps_record.clear();
  _pps_diff_record.clear();
  if(_t_step>1) _pps_record.push_back(_pps_value);
  _pps_good=true;
  return getParam<Real>("dt");
}

Real
PostprocessorAdaptiveDT::computeDT()
{
  //We assume that if we enter computeDT(), converge() return true;
  if(!_update_step) return getCurrentDT();
  if(_on_trial){
    //previous trial is successful.
    _on_trial=false;
    _counter=1;
  }else{
    ++_counter;
  }

  if(_t_step>2){
    _pps_diff_record.push_back(_pps_value-_pps_record.back()); //pps_record should not be empty if we arrive here
    Moose::out<<"_pps_diff:\n";
    for(size_t i=0;i<_pps_diff_record.size();i++) {
	     Moose::out<<_pps_diff_record[i]<<" ";
     }
    Moose::out<<"\n";
  }
  Moose::out<<"In computeDT, _t_step="<<_t_step<<":"<<_pps_value<<"\n";
  _pps_record.push_back(_pps_value);
  if(_pps_record.size()>_max_record_size){
    _pps_record.pop_front();
    _pps_diff_record.pop_front();
  }
  Real new_dt;

  if(_counter % _update_step==0){
    _old_counter=_counter;
    _old_dt=getCurrentDT();
    new_dt=getCurrentDT()*_increase_rate;
    Moose::out<<"PostprocessorAdaptiveDT: Try to enlarge time step to "<<new_dt<<std::endl;
    _counter=0;
    _on_trial=true;
  }else new_dt=getCurrentDT();
  return new_dt;
}

void
PostprocessorAdaptiveDT::rejectStep()
{
  if(!_converged) Moose::out<<"PostprocessorAdaptiveDT:solver not converged"<<std::endl;
  else {
    Moose::out<<"PostprocessorAdaptiveDT: pps value is bad: "<<_pps_value<<std::endl;
    Moose::out << "PostprocessorAdaptiveDT: Marking last solve not converged"<< std::endl;
  }
  TransientNonlinearImplicitSystem & nl_sys = _fe_problem.getNonlinearSystem().sys();
  TransientExplicitSystem & aux_sys = _fe_problem.getAuxiliarySystem().sys();

  // recover initial state
  *nl_sys.current_local_solution = *_u_saved;
  *nl_sys.old_local_solution = *_u_saved;
  *nl_sys.older_local_solution = *_u_older_saved;

  *aux_sys.current_local_solution = *_aux_saved;
  *aux_sys.old_local_solution = *_aux_saved;
  *aux_sys.older_local_solution = *_aux_older_saved;

  nl_sys.solution->close();
  nl_sys.old_local_solution->close();
  nl_sys.older_local_solution->close();
  aux_sys.solution->close();
  aux_sys.old_local_solution->close();
  aux_sys.older_local_solution->close();
  TimeStepper::rejectStep();
}

Real
PostprocessorAdaptiveDT::computeFailedDT(){
   if (_dt <= _dt_min)
    mooseError("Solve failed and timestep already at or below dtmin, cannot continue!");

   Real new_dt;

   if(_on_trial){
     new_dt=_old_dt;
     _counter=_old_counter;
     _on_trial=false;
     Moose::out<<"trial step fail: back to old dt:"<<std::endl;
   }else{
      Moose::out<<"shrink step size:"<<std::endl;
     new_dt=getCurrentDT()*_decrease_rate;
     _counter=0;
   }
   _pps_good=true;
   Moose::out<<"new step size:"<<new_dt<<std::endl;
   return new_dt;
}

bool
PostprocessorAdaptiveDT::converged()
{
  return _converged &&_pps_good;
}

void
PostprocessorAdaptiveDT::postSolve(){
  if(_t_step<2) return;
  if(_converged){
     Moose::out<<"In postSolve, _t_step="<<_t_step<<":"<<_pps_value<<"\n";
    PostprocessorValue diff=_pps_value-_pps_record.back();
    if(diff>0.0)
     _pps_good=false; //Change me back
    //pps_good=true
    else{
      //_pps_good=true;
      if(!_on_trial) _pps_good=true;
      else{
	_pps_good=false;
	//PostprocessorValue max_diff=-1.0*std::numeric_limits<PostprocessorValue>::max();
	// PostprocessorValue max_diff=-1e7;
      	//  size_t n=_pps_diff_record.size();
      	//  if(n>0){
      	//    for(size_t i=0;i<n;i++) {
	//      Moose::out<<_pps_diff_record[i]<<" ";
	//      max_diff=std::max(max_diff,_pps_diff_record[i]);
	//    }
	//    Moose::out<<std::endl;
	//    Moose::out<<"PostprocessorAdaptiveDT: diff="<<diff<<" "<<"history max="<<max_diff<<std::endl;
      	//    if(diff<max_diff) _pps_good=true;
      	//    else {
	//      _pps_good=false;
	//    }
      	//  }else _pps_good=true;
      }
    }
  }
}
