/**
 * @file   PostprocessorAdaptiveDT.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Fri Nov 15 13:17:57 2013
 *
 * @brief  Use pps for adaptive time stepping
 *
 *
 */
#ifndef PostprocessorAdaptiveDT_H_
#define PostprocessorAdaptiveDT_H_

#include "TimeStepper.h"
#include "PostprocessorInterface.h"
#include "libmesh/numeric_vector.h"

#include<deque>
class PostprocessorAdaptiveDT;

template<>
InputParameters validParams<PostprocessorAdaptiveDT>();

/**
 *
 */
class PostprocessorAdaptiveDT :
  public TimeStepper,
  public PostprocessorInterface
{
public:
  PostprocessorAdaptiveDT(const std::string & name, InputParameters parameters);
  virtual void preExecute();
  virtual void preSolve();
  virtual void postSolve();
  virtual void rejectStep();
  virtual bool converged();

protected:
  virtual Real computeInitialDT();
  virtual Real computeDT();
  virtual Real computeFailedDT();
  bool _on_trial, _pps_good;
  Real _increase_rate, _decrease_rate;
  unsigned int _update_step;
  const PostprocessorValue & _pps_value;
  //history record of pps
  std::deque<PostprocessorValue> _pps_record;
  std::deque<PostprocessorValue> _pps_diff_record;
  size_t _max_record_size;
  //for saving states
  NumericVector<Number> * _u_saved, * _u_older_saved;
  NumericVector<Number> * _aux_saved, * _aux_older_saved;
  size_t _counter, _old_counter;
  Real _old_dt;
};

#endif /* PostprocessorAdaptiveDT_H_ */
