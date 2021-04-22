#include <TMB.hpp>
template<class Type>
Type objective_function <Type>::operator() () {

  DATA_VECTOR(x1);

  PARAMETER_VECTOR(u) //random

  T
