
  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() () { DATA_VECTOR(x); DATA_VECTOR(y); PARAMETER(b0); PARAMETER(b1); Type nll = 0.0; nll = -sum(dnorm(+, true)); return nll; }

  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() () { DATA_VECTOR(x); DATA_VECTOR(y); PARAMETER(b0); PARAMETER(b1); Type nll = 0.0; nll = -sum(dnorm(b0, true)); return nll; }

  #include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() () { DATA_VECTOR(x); DATA_VECTOR(y); PARAMETER(b0); PARAMETER(b1); Type nll = 0.0; nll = -sum(dnorm(b1 * x, true)); return nll; }
