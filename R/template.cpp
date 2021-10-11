#include <TMB.hpp> 
template<class Type> 
Type objective_function<Type>::operator() () { 
DATA_VECTOR(y); 
DATA_VECTOR(x); 
PARAMETER(b0); 
PARAMETER(b1); 
PARAMETER(log_sigma); 
Type sigma = exp(log_sigma); 
Type nll = 0.0;
nll -= sum(dnorm(y, b0 + b1 * x,log_sigma, true));

return nll;}

