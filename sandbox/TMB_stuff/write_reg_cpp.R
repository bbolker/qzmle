tmb_model <- "
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

DATA_VECTOR(x);
DATA_VECTOR(y);

PARAMETER(b0);
PARAMETER(b1);
//PARAMETER(log_sigma);

//Type sigma = exp(log_sigma);

Type nll = 0.0;

nll = -sum(dnorm(y, b0 + b1 * x, sd, true)); //sigma

return nll;

}"
write(tmb_model, file = "reg.cpp")
