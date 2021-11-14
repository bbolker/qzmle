#include <TMB.hpp> 
template<class Type> 
Type objective_function<Type>::operator() () { 
DATA_VECTOR(surv); 
DATA_VECTOR(density); 
DATA_MATRIX(X_log_a); 
PARAMETER_VECTOR(log_a_param); 
PARAMETER(h); 
vector <Type> log_a = X_log_a * log_a_param; 
vector <Type> a = exp(log_a); 
Type nll = 0.0;
nll -= sum(dbinom(surv, density,vector<Type>(a/(1 + a * h * density)), true));

return nll;}

