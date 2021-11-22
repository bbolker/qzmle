#include <TMB.hpp> 
template<class Type> 
Type objective_function<Type>::operator() () { 
DATA_VECTOR(killed); 
DATA_VECTOR(density); 
DATA_MATRIX(X_logit_a); 
PARAMETER_VECTOR(logit_a_param); 
PARAMETER(log_h); 
vector <Type> logit_a = X_logit_a * logit_a_param; 
vector <Type> a = invlogit(logit_a); 
Type h = exp(log_h); 
Type nll = 0.0;
nll -= sum(dbinom(killed, density,vector<Type>(a/(1 + a * h * density)), true));

return nll;}

