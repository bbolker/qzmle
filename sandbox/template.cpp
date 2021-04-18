#include <TMB.hpp> 
template<class Type> 
Type objective_function<Type>::operator() () { 
DATA_VECTOR(killed); 
DATA_VECTOR(density); 
DATA_MATRIX(X_logit_a); 
DATA_SPARSE_MATRIX(Z_logit_a); 
PARAMETER_VECTOR(logit_a_param); 
PARAMETER(log_h); 
PARAMETER_VECTOR(logit_a_rand); 
PARAMETER(logit_a_logsd); 
vector <Type> logit_a = X_logit_a * logit_a_param; 
logit_a += exp(logit_a_logsd) * (Z_logit_a * logit_a_rand);
vector <Type> a = invlogit(logit_a); 
Type h = exp(log_h); 
Type nll = 0.0;
nll -= sum(dbinom(killed, density,vector<Type>(a/(1 + a * h * density)), true));
nll -= sum(dnorm(logit_a_rand, Type(0), Type(1), true));

return nll;}

