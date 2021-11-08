#include <TMB.hpp> 
template<class Type> 
Type objective_function<Type>::operator() () { 
DATA_VECTOR(y); 
DATA_VECTOR(x); 
DATA_MATRIX(X_lymax); 
DATA_MATRIX(X_lhalf); 
PARAMETER_VECTOR(lymax_param); 
PARAMETER_VECTOR(lhalf_param); 
vector <Type> lymax = X_lymax * lymax_param; 
vector <Type> lhalf = X_lhalf * lhalf_param; 
Type nll = 0.0;
nll -= sum(dpois(y, vector<Type>(exp(lymax)/(1 + x/exp(lhalf))), true));

return nll;}

