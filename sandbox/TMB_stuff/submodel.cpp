#include <TMB.hpp>
template<class Type>
Type objective_function <Type>::operator() () {
  DATA_VECTOR(density);
  DATA_VECTOR(surv);
  DATA_MATRIX(X_log_a);
  PARAMETER_VECTOR(log_a_param);
  PARAMETER(h);

  vector<Type> log_a = X_log_a*log_a_param;
  //vector<Type> a = exp(log_a);
  vector<Type> prob = exp(log_a)/(1 + exp(log_a)*h*density);
  Type nll = 0.0;
  //nll = -sum(dbinom_robust(surv, density, prob, true));
  for (int i=0; i<surv.size(); i++) {
       nll -= dbinom_robust(surv[i], density[i], prob[i], true);
   }
  return nll;
}
