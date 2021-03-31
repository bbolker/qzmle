#include <TMB.hpp>
template<class Type>
Type objective_function <Type>::operator() () {
  DATA_VECTOR(density);
  DATA_VECTOR(surv);
  DATA_MATRIX(X_log_a);
  PARAMETER_VECTOR(log_a_param);
  PARAMETER(h);


  //vector<Type> log_a = X_log_a*log_a_param;
  //vector<Type> prob = exp(log_a)/(1 + exp(log_a)*h*density);
  //ADREPORT(prob);
  Type nll = 0.0;
  vector<Type> log_a (surv.size());
  //nll -= dbinom_robust(surv, density, vector<Type>(exp(log_a)/(1 + exp(log_a)*h*density)), true).sum();
  for (int i=0; i<surv.size(); i++) {
    log_a(i) = (vector<Type>(X_log_a.row(i))*log_a_param).sum();
    nll -= dbinom_robust(surv[i], density[i], log_a(i)/(1 + log_a(i)*h*density[i]), true);
   }
  return nll;
}
