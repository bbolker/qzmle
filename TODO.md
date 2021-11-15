* ideal behaviour for unspecified method: default is R. If not missing() and explicitly specified as "R", then warn when random effects are in the model that we're switching to TMB

* figure out why re_rand_val breaks stuff when we're running a TMB model *without* REs (and fix it)

* for TMB models, maybe `names(coefficients) <- gsub("_param\\.", "", names(coefficients))` (after everything is fitted) ? (ideally, we'd like to have TMB and R output be as close to identical as possible)

* try out some examples from lme4 and glmmTMB with the models specified explicitly

e.g.

lmer(Reaction ~ Days + (1|Subject), sleepstudy)  
mle(Reaction ~ dnorm(mean, exp(logsd)),
   parameters = list(mean ~ Days + (1|Subject)))

