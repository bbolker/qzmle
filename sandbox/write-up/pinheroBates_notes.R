data(Oats)
str(Oats)
plot(Oats)
## 6 blocks, yields as function of nitrogen cont.


## linear model - hard to interpret
model=lm(yield~Variety*nitro, data=Oats)
summary(model)

## Intercept - Golden rain yield when nitrogen = 0
## VarietyM/V - differ from goldenrain (nitrogen = 0)
## nitro - (intercept) nitrogen effect on golden rain
## VarietyM/V::nitro - differ golden rain

data(Oats)
library(nlme)
model2=lme(yield~Variety*nitro, ## fixed effects
           data=Oats,
           random=~1|Block/Variety) ## random effects (nested)

summary(model2)
## parameter estimate dont change but std.err is impacted
coef(model)
coef(model2)

## every block has a different intercept
## mean of the variety intercept == group intercept of the model
mean(coef(model2)[[1]])
coef(model)[[1]]

## visualize that the random effect are randomly distributed
plot(ranef(model2))
plot(model2)

model3=lme(yield~Variety*nitro, ## fixed effects
           data=Oats,
           random=~ 1|Block)



