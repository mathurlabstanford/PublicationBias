# calculate effect sizes from example dataset in metafor
require(metafor)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)

##### Fixed-Effects Specification #####
# S-values and worst-case meta-analysis under fixed-effects specification
svals_fixed_0 <- pubbias_svalue(yi = dat$yi,
                                vi = dat$vi,
                                q = 0,
                                model_type = "fixed",
                                favor_positive = FALSE)

# publication bias required to shift point estimate to 0
svals_fixed_0$stats$sval_est

# and to shift CI to include 0
svals_fixed_0$stats$sval_ci

# now try shifting to a nonzero value (RR = 0.90)
svals_fixed_q <- pubbias_svalue(yi = dat$yi,
                                vi = dat$vi,
                                q = log(.9),
                                model_type = "fixed",
                                favor_positive = FALSE)

# publication bias required to shift point estimate to RR = 0.90
svals_fixed_q$stats$sval_est

# and to shift CI to RR = 0.90
svals_fixed_q$stats$sval_ci

##### Robust Clustered Specification #####
svals <- pubbias_svalue(yi = dat$yi,
                        vi = dat$vi,
                        q = 0,
                        model_type = "robust",
                        favor_positive = FALSE)
summary(svals)
