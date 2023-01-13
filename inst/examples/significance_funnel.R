##### Make Significance Funnel #####
# compute meta-analytic effect sizes for an example dataset
require(metafor)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)

# favor_positive = FALSE since we think publication bias operates in favor of negative
significance_funnel(yi = dat$yi, vi = dat$vi, favor_positive = FALSE)
