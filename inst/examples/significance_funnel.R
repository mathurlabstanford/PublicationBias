##### Make Significance Funnel with User-Specified Pooled Estimates #####
# compute meta-analytic effect sizes for an example dataset
require(metafor)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)

# favor_positive = FALSE since we think publication bias operates in favor of negative
significance_funnel(yi = dat$yi, vi = dat$vi, favor_positive = FALSE)

# by default, the meta-analytic estimates plotted are from metafor::rma.uni()
# you can instead supply estimates from another meta-analysis model
# for example, using the robust independent specification since the
# point estimates seem to be from unique papers
# require(robumeta)
# meta_all <- robumeta::robu(yi ~ 1,
#                            studynum = 1:nrow(dat),
#                            data = dat,
#                            var.eff.size = vi,
#                            small = TRUE)
#
# # worst-case meta-analysis (non-affirmative studies)
# dat$pval <- 2 * (1 - pnorm(abs(dat$yi / sqrt(dat$vi))))  # two-tailed p-value
# dat$affirm <- (dat$yi > 0) & (dat$pval < 0.05)  # is study affirmative?
# meta_worst <- robumeta::robu(yi ~ 1,
#                              studynum = 1:nrow( dat[ dat$affirm == TRUE, ] ),
#                              data = dat[ dat$affirm == TRUE, ],
#                              var.eff.size = vi,
#                              small = TRUE)
#
# meta_all <- pubbias_meta(yi = dat$yi, vi = dat$vi, selection_ratio = 1) # etc
# meta_all$fits$robust$b.r
# meta_all$data
# pubbias_svalue(yi = dat$yi, vi = dat$vi, favor_positive = FALSE, return_worst_meta = TRUE)
#
# significance_funnel(yi = dat$yi, vi = dat$vi, est_all = meta_all$b.r,
#                     est_worse = meta_worst$b.r, favor_positive = FALSE)
