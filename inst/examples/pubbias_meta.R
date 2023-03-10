# calculate effect sizes from example dataset in metafor
require(metafor)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)

# first fit fixed-effects model without any bias correction
# since the point estimate is negative here, we'll assume publication bias
# favors negative log-RRs rather than positive ones
metafor::rma(yi, vi, data = dat, method = "FE")

# warmup
# note that passing selection_ratio = 1 (no publication bias) yields the naive
# point estimate from rma above, which makes sense
meta <- pubbias_meta(yi = dat$yi,
                     vi = dat$vi,
                     selection_ratio = 1,
                     model_type = "fixed",
                     favor_positive = FALSE)
summary(meta)

# assume a known selection ratio of 5
# i.e., affirmative results are 5x more likely to be published than
# nonaffirmative ones
meta <- pubbias_meta(yi = dat$yi,
                     vi = dat$vi,
                     selection_ratio = 5,
                     model_type = "fixed",
                     favor_positive = FALSE)
summary(meta)

# same selection ratio, but now account for heterogeneity and clustering via
# robust specification
meta <- pubbias_meta(yi = dat$yi,
                     vi = dat$vi,
                     cluster = dat$author,
                     selection_ratio = 5,
                     model_type = "robust",
                     favor_positive = FALSE)
summary(meta)

##### Make sensitivity plot as in Mathur & VanderWeele (2020) #####
# range of parameters to try (more dense at the very small ones)
selection_ratios <- c(200, 150, 100, 50, 40, 30, 20, seq(15, 1))

# compute estimate for each value of selection_ratio
estimates <- lapply(selection_ratios, function(e) {
  pubbias_meta(yi = dat$yi, vi = dat$vi, cluster = dat$author,
               selection_ratio = e, model_type = "robust",
               favor_positive = FALSE)$stats
})
estimates <- dplyr::bind_rows(estimates)
estimates$selection_ratio <- selection_ratios

require(ggplot2)
ggplot(estimates, aes(x = selection_ratio, y = estimate)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "gray") +
  geom_line(lwd = 1.2) +
  labs(x = bquote(eta), y = bquote(hat(mu)[eta])) +
  theme_classic()
