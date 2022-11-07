##### Make Significance Funnel with User-Specified Pooled Estimates #####

# compute meta-analytic effect sizes for an example dataset
require(metafor)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)

# flip signs since we think publication bias operates in favor of negative
# effects alternatively, if not flipping signs, could pass
# favor_positive = FALSE to significance_funnel
dat$yi <- -dat$yi

# optional: regular meta-analysis of all studies (for the black diamond)
# for flexibility, you can use any choice of meta-analysis model here
# in this case, we'll use the robust independent specification since the
# point estimates seem to be from unique papers
# thus, each study gets its own studynum
require(robumeta)
meta_all <- robumeta::robu(yi ~ 1,
                           studynum = 1:nrow(dat),
                           data = dat,
                           var.eff.size = vi,
                           small = TRUE)

# optional: calculate worst-case estimate (for the gray diamond) by analyzing
# only the nonaffirmative studies
dat$pval <- 2 * (1 - pnorm(abs(dat$yi / sqrt(dat$vi))))  # two-tailed p-value
dat$affirm <- (dat$yi > 0) & (dat$pval < 0.05)  # is study affirmative?
meta_worst <- robumeta::robu(yi ~ 1,
                             studynum = 1:nrow( dat[ dat$affirm == TRUE, ] ),
                             data = dat[ dat$affirm == TRUE, ],
                             var.eff.size = vi,
                             small = TRUE)

##### Make Significance Funnel with Alpha = 0.50 and Default Pooled Estimates #####
# change alpha to 0.50 just for illustration
# now the pooled estimates are from the fixed-effect specification because
# they are not provided by the user
significance_funnel(yi = dat$yi,
                    vi = dat$vi,
                    favor_positive = TRUE,
                    alpha_select = 0.50,
                    plot_pooled = TRUE)
