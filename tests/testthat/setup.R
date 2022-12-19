library(metafor)

# helper fn for simulating meta-analysis with publication bias

# # p: row of parameters dataframe
sim_data <- function(p) {
  withr::with_seed(42, {

    nn <- p$k * p$per.cluster

    # generate cluster random intercepts
    gam1 <- rnorm(n = p$k, mean = 0, sd = sqrt(p$V.gam))
    gam1i <- rep(gam1, each = p$per.cluster)

    # generate individual-study random intercepts
    gam2i <- rnorm(n = nn, mean = 0, sd = sqrt(p$V - p$V.gam))

    # individual study means
    mui <- p$mu + gam1i + gam2i
    sei <- runif(n = nn, min = p$sei.min, max = p$sei.max)
    yi <- rnorm(n = nn, mean = mui, sd = sei)

    d <- data.frame(cluster = rep(1:p$k, each = p$per.cluster),
                    Study.name = 1:nn,
                    yi = yi,
                    sei = sei,
                    vi = sei^2,
                    pval = 2 * (1 - pnorm(abs(yi) / sei)))

    # 1-tailed publication bias
    signif <- d$pval < 0.05 && d$yi > 0
    publish <- rep(1, nrow(d))
    publish[!signif] <- rbinom(n = sum(!signif), size = 1,
                               prob = 1 / p$selection_ratio)

    d$weight <- 1
    d$weight[signif == 0] <- p$selection_ratio
    d <- d[publish == 1, ]

    return(d)
  })
}
