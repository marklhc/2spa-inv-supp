# Simulation code for comparing 2S-PA and SEM in latent regression
# with a fallible predictor

library(MASS)
library(SimDesign)
library(psych)
library(OpenMx)


set.seed(1221)

# Note: the standard error in standardized coefficients assume
#       that the latent predictor is stochastic

# Design factors:
DESIGNFACTOR <- createDesign(
    N = c(50, 100, 1000),
    b1 = c(0, 0.5)
)

# Function to generate covariance matrix for minor factors
get_ucov <- function(p, scale = sqrt(3.56 * 0.1), seed = 201026, n = 5) {
    seed_state <- .GlobalEnv$.Random.seed
    set.seed(seed)
    W <- matrix(rnorm(p * n), nrow = n)
    .GlobalEnv$.Random.seed <- seed_state
    WtW <- crossprod(W)
    D <- diag(1 / sqrt(diag(WtW))) * scale
    D %*% WtW %*% D
}

FIXED <- list(
    N_ratio = .6,
    num_items = 6,
    lambdax = 1,  # standardized loadings
    dlambdax = .5,
    nux = 0,
    dnux = list(c(-.80, .50)),
    dalphax = -0.3,
    dpsix = 0.2,
    thetax = 2.56,
    dthetax = list(c(0.25, -0.5))
)
FIXED <- within(FIXED, {
    ucov1 <- get_ucov(max(FIXED$num_items))
    ucov2 <- get_ucov(max(FIXED$num_items))
})

DESIGNFACTOR$beta1 <- with(FIXED, {
    grand_mean <- dalphax * N_ratio / (1 + N_ratio)
    pooled_var <- weighted.mean(
        c(1, 1 + dpsix)^2 +
            (c(0, dalphax) - grand_mean)^2,
        w = c(1, N_ratio)
    )
    DESIGNFACTOR$b1 * sqrt(pooled_var)
})

# Data Generation ---------------------------------------------------------

# Helper functions
seq_lam <- function(lambda, num_items) {
    seq(1.25, 0.75, length.out = num_items) * lambda
}

add_vec <- function(x, pos, dx) {
    x[pos] <- x[pos] + dx
    x
}

GenData <- function(condition, fixed_objects = NULL) {
    num_items <- fixed_objects$num_items
    num_obs <- condition$N * c(1, fixed_objects$N_ratio)
    b1 <- condition$b1
    lam1 <- seq_lam(fixed_objects$lambdax, num_items = num_items)
    lam2 <- add_vec(lam1, c(2, 5), fixed_objects$dlambdax[[1]])
    nu1 <- rep(fixed_objects$nux, num_items)
    nu2 <- add_vec(nu1, c(4, 5), fixed_objects$dnux[[1]])
    th1 <- sqrt(fixed_objects$thetax + fixed_objects$lambdax - lam1^2)
    th2 <- add_vec(th1, c(4, 6), fixed_objects$dthetax[[1]])
    psi2 <- fixed_objects$dpsix + 1  # standard deviation
    alpha2 <- fixed_objects$dalphax
    etax1 <- rnorm(num_obs[1])
    etax2 <- rnorm(num_obs[2],
        mean = alpha2,
        sd = psi2
    )
    Theta1 <- diag(th1^2 - diag(fixed_objects$ucov1)) +
        fixed_objects$ucov1
    Theta2 <- diag(th2^2 - diag(fixed_objects$ucov2)) +
        fixed_objects$ucov2
    x1 <- tcrossprod(lam1, etax1) + nu1 +
        t(mvrnorm(num_obs[1],
            mu = rep(0, num_items),
            Sigma = Theta1
        ))
    x2 <- tcrossprod(lam2, etax2) + nu2 +
        t(mvrnorm(num_obs[2],
            mu = rep(0, num_items),
            Sigma = Theta1
        ))
    y <- c(etax1, etax2) * b1 +
        rnorm(sum(num_obs), sd = sqrt(1 - b1^2))
    df <- data.frame(t(cbind(x1, x2)), c(etax1, etax2), y,
        group = rep(1:2, num_obs)
    )
    colnames(df) <- c(paste0("x", seq_len(num_items)), "etax", "y", "group")
    df$xsum <- rowSums(t(cbind(x1, x2)))
    df
}
# Test: generate data
# test_dat <- GenData(DESIGNFACTOR[1, ], fixed_objects = FIXED)

# Analysis function and subfunctions --------------------------------------

ExtractMx <- function(model, par = "b1") {
    se <- mxSE(par, model, silent = TRUE)
    ci <- model@output$confidenceIntervals[par, ]
    c(est = ci[["estimate"]], se = c(se),
      ll = ci[["lbound"]], ul = ci[["ubound"]])
}

AnalyseRegMx <- function(condition, dat, fixed_objects = NULL) {
    b1 <- condition$beta1
    sd_xsum <- sd(dat$xsum)
    reg_mx <- mxModel("REG",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("xsum", "y"),
        latentVars = "fx",
        # Factor loadings
        mxPath(from = "fx", to = "xsum", free = TRUE,
               values = sd_xsum),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE, values = 1),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Unique variances
        # Mean
        mxPath(from = "one", to = "xsum", values = 0, free = TRUE),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a1"),
        mxCI(c("b1"))
    )
    reg_fit <- mxRun(reg_mx, intervals = TRUE, silent = TRUE)
    if (!reg_fit@output$status$status == 0) {
        stop("Reg did not converge")
    }
    ExtractMx(reg_fit)
}
# Test: Composite score
# AnalyseRegMx(DESIGNFACTOR[1, ], test_dat, fixed_objects = FIXED)

RunSemMx <- function(dat, lambdax, dlambdax, b1, num_items = 6) {
    item_seq <- seq_len(num_items)
    ind_names <- paste0("x", item_seq)
    num_obs <- as.numeric(table(dat$group))
    lam1 <- seq_lam(lambdax, num_items = num_items)
    lam2 <- add_vec(lam1, c(2, 5), dlambdax[[1]])
    lam1_labs <- lam2_labs <- paste0("l", item_seq)
    lam2_labs[c(2, 5)] <- paste0(lam2_labs[c(2, 5)], "2")
    nu1_labs <- nu2_labs <- paste0("n", item_seq)
    nu2_labs[c(2, 4, 5)] <- paste0(nu2_labs[c(2, 4, 5)], "2")
    sem1_mx <- mxModel("SEM1",
        type = "RAM",
        mxData(observed = dat[dat$group == 1, ], type = "raw"),
        manifestVars = c(ind_names, "y"),
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "rel_n1",
            values = num_obs[1] / sum(num_obs), free = FALSE
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dalpha",
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dpsi",
        ),
        mxAlgebra(dalpha / rel_n1, name = "alpha1"),
        mxAlgebra(dpsi / rel_n1 - alpha1^2, name = "psi1"),
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam1, labels = lam1_labs),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi1[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu1_labs),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a1"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha1[1,1]")
    )
    sem2_mx <- mxModel("SEM2",
        type = "RAM",
        mxData(observed = dat[dat$group == 2, ], type = "raw"),
        manifestVars = c(ind_names, "y"),
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "rel_n2",
            values = num_obs[2] / sum(num_obs), free = FALSE
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dalpha",
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dpsi",
        ),
        mxAlgebra(- dalpha / rel_n2, name = "alpha2"),
        mxAlgebra((1 - dpsi) / rel_n2 - alpha2^2, name = "psi2"),
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam2, labels = lam2_labs),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi2[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu2_labs),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a2"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha2[1,1]")
    )
    mgsem_mx <- mxModel(
        "MultipleGroupCFA",
        sem1_mx, sem2_mx,
        mxFitFunctionMultigroup(c("SEM1", "SEM2")),
        mxCI(c("b1"))
    )
    mxRun(mgsem_mx, intervals = TRUE, silent = TRUE)
}
# Test: Full SEM (Mx)
# RunSemMx(test_dat, lambdax = 1, dlambdax = .5, b1 = 0.5)

AnalyseSem <- function(condition, dat, fixed_objects = NULL) {
    sem_fit <- RunSemMx(dat,
        lambdax = fixed_objects$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        num_items = fixed_objects$num_items
    )
    if (!sem_fit@output$status$status == 0) {
        stop("full sem did not converge")
    }
    ExtractMx(sem_fit)
}
# Test: Full SEM Analyses
# AnalyseSem(DESIGNFACTOR[1, ],
#     test_dat,
#     fixed_objects = FIXED
# )

GetFscoresMx <- function(dat, lambdax, dlambdax, num_items = 6) {
    item_seq <- seq_len(num_items)
    ind_names <- paste0("x", item_seq)
    num_obs <- as.numeric(table(dat$group))
    lam1 <- seq_lam(lambdax, num_items = num_items)
    lam2 <- add_vec(lam1, c(2, 5), dlambdax[[1]])
    lam1_labs <- lam2_labs <- paste0("l", item_seq)
    lam2_labs[c(2, 5)] <- paste0(lam2_labs[c(2, 5)], "2")
    nu1_labs <- nu2_labs <- paste0("n", item_seq)
    nu2_labs[c(2, 4, 5)] <- paste0(nu2_labs[c(2, 4, 5)], "2")
    cfa1_mx <- mxModel("CFA1",
        type = "RAM",
        mxData(observed = dat[dat$group == 1, ], type = "raw"),
        manifestVars = ind_names,
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "rel_n1",
            values = num_obs[1] / sum(num_obs), free = FALSE
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dalpha",
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dpsi",
        ),
        mxAlgebra(dalpha / rel_n1, name = "alpha1"),
        mxAlgebra(dpsi / rel_n1 - alpha1^2, name = "psi1"),
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam1, labels = lam1_labs),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi1[1,1]"),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu1_labs),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha1[1,1]")
    )
    cfa2_mx <- mxModel("CFA2",
        type = "RAM",
        mxData(observed = dat[dat$group == 2, ], type = "raw"),
        manifestVars = ind_names,
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "rel_n2",
            values = num_obs[2] / sum(num_obs), free = FALSE
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dalpha",
        ),
        mxMatrix("Full",
            nrow = 1, ncol = 1,
            values = 0, free = TRUE, labels = "dpsi",
        ),
        mxAlgebra(- dalpha / rel_n2, name = "alpha2"),
        mxAlgebra((1 - dpsi) / rel_n2 - alpha2^2, name = "psi2"),
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam2, labels = lam2_labs),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi2[1,1]"),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu2_labs),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha2[1,1]")
    )
    mgcfa_mx <- mxModel(
        "MultipleGroupCFA",
        cfa1_mx, cfa2_mx,
        mxFitFunctionMultigroup(c("CFA1", "CFA2"))
    )
    mgcfa_fit <- mxRun(mgcfa_mx, silent = TRUE)
    if (!mgcfa_fit@output$status$status == 0) {
        stop("CFA model not converged")
    }
    fscore <- mxFactorScores(mgcfa_fit, type = "Regression")
    psi_names <- c("CFA1.psi1", "CFA2.psi2")
    do.call(
        rbind,
        lapply(seq_along(fscore), function(i) {
            data.frame(
                fx_fs = c(fscore[[i]][, , "Scores"]),
                rel = 1 - fscore[[i]][, , "StandardErrors"]^2 /
                    c(mgcfa_fit$output$algebras[[psi_names[i]]])
            )
        })
    )
}
# Test: Compute factor scores (Mx)
# GetFscoresMx(test_dat, lambdax = .53, dlambdax = 0.3)

Run2spa <- function(dat, lambdax, dlambdax, b1, num_items = 6) {
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        num_items = num_items
    )
    dat <- cbind(dat, fs)
    if (any(dat$rel <= 0)) {
        dat$rel <- ifelse(dat$rel <= 0, yes = 0.01, no = dat$rel)
        warning("Negative reliability estimates")
    }
    form_mx <- mxModel("2SPA",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fx_fs", "y"),
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "v_fs",
            values = 1, free = TRUE
        ),
        mxAlgebra(sqrt(v_fs * data.rel), name = "ld"),
        mxAlgebra(v_fs * (1 - data.rel), name = "ev_fs"),
        # Factor loadings
        mxPath(from = "fx", to = "fx_fs", free = FALSE,
               values = sqrt(mean(dat$rel)), labels = "ld[1,1]"),
        # Error variance
        mxPath(from = "fx", arrows = 2, free = FALSE, values = 1),
        mxPath(from = "fx_fs", arrows = 2, free = FALSE,
               values = 1 - mean(dat$rel), labels = "ev_fs[1,1]"),
        mxPath(from = "fx", to = "y", values = b1, labels = "b1"),
        mxPath(from = "y", arrows = 2, values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = c("fx_fs", "y"), values = 0, free = TRUE),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE),
        mxCI(c("b1"))
    )
    mxRun(form_mx, intervals = TRUE, silent = TRUE)
}
# Test: 2S-PA
# Run2spa(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5)

Analyse2spa <- function(condition, dat, fixed_objects = NULL) {
    tspa_fit <- Run2spa(dat,
        lambdax = fixed_objects$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        num_items = fixed_objects$num_items
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spa did not converge")
    }
    ExtractMx(tspa_fit)
}
# Test: Analyse function
# Analyse2spa(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

# Evaluate ----------------------------------------------------------------

Evaluate <- function(condition, results, fixed_objects = NULL) {
    b1_pop <- condition$beta1
    est <- results[, c("reg.est", "sem.est", "2spa.est")]
    c(
        bias = bias(est, b1_pop),
        rmse = RMSE(est, b1_pop),
        rsb = bias(
            sweep(results[, c("reg.se", "sem.se", "2spa.se")],
                MARGIN = 2, STATS = apply(est, 2, sd), FUN = "/"),
            parameter = 1,
            type = "relative"
        ),
        coverage = ECR(
            results[, c(
                "reg.ll", "reg.ul",
                "sem.ll", "sem.ul",
                "2spa.ll", "2spa.ul"
            )],
            parameter = b1_pop
        )
    )
}

res <- runSimulation(
    design = DESIGNFACTOR,
    replications = 2500,
    generate = GenData,
    analyse = list(reg = AnalyseRegMx,
                   sem = AnalyseSem,
                   `2spa` = Analyse2spa),
    summarise = Evaluate,
    fixed_objects = FIXED,
    packages = c("OpenMx", "psych"),
    filename = "simulation/error-in-x-results",
    seed = rep(15300812, nrow(DESIGNFACTOR)),
    parallel = TRUE,
    ncores = 10L,
    save = TRUE,
    # save_seeds = TRUE,
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "error-in-x-results_"
    )
)
