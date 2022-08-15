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
    N = c(40, 80, 400, 800),
    b1 = c(0, 0.5),
    N_ratio = c(1, .6),
    lambdax = c(.7, 1)
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
    grand_mean <- dalphax * DESIGNFACTOR$N_ratio / (1 + DESIGNFACTOR$N_ratio)
    w <- rbind(1, DESIGNFACTOR$N_ratio)
    vars <- (c(1, 1 + dpsix)^2 + outer(c(0, dalphax), grand_mean, FUN = "-")^2)
    pooled_var <- colSums(vars * w) / colSums(w)
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

with(FIXED, {
    lam1 <- sapply(DESIGNFACTOR$lambdax, FUN = seq_lam, num_items = num_items)
    lam2 <- apply(lam1, MARGIN = 2, add_vec, pos = c(2, 5),
                  dx = dlambdax[[1]])
    tvar1 <- colSums(lam1)^2
    tvar2 <- colSums(lam2)^2
    th1 <- sqrt(thetax + lambdax - lam1^2)
    th2 <- apply(th1, MARGIN = 2, add_vec, pos = c(4, 6),
                 dx = dthetax[[1]])
    psi2 <- dpsix + 1  # standard deviation
    evar1 <- colSums(th1^2) + sum(ucov1) - sum(diag(ucov1))
    evar2 <- colSums(th2^2) + sum(ucov2) - sum(diag(ucov2))
    DESIGNFACTOR$omega1 <<- tvar1 / (tvar1 + evar1)
    DESIGNFACTOR$omega2 <<- tvar2 / (tvar2 + evar2)
})

# Data Generation ---------------------------------------------------------

GenData <- function(condition, fixed_objects = NULL) {
    num_items <- fixed_objects$num_items
    num_obs <- condition$N * 2 * c(1, condition$N_ratio) / (1 + condition$N_ratio)
    b1 <- condition$b1
    lam1 <- seq_lam(condition$lambdax, num_items = num_items)
    lam2 <- add_vec(lam1, c(2, 5), fixed_objects$dlambdax[[1]])
    nu1 <- rep(fixed_objects$nux, num_items)
    nu2 <- add_vec(nu1, c(4, 5), fixed_objects$dnux[[1]])
    th1 <- sqrt(fixed_objects$thetax + condition$lambdax - lam1^2)
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
# test_dat <- GenData(DESIGNFACTOR[7, ], fixed_objects = FIXED)

# Analysis function and subfunctions --------------------------------------

ExtractMx <- function(model, par = "b1", wald_ci = FALSE) {
    se <- mxSE(par, model, silent = TRUE)
    if (wald_ci) {
        est = model@output$algebras[[par]]
        c(est = c(est), se = c(se),
          ll = c(est - qnorm(.975) * se),
          ul = c(est + qnorm(.975) * se))
    } else {
        ci <- model@output$confidenceIntervals[par, ]
        c(est = ci[["estimate"]], se = c(se),
          ll = ci[["lbound"]], ul = ci[["ubound"]])
    }
}

AnalyseRegMx <- function(condition, dat, fixed_objects = NULL) {
    b1 <- condition$beta1
    sd_xsum <- sd(dat$xsum)
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    reg_mx <- mxModel("REG",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("xsum", "y", "group"),
        latentVars = "fx",
        mxAlgebra(1 - vg * a1^2, name = "evfx"),
        # Factor loadings
        mxPath(from = "fx", to = "xsum", free = TRUE,
               values = sd_xsum),
        # Path
        mxPath(from = "group", to = c("fx", "y"), free = TRUE,
               values = 0, labels = c("a1", "a2")),
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "group", arrows = 2, free = FALSE, values = vg,
               labels = "vg"),
        mxPath(from = "fx", arrows = 2, free = FALSE, labels = "evfx[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE, values = 1),
        # Unique variances
        # Mean
        mxPath(from = "one", to = "group", values = 1.5, free = FALSE),
        mxPath(from = "one", to = "xsum", values = 0, free = TRUE),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "b0"),
        mxCI(c("b1"))
    )
    reg_fit <- mxRun(reg_mx, intervals = TRUE, silent = TRUE)
    if (!reg_fit@output$status$status == 0) {
        stop("Reg did not converge")
    }
    out <- ExtractMx(reg_fit)
    if (anyNA(out)) {
        stop("Reg did not obtain SE or CI")
    }
    out
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu2_labs),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a2"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha2[1,1]")
    )
    mgsem_mx <- mxModel(
        "MultipleGroupSEM",
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
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        num_items = fixed_objects$num_items
    )
    if (!sem_fit@output$status$status == 0) {
        stop("full sem did not converge")
    }
    out <- ExtractMx(sem_fit)
    if (anyNA(out)) {
        stop("full sem did not obtain SE or CI")
    }
    out
}
# Test: Full SEM Analyses
# AnalyseSem(DESIGNFACTOR[1, ], test_dat, fixed_objects = FIXED)

RunSemMx2 <- function(dat, lambdax, dlambdax, b1, num_items = 6) {
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
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam1, labels = lam1_labs),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi1"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu1_labs),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a1"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha1")
    )
    sem2_mx <- mxModel("SEM2",
        type = "RAM",
        mxData(observed = dat[dat$group == 2, ], type = "raw"),
        manifestVars = c(ind_names, "y"),
        latentVars = "fx",
        # Factor loadings
        mxPath(from = "fx", to = ind_names, free = TRUE,
               values = lam2, labels = lam2_labs),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = TRUE,
               values = 1, labels = "psi2"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Unique variances
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu2_labs),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a2"),
        mxPath(from = "one", to = c("fx"), values = 0, free = TRUE,
               labels = "alpha2")
    )
    mgsem_mx <- mxModel(
        "MultipleGroupSEM",
        mxMatrix("Full",
            nrow = 1, ncol = 2, name = "rel_n",
            values = num_obs / sum(num_obs), free = FALSE
        ),
        sem1_mx, sem2_mx,
        mxFitFunctionMultigroup(c("SEM1", "SEM2")),
        mxAlgebra((rel_n[1, 1] * alpha1 + rel_n[1, 2] * alpha2) / sum(rel_n),
                  name = "alpha"),
        mxAlgebra(
            (rel_n[1, 1] * (psi1 + (alpha1 - alpha)^2) +
                rel_n[1, 2] * (psi2 + (alpha2 - alpha)^2)) / sum(rel_n),
            name = "psi"
        ),
        mxAlgebra(b1 * sqrt(psi), name = "beta1")
    )
    mxRun(mgsem_mx, intervals = TRUE, silent = TRUE)
}
# Test: Full SEM (Mx)
# RunSemMx2(test_dat, lambdax = 1, dlambdax = .5, b1 = 0.5)

AnalyseSem2 <- function(condition, dat, fixed_objects = NULL) {
    sem_fit <- RunSemMx2(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        num_items = fixed_objects$num_items
    )
    if (!sem_fit@output$status$status == 0) {
        stop("full sem2 did not converge")
    }
    out <- ExtractMx(sem_fit, par = "MultipleGroupSEM.beta1",
                     wald_ci = TRUE)
    if (anyNA(out)) {
        stop("full sem2 did not obtain SE or CI")
    }
    out
}
# Test: Full SEM Analyses (standardization afterwards)
# AnalyseSem2(DESIGNFACTOR[1, ], test_dat, fixed_objects = FIXED)

GetFscoresMx <- function(dat, lambdax, dlambdax,
                         unit_weight = FALSE, return_croon = FALSE,
                         num_items = 6) {
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu1_labs),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha1[1,1]"),
        # Composite reliability
        mxAlgebra(sum(CFA1.A[1:6, 7]), name = "sum_lam1"),
        mxAlgebra(sum_lam1^2 * psi1, name = "tvar1"),
        mxAlgebra(sum(CFA1.S[1:6, 1:6]), name = "evar1"),
        mxAlgebra(tvar1 / (tvar1 + evar1), name = "omega1")
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
        # Mean
        mxPath(from = "one", to = ind_names, values = 0, free = TRUE,
               labels = nu2_labs),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha2[1,1]"),
        # Composite reliability
        mxAlgebra(sum(CFA2.A[1:6, 7]), name = "sum_lam2"),
        mxAlgebra(sum_lam2^2 * psi2, name = "tvar2"),
        mxAlgebra(sum(CFA2.S[1:6, 1:6]), name = "evar2"),
        mxAlgebra(tvar2 / (tvar2 + evar2), name = "omega2")
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
    if (unit_weight) {
        out <- rbind(
            cbind(
                dat[dat$group == 1, "xsum", drop = FALSE] / num_items,
                rel = mgcfa_fit$output$algebras$CFA1.omega1,
                fs_ld = mgcfa_fit$output$algebras$CFA1.sum_lam1
            ),
            cbind(
                dat[dat$group == 2, "xsum", drop = FALSE] / num_items,
                rel = mgcfa_fit$output$algebras$CFA2.omega2,
                fs_ld = mgcfa_fit$output$algebras$CFA2.sum_lam2
            )
        )
        names(out)[1] <- "fx_fs"
    } else {
        fscore <- mxFactorScores(mgcfa_fit, type = "Regression")
        psi_names <- c("CFA1.psi1", "CFA2.psi2")
        lambda_names <- c("CFA1.A", "CFA2.A")
        fit_names <- c("CFA1.fitfunction", "CFA2.fitfunction")
        out <- do.call(
            rbind,
            lapply(seq_along(fscore), function(i) {
                psi <- mgcfa_fit$output$algebras[[psi_names[i]]]
                lambda_vec <- mgcfa_fit$output$matrices[[lambda_names[i]]][
                    seq_len(num_items),
                    num_items + 1
                ]
                implied_cov <- attr(
                    mgcfa_fit$output$algebras[[fit_names[i]]],
                    "expCov"
                )
                fs_mat <- solve(implied_cov, lambda_vec %*% psi)
                Alam <- crossprod(fs_mat, lambda_vec)
                data.frame(
                    fx_fs = c(fscore[[i]][, , "Scores"]),
                    rel = 1 - fscore[[i]][, , "StandardErrors"]^2 /
                        c(psi),
                    rel2 = c(Alam %*% psi %*% t(Alam)) /
                        c(crossprod(fs_mat, cov(dat[dat$group == i, ind_names])) %*%
                            fs_mat)
                    # rel = var(fscore[[i]][, , "Scores"]) /
                    #     c(mgcfa_fit$output$algebras[[psi_names[i]]])
                )
            })
        )
        # Compute factor score matrix (same as reliability)
        if (return_croon) {
            # lambda_names <- c("CFA1.A", "CFA2.A")
            theta_names <- c("CFA1.S", "CFA2.S")
            # fit_names <- c("CFA1.fitfunction", "CFA2.fitfunction")
            out_fsm <- lapply(seq_along(fscore), function(i) {
                psi <- mgcfa_fit$output$algebras[[psi_names[i]]]
                lambda_vec <- mgcfa_fit$output$matrices[[lambda_names[i]]][
                    seq_len(num_items),
                    num_items + 1
                ]
                implied_cov <- attr(
                    mgcfa_fit$output$algebras[[fit_names[i]]],
                    "expCov"
                )
                theta_mat <- mgcfa_fit$output$matrices[[theta_names[i]]][
                    seq_len(num_items),
                    seq_len(num_items)
                ]
                fs_mat <- solve(implied_cov, lambda_vec %*% psi)
                list(
                    Alam = crossprod(fs_mat, lambda_vec),
                    AThetaAt = crossprod(fs_mat, theta_mat) %*% fs_mat
                )
            })
            out <- list(fs = out, croon = out_fsm)
        }
    }
    out
}
# Test: Compute factor scores (Mx)
# GetFscoresMx(test_dat, lambdax = .53, dlambdax = 0.3)

GetFscoresBartMx <- function(dat, lambdax, dlambdax, num_items = 6) {
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
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
        mxPath(from = ind_names, arrows = 2, free = TRUE, values = 2.56,
               lbound = 0.001),
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
    fscore <- mxFactorScores(mgcfa_fit, type = "ML")
    psi_names <- c("CFA1.psi1", "CFA2.psi2")
    lambda_names <- c("CFA1.A", "CFA2.A")
    theta_names <- c("CFA1.S", "CFA2.S")
    do.call(
        rbind,
        lapply(seq_along(fscore), function(i) {
            psi <- mgcfa_fit$output$algebras[[psi_names[i]]]
            lambda_vec <- mgcfa_fit$output$matrices[[lambda_names[i]]][
                seq_len(num_items),
                num_items + 1
            ]
            theta_mat <- mgcfa_fit$output$matrices[[theta_names[i]]][
                seq_len(num_items),
                seq_len(num_items)
            ]
            lam_thinv <- solve(theta_mat, lambda_vec)
            fs_mat <- solve(crossprod(lambda_vec, lam_thinv), t(lam_thinv))
            data.frame(
                fx_fs = c(fscore[[i]][, , "Scores"]),
                rel = c(psi) /
                    c(fscore[[i]][, , "StandardErrors"]^2),
                rel2 = c(psi) /
                    c(fs_mat %*% cov(dat[dat$group == i, ind_names]) %*% t(fs_mat))
            )
        })
    )
}
# Test: Compute factor scores (Mx)
# GetFscoresBartMx(test_dat, lambdax = .53, dlambdax = 0.3)

RunFspa <- function(dat, lambdax, dlambdax, b1, vg, num_items = 6) {
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        num_items = num_items
    )
    dat <- cbind(dat, fs)
    form_mx <- mxModel("FSPA",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fx_fs", "y", "group"),
        latentVars = "fx",
        mxAlgebra(1 - vg * a1^2, name = "evfx"),
        # Factor loadings
        mxPath(from = "fx", to = "fx_fs", free = TRUE,
               values = sqrt(mean(dat$rel))),
        # Path
        mxPath(from = "group", to = c("fx", "y"), free = TRUE,
               values = 0, labels = c("a1", "a2")),
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "group", arrows = 2, free = FALSE, values = vg,
               labels = "vg"),
        mxPath(from = "fx", arrows = 2, free = FALSE, labels = "evfx[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE, values = 1),
        # Unique variances
        # Mean
        mxPath(from = "one", to = "group", values = 1.5, free = FALSE),
        mxPath(from = "one", to = "fx_fs", values = 0, free = TRUE),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "b0"),
        mxCI(c("b1"))
    )
    mxRun(form_mx, intervals = TRUE, silent = TRUE)
}
# Test: FS-PA
# RunFspa(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5)

AnalyseFspa <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    fspa_fit <- RunFspa(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items
    )
    if (!fspa_fit@output$status$status == 0) {
        stop("fspa did not converge")
    }
    out <- ExtractMx(fspa_fit)
    if (anyNA(out)) {
        stop("fspa did not obtain SE or CI")
    }
    out
}
# Test: FS-PA Analyses
# AnalyseFspa(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

RunCroon <- function(dat, lambdax, dlambdax, b1, num_items = 6) {
    num_obs <- as.numeric(table(dat$group))
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        return_croon = TRUE,
        num_items = num_items
    )
    dat <- cbind(dat, fs$fs)
    cov_y_xi <- lapply(1:2, function(g) {
        out <- cov(dat[dat$group == g, c("y", "fx_fs")])
        out[1, 2] <- out[2, 1] <- out[2, 1] / as.numeric(fs$croon[[g]]$Alam)
        out[2, 2] <- (out[2, 2] - fs$croon[[g]]$AThetaAt) /
        as.numeric(crossprod(fs$croon[[g]]$Alam))
        out
    })
    pa1_mx <- mxModel("PA1",
        type = "RAM",
        mxData(cov_y_xi[[1]], type = "cov",
               means = colMeans(dat[dat$group == 1, c("y", "fx_fs")]),
               numObs = sum(dat$group == 1)),
        manifestVars = c("y", "fx_fs"),
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
        mxPath(from = "fx", to = "fx_fs", free = TRUE,
               values = 1, labels = "lam"),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi1[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "fx_fs", values = 0, free = TRUE,
               labels = "nu"),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a1"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha1[1,1]")
    )
    pa2_mx <- mxModel("PA2",
        type = "RAM",
        mxData(cov_y_xi[[2]], type = "cov",
               means = colMeans(dat[dat$group == 2, c("y", "fx_fs")]),
               numObs = sum(dat$group == 2)),
        manifestVars = c("y", "fx_fs"),
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
        mxPath(from = "fx", to = "fx_fs", free = TRUE,
               values = 1, labels = "lam"),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "psi2[1,1]"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "fx_fs", values = 0, free = TRUE,
               labels = "nu"),
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a2"),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE,
               labels = "alpha2[1,1]")
    )
    mgpa_mx <- mxModel(
        "MultiplePathAnalysis",
        pa1_mx, pa2_mx,
        mxFitFunctionMultigroup(c("PA1", "PA2")),
        mxCI(c("b1"))
    )
    mxRun(mgpa_mx, intervals = TRUE, silent = TRUE)
}
# Test: FS-PA with Croon
# RunCroon(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5)

AnalyseCroon <- function(condition, dat, fixed_objects = NULL) {
    croon_fit <- RunCroon(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        num_items = fixed_objects$num_items
    )
    if (!croon_fit@output$status$status == 0) {
        stop("croon did not converge")
    }
    out <- ExtractMx(croon_fit)
    if (anyNA(out)) {
        stop("croon did not obtain SE or CI")
    }
    out
}
# Test: FS-PA with Croon Analyses
# AnalyseCroon(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

RunCroon2 <- function(dat, lambdax, dlambdax, b1, num_items = 6) {
    num_obs <- as.numeric(table(dat$group))
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        return_croon = TRUE,
        num_items = num_items
    )
    dat <- cbind(dat, fs$fs)
    dat$fx <- dat$fx_fs
    cov_y_xi <- lapply(1:2, function(g) {
        out <- cov(dat[dat$group == g, c("y", "fx")])
        out[1, 2] <- out[2, 1] <- out[2, 1] / as.numeric(fs$croon[[g]]$Alam)
        out[2, 2] <- (out[2, 2] - fs$croon[[g]]$AThetaAt) /
        as.numeric(crossprod(fs$croon[[g]]$Alam))
        out
    })
    pa1_mx <- mxModel("PA1",
        type = "RAM",
        mxData(cov_y_xi[[1]], type = "cov",
               means = colMeans(dat[dat$group == 1, c("y", "fx")]),
               numObs = sum(dat$group == 1)),
        manifestVars = c("y", "fx"),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = TRUE,
               values = 1, labels = "psi1"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a1"),
        mxPath(from = "one", to = c("fx"), values = 0, free = TRUE,
               labels = "alpha1")
    )
    pa2_mx <- mxModel("PA2",
        type = "RAM",
        mxData(cov_y_xi[[2]], type = "cov",
               means = colMeans(dat[dat$group == 2, c("y", "fx")]),
               numObs = sum(dat$group == 2)),
        manifestVars = c("y", "fx"),
        # Path
        mxPath(from = "fx", to = "y", free = TRUE,
               values = b1, labels = "b1"),
        # Variance
        mxPath(from = "fx", arrows = 2, free = TRUE,
               values = 1, labels = "psi2"),
        mxPath(from = "y", arrows = 2, free = TRUE,
               values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "y", values = 0, free = TRUE,
               labels = "a2"),
        mxPath(from = "one", to = c("fx"), values = 0, free = TRUE,
               labels = "alpha2")
    )
    mgpa_mx <- mxModel(
        "MultiplePathAnalysis",
        mxMatrix("Full",
            nrow = 1, ncol = 2, name = "rel_n",
            values = num_obs / sum(num_obs), free = FALSE
        ),
        pa1_mx, pa2_mx,
        mxFitFunctionMultigroup(c("PA1", "PA2")),
        mxAlgebra((rel_n[1, 1] * alpha1 + rel_n[1, 2] * alpha2) / sum(rel_n),
                  name = "alpha"),
        mxAlgebra(
            (rel_n[1, 1] * (psi1 + (alpha1 - alpha)^2) +
                rel_n[1, 2] * (psi2 + (alpha2 - alpha)^2)) / sum(rel_n),
            name = "psi"
        ),
        mxAlgebra(b1 * sqrt(psi), name = "beta1"),
        mxCI("beta1")
    )
    mxRun(mgpa_mx, intervals = TRUE, silent = TRUE)
}
# Test: FS-PA with Croon
# RunCroon2(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5)

AnalyseCroon2 <- function(condition, dat, fixed_objects = NULL) {
    num_obs <- as.numeric(table(dat$group))
    total_n <- sum(num_obs)
    fs <- GetFscoresMx(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        return_croon = TRUE,
        num_items = fixed_objects$num_items
    )
    dat <- cbind(dat, fs$fs)
    dat$fx <- dat$fx_fs
    cov_y_fxfs <- lapply(1:2, function(g) {
        cov(dat[dat$group == g, c("y", "fx_fs")])
    })
    cov_y_xi <- lapply(1:2, function(g) {
        out <- cov(dat[dat$group == g, c("y", "fx")])
        out[1, 2] <- out[2, 1] <- out[2, 1] / as.numeric(fs$croon[[g]]$Alam)
        out[2, 2] <- (out[2, 2] - fs$croon[[g]]$AThetaAt) /
        as.numeric(crossprod(fs$croon[[g]]$Alam))
        out
    })
    cov_y_fxfs_pooled <- (cov_y_fxfs[[1]] * (num_obs[1] - 1) +
        cov_y_fxfs[[2]] * (num_obs[2] - 1)) / (total_n - 2)
    cov_y_xi_pooled <- (cov_y_xi[[1]] * (num_obs[1] - 1) +
        cov_y_xi[[2]] * (num_obs[2] - 1)) / (total_n - 2)
    est_gamma <- cov_y_xi_pooled[2, 1] / cov_y_xi_pooled[2, 2]
    err_var <- cov_y_xi_pooled[1, 1] - est_gamma ^ 2 * cov_y_fxfs_pooled[2, 2]
    est_se <- sqrt(err_var / cov_y_fxfs_pooled[2, 2] / (total_n - 3))
    c(est = est_gamma, se = est_se,
      ll = est_gamma - qt(.975, df = total_n - 3) * est_se,
      ul = est_gamma + qt(.975, df = total_n - 3) * est_se)
}
# Test: FS-PA with Croon Analyses (standardization afterwards)
# AnalyseCroon2(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

Run2spa <- function(dat, lambdax, dlambdax, b1, vg, unit_weight = FALSE,
                    num_items = 6, emp_rel = FALSE) {
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        num_items = num_items, unit_weight = unit_weight
    )
    if (emp_rel) {
        fs$rel <- fs$rel2
    }
    fs$rel_ld <- fs$rel / fs$rel[1]
    dat <- cbind(dat, fs)
    form_mx <- mxModel("2SPA",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fx_fs", "y", "group"),
        latentVars = "fx",
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "ld1",
            values = 1, free = TRUE
        ),
        mxAlgebra(1 - vg * a1^2, name = "ev_fx"),
        mxAlgebra(ld1 * data.rel_ld, name = "ld"),
        mxAlgebra(ld^2 / data.rel * (1 - data.rel), name = "ev_fs"),
        # Factor loadings
        mxPath(from = "fx", to = "fx_fs", free = FALSE,
               values = sqrt(mean(dat$rel)), labels = "ld[1,1]"),
        # Error variance
        mxPath(from = "group", arrows = 2, free = FALSE, values = vg,
               labels = "vg"),
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "ev_fx[1,1]"),
        mxPath(from = "fx_fs", arrows = 2, free = FALSE,
               values = 1 - mean(dat$rel), labels = "ev_fs[1,1]"),
        mxPath(from = "group", to = c("fx", "y"), free = TRUE,
               values = 0, labels = c("a1", "a2")),
        mxPath(from = "fx", to = "y", values = b1, labels = "b1"),
        mxPath(from = "y", arrows = 2, values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "group", values = 1.5, free = FALSE),
        mxPath(from = "one", to = c("fx_fs", "y"), values = 0, free = TRUE),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE),
        mxCI(c("b1"))
    )
    mxRun(form_mx, intervals = TRUE, silent = TRUE)
}
# Test: 2S-PA
# Run2spa(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5, vg = .25)

Analyse2spa <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    tspa_fit <- Run2spa(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spa did not converge")
    }
    out <- ExtractMx(tspa_fit)
    if (anyNA(out)) {
        stop("2spa did not obtain SE or CI")
    }
    out
}
# Test: Analyse function for 2S-PA
# Analyse2spa(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

Analyse2spaEmp <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    tspa_fit <- Run2spa(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items,
        emp_rel = TRUE
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spa did not converge")
    }
    out <- ExtractMx(tspa_fit)
    if (anyNA(out)) {
        stop("2spa did not obtain SE or CI")
    }
    out
}
# Test: Analyse function for 2S-PA (empirical reliability)
# Analyse2spaEmp(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

Run2spaBart <- function(dat, lambdax, dlambdax, b1, vg,
                        num_items = 6, emp_rel = FALSE) {
    fs <- GetFscoresBartMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        num_items = num_items
    )
    if (emp_rel) {
        fs$rel <- fs$rel2
    }
    dat <- cbind(dat, fs)
    # if (any(dat$rel <= 0)) {
    #     dat$rel <- ifelse(dat$rel <= 0, yes = 0.01, no = dat$rel)
    #     warning("Negative reliability estimates")
    # }
    form_mx <- mxModel("2SPA",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fx_fs", "y", "group"),
        latentVars = c("fx"),
        mxAlgebra(1 - vg * a1^2, name = "ev_fx"),
        mxAlgebra(ld^2 / data.rel * (1 - data.rel), name = "ev_fs"),
        # Factor loadings
        mxPath(from = "fx", to = "fx_fs", free = FALSE,
               values = 1, labels = "ld"),
        # Error variance
        mxPath(from = "group", arrows = 2, free = FALSE, values = vg,
               labels = "vg"),
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "ev_fx[1,1]"),
        mxPath(from = "fx_fs", arrows = 2, free = FALSE,
               values = 1 - mean(dat$rel), labels = "ev_fs[1,1]"),
        mxPath(from = "group", to = c("fx", "y"), free = TRUE,
               values = 0, labels = c("a1", "a2")),
        mxPath(from = "fx", to = "y", values = b1, labels = "b1"),
        mxPath(from = "y", arrows = 2, values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "group", values = 1.5, free = FALSE),
        mxPath(from = "one", to = c("fx_fs", "y"), values = 0, free = TRUE),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE),
        mxCI(c("b1"))
    )
    mxRun(form_mx, intervals = TRUE, silent = TRUE)
}
# Test: 2S-PA (Bartlett)
# Run2spaBart(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5, vg = .25)

Analyse2spaBart <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    tspa_fit <- Run2spaBart(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spa-Bart did not converge")
    }
    out <- ExtractMx(tspa_fit)
    if (anyNA(out)) {
        stop("2spa-Bart did not obtain SE or CI")
    }
    out
}
# Test: Analyse function for 2S-PA
# Analyse2spaBart(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

Analyse2spaBartEmp <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    tspa_fit <- Run2spaBart(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items,
        emp_rel = TRUE
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spa-Bart did not converge")
    }
    out <- ExtractMx(tspa_fit)
    if (anyNA(out)) {
        stop("2spa-Bart did not obtain SE or CI")
    }
    out
}
# Test: Analyse function for 2S-PA
# Analyse2spaBartEmp(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

Run2spaUW <- function(dat, lambdax, dlambdax, b1, vg,
                      num_items = 6) {
    fs <- GetFscoresMx(dat,
        lambdax = lambdax, dlambdax = dlambdax,
        num_items = num_items, unit_weight = TRUE
    )
    fs$rel_ld <- fs$fs_ld / fs$fs_ld[1]
    dat <- cbind(dat, fs)
    # if (any(dat$rel <= 0)) {
    #     dat$rel <- ifelse(dat$rel <= 0, yes = 0.01, no = dat$rel)
    #     warning("Negative reliability estimates")
    # }
    form_mx <- mxModel("2SPA",
        type = "RAM",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fx_fs", "y", "group"),
        latentVars = c("fx"),
        mxMatrix("Full",
            nrow = 1, ncol = 1, name = "ld1",
            values = 1, free = TRUE
        ),
        mxAlgebra(1 - vg * a1^2, name = "evfx"),
        mxAlgebra(ld1 * data.rel_ld, name = "ld"),
        mxAlgebra(ld^2 / data.rel * (1 - data.rel), name = "ev_fs"),
        # Factor loadings
        mxPath(from = "fx", to = "fx_fs", free = FALSE,
               values = sqrt(mean(dat$rel)), labels = "ld[1,1]"),
        # Error variance
        mxPath(from = "group", arrows = 2, free = FALSE, values = vg,
               labels = "vg"),
        mxPath(from = "fx", arrows = 2, free = FALSE,
               values = 1, labels = "evfx[1,1]"),
        mxPath(from = "fx_fs", arrows = 2, free = FALSE,
               values = 1 - mean(dat$rel), labels = "ev_fs[1,1]"),
        mxPath(from = "group", to = c("fx", "y"), free = TRUE,
               values = 0, labels = c("a1", "a2")),
        mxPath(from = "fx", to = "y", values = b1, labels = "b1"),
        mxPath(from = "y", arrows = 2, values = 1 - b1^2),
        # Mean
        mxPath(from = "one", to = "group", values = 1.5, free = FALSE),
        mxPath(from = "one", to = c("fx_fs", "y"), values = 0, free = TRUE),
        mxPath(from = "one", to = c("fx"), values = 0, free = FALSE),
        mxCI(c("b1"))
    )
    mxRun(form_mx, intervals = TRUE, silent = TRUE)
}
# Test: 2S-PA
# Run2spaUW(test_dat, lambdax = .53, dlambdax = 0.3, b1 = 0.5)

Analyse2spaUW <- function(condition, dat, fixed_objects = NULL) {
    vg <- condition$N_ratio / (condition$N_ratio + 1) -
        (condition$N_ratio / (condition$N_ratio + 1))^2
    tspa_fit <- Run2spaUW(dat,
        lambdax = condition$lambdax,
        dlambdax = fixed_objects$dlambdax,
        b1 = condition$beta1,
        vg = vg,
        num_items = fixed_objects$num_items
    )
    if (!tspa_fit@output$status$status == 0) {
        stop("2spauw did not converge")
    }
    out <- ExtractMx(tspa_fit)
    if (anyNA(out)) {
        stop("2spauw did not obtain SE or CI")
    }
    out
}
# Test: Analyse function for 2S-PA (with unit weight)
# Analyse2spaUW(DESIGNFACTOR[1, ], dat = test_dat, fixed_objects = FIXED)

# Evaluate ----------------------------------------------------------------

Evaluate <- function(condition, results, fixed_objects = NULL) {
    b1_pop <- condition$beta1
    meth <- c("reg", "fspa", "croon2", "sem2",
              "2spa", "2spae", "2spabart", "2spabarte", "2spauw")
    est <- results[, paste0(meth, ".est")]
    se <- results[, paste0(meth, ".se")]
    ci <- results[, t(outer(meth, Y = c(".ll", ".ul"), FUN = paste0))]
    c(
        bias = bias(est, b1_pop),
        rmse = RMSE(est, b1_pop),
        rsb = bias(
            sweep(se,
                MARGIN = 2, STATS = apply(est, 2, sd), FUN = "/"),
            parameter = 1,
            type = "relative"
        ),
        coverage = ECR(
            ci,
            parameter = b1_pop
        )
    )
}

res <- runSimulation(
    design = DESIGNFACTOR,
    replications = 2500,
    generate = GenData,
    analyse = list(reg = AnalyseRegMx,
                   fspa = AnalyseFspa,
                #    croon = AnalyseCroon,
                   croon2 = AnalyseCroon2,
                #    sem = AnalyseSem,
                   sem2 = AnalyseSem2,
                   `2spa` = Analyse2spa,
                   `2spae` = Analyse2spaEmp,
                   `2spabart` = Analyse2spaBart,
                   `2spabarte` = Analyse2spaBartEmp,
                   `2spauw` = Analyse2spaUW),
    summarise = Evaluate,
    fixed_objects = FIXED,
    packages = c("OpenMx", "psych"),
    filename = "simulation/error-in-x-results-emp",
    seed = rep(198300, nrow(DESIGNFACTOR)),
    parallel = TRUE,
    ncores = 20L,
    save = TRUE,
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "error-in-x-results-rev_"
    )
)
