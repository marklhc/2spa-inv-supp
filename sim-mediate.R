library(MASS)
library(SimDesign)
library(lavaan)
library(OpenMx)
library(mirt)

set.seed(1221)

# Design factors:
DESIGNFACTOR <- createDesign(
    N = c(50, 100, 300),
    beta1 = c(0, 0.5),  # treat --> fm
    beta2 = 0.3,  # treat --> fy
    beta3 = c(0, 0.4)  # fm --> fy
)

# Function to generate covariance matrix for minor factors
get_ucov <- function(p, scale = sqrt(0.1), seed = 201026, n = 26) {
    seed_state <- .GlobalEnv$.Random.seed
    set.seed(seed)
    W <- matrix(rnorm(p * n), nrow = n)
    .GlobalEnv$.Random.seed <- seed_state
    WtW <- crossprod(W)
    D <- diag(1 / sqrt(diag(WtW))) * scale
    D %*% WtW %*% D
}

# Function to compute error variance
compute_psi <- function(beta_mat, fixed_psi) {
    p <- nrow(beta_mat)
    psi_mat <- matrix(0, nrow = p, ncol = p)
    pfixed <- nrow(fixed_psi)
    psi_mat[seq_len(pfixed), seq_len(pfixed)] <- fixed_psi
    IminusBinv <- solve(diag(p) - beta_mat)
    for (i in (pfixed + 1):p) {
        psi_mat[i, i] <- 1 - IminusBinv[i, ] %*% psi_mat %*% IminusBinv[i, ]
    }
    diag(psi_mat)[-seq_len(pfixed)]
}

get_psimy <- function(beta1, beta3) {
    out <- mapply(function(b1, b3) {
        compute_psi(
            beta_mat = rbind(0, 0, c(0.2, b1, 0, 0), c(0.3, 0.3, b3, 0)),
            fixed_psi = matrix(c(0.25, 0, 0, 0.25), nrow = 2)
        )
    },
    b1 = beta1, b3 = beta3,
    SIMPLIFY = FALSE
    )
    out <- do.call(rbind, out)
    colnames(out) <- c("psim", "psiy")
    out
}
DESIGNFACTOR <- cbind(DESIGNFACTOR,
                      with(DESIGNFACTOR, sqrt(get_psimy(beta1, beta3))))

FIXED <- list(
    num_items_y = 20,
    num_items_m = 6,
    lambdam = 1,
    dlambdam = -0.5,
    num = 0,
    dnum = list(c(-.80, .50)),
    thetam = 2.56,
    dthetam = list(c(0.25, -0.5)),
    alpham1 = 0,
    dalpham = 0.2,
    # psim = sqrt(0.9275),
    discrim1 = c(
        1.578, 1.154, 0.958, 1.067, 1.290,
        1.446, 0.791, 1.628, 1.334, 0.942,
        1.236,
        0.816, 0.932, 1.227, 1.179, 1.052,
        1.144, 0.496, 0.692, 0.715
    ) / 1.7,
    diff1 = c(
        -2.118, 0.032, -1.723, -0.466, -1.114,
        -0.700, 0.579, -0.796, 0.474, 0.441,
        -0.820,
        1.061, 0.910, 0.531, 0.035, 0.435,
        1.070, -1.303, -1.097, 0.933
    ),
    ddiscrimy = list(c(0.5, 0.3, -0.2) / 1.7),
    ddiffy = list(c(-0.3, 0.5, -0.3)),
    # psiy = sqrt(0.7666),
    alphay1 = 0,
    dalphay = 0.22
)
FIXED <- within(FIXED, {
    num_items <- num_items_y + num_items_m
    ucov1 <- get_ucov(num_items)
    ucov2 <- get_ucov(num_items, seed = 201027)
})

# Data Generation ---------------------------------------------------------

# Helper functions
seq_lam <- function(lambda, num_items, from = 1.25, to = 0.75) {
    seq(from, to = to, length.out = num_items) * lambda
}

add_vec <- function(x, pos, dx) {
    x[pos] <- x[pos] + dx
    x
}

GenData <- function(condition, fixed_objects = NULL) {
    # Extract parameters from condition/fixed_objects
    num_items_m <- fixed_objects$num_items_m
    num_items_y <- fixed_objects$num_items_y
    num_items <- num_items_m + num_items_y
    num_obs <- rep(condition$N, 2)
    b1 <- condition$beta1
    b2 <- condition$beta2
    b3 <- condition$beta3
    lambdam <- fixed_objects$lambdam
    alpham1 <- fixed_objects$alpham1
    psim <- condition$psim
    alphay1 <- fixed_objects$alphay1
    psiy <- condition$psiy
    discrim1 <- fixed_objects$discrim1
    diff1 <- fixed_objects$diff1
    # Compute population parameters
    lam1 <- seq_lam(lambdam, num_items = num_items_m)
    lam2 <- add_vec(lam1, c(2, 5), fixed_objects$dlambdam[[1]])
    nu1 <- rep(fixed_objects$num, num_items_m)
    nu2 <- add_vec(nu1, c(4, 5), fixed_objects$dnum[[1]])
    th1 <- sqrt(fixed_objects$thetam + lambdam - lam1^2)
    th2 <- add_vec(th1, c(4, 6), fixed_objects$dthetam[[1]])
    alpham2 <- alpham1 + fixed_objects$dalpham
    discrim2 <- add_vec(
        discrim1, c(1, 5, 9), fixed_objects$ddiscrimy[[1]]
    )
    diff2 <- add_vec(
        diff1, c(2, 5, 8), fixed_objects$ddiffy[[1]]
    )
    alphay2 <- alphay1 + fixed_objects$dalphay
    # Simulate errors
    u1 <- mvrnorm(num_obs[1],
        mu = rep(0, num_items),
        Sigma = fixed_objects$ucov1 + diag(0.9, num_items)
    )
    u2 <- mvrnorm(num_obs[2],
        mu = rep(0, num_items),
        Sigma = fixed_objects$ucov2 + diag(0.9, num_items)
    )
    # Simulate group
    group <- rep(1:2, num_obs)
    # Simulate treatment
    treat <- rep(c(0, 1, 0, 1), rep(num_obs / 2, 2))
    # Simulate M
    etam1 <- alpham1 + b1 * treat[group == 1] +
        rnorm(num_obs[1], sd = psim)
    etam2 <- alpham2 + b1 * treat[group == 2] +
        rnorm(num_obs[2], sd = psim)
    m1 <- t(tcrossprod(lam1, etam1) + nu1) +
        u1[, seq_len(num_items_m)] %*% diag(th1)
    m2 <- t(tcrossprod(lam2, etam2) + nu2) +
        u2[, seq_len(num_items_m)] %*% diag(th2)
    # Simulate Y
    etay1 <- alphay1 + b2 * treat[group == 1] + b3 * etam1 +
        rnorm(num_obs[1], sd = psiy)
    etay2 <- alphay2 + b2 * treat[group == 2] + b3 * etam2 +
        rnorm(num_obs[2], sd = psiy)
    p_y1 <- pnorm(discrim1 * -outer(diff1, etay1, "-") +
        t(u1[, seq_len(num_items_y) + num_items_m]))
    y1 <- p_y1
    y1[] <- rbinom(p_y1, size = 1, prob = p_y1)
    p_y2 <- pnorm(discrim2 * -outer(diff2, etay2, "-") +
        t(u2[, seq_len(num_items_y) + num_items_m]))
    y2 <- p_y2
    y2[] <- rbinom(p_y2, size = 1, prob = p_y2)
    df <- data.frame(treat,
        rbind(m1, m2),
        t(cbind(y1, y2)),
        c(etam1, etam2),
        c(etay1, etay2),
        group = group
    )
    colnames(df) <- c(
        "treat",
        paste0("m", seq_len(num_items_m)),
        paste0("y", seq_len(num_items_y)),
        "etam", "etay",
        "group"
    )
    df
}
# Test: generate data
# test_dat <- GenData(DESIGNFACTOR[1,], fixed_objects = FIXED)

RunIrtY <- function(dat, num_items, pars = NULL) {
    ninv_a <- paste(setdiff(seq_len(num_items), c(1, 5, 9)), collapse = ", ")
    ninv_d <- paste(setdiff(seq_len(num_items), c(2, 5, 8)), collapse = ", ")
    # Factor scores adjusting for invariance
    difmodel <- paste0(c(
        paste0("FY = 1-", num_items),
        paste0("CONSTRAINB = (", ninv_a, ", a1),"),
        paste0("             (", ninv_d, ", d)")),
        collapse = "\n"
    )
    ynames <- paste0("y", seq_len(num_items))
    multipleGroup(dat[, ynames],
        model = difmodel,
        group = as.character(dat$group),
        invariance = c("free_mean", "free_var"),
        pars = pars,
        technical = list(NCYCLES = 5000)
    )
}

# Get starting values based on a large data
FIXED$irt_pars <- local({
    large_condition <- cbind(N = 100000, DESIGNFACTOR[1, -1])
    large_dat <- GenData(large_condition, fixed_objects = FIXED)
    mirt::mod2values(RunIrtY(large_dat, num_items = FIXED$num_items_y))
})

GetFscoresM <- function(dat, lambdam, dlambdam, num_items = 6) {
    item_seq <- seq_len(num_items)
    lam1 <- seq_lam(lambdam, num_items = num_items)
    lam2 <- add_vec(lam1, c(2, 5), dlambdam[[1]])
    start_lds <- paste0(
        "start(c(", lam1 * 1.005, ", ", lam2 * 1.005, ")) * m", item_seq,
        collapse = " + "
    )
    # Factor scores adjusting for invariance
    mgcfa_mod <- paste0(
        " fm =~ NA * m1 + ", start_lds, " \n ",
        " fm ~ c(0, NA) * 1
          fm ~~ c(1, NA) * fm "
    )
    mgcfa_fit <- cfa(mgcfa_mod,
        data = dat,
        group = "group",
        group.equal = c("loadings", "intercepts"),
        group.partial = c(
            "fm=~m2", "fm=~m5",
            "m2~1", "m4~1", "m5~1"
        ),
        se = "none",
        test = "none",
        bounds = "standard"
    )
    fscore_m <- lavPredict(mgcfa_fit, se = TRUE)
    do.call(
        rbind,
        lapply(seq_along(fscore_m), function(i) {
            data.frame(
                fm_fs = c(fscore_m[[i]]),
                fm_fs_se = c(attr(fscore_m, "se")[[i]]),
                fm_fs_rel = 1 - c(attr(fscore_m, "se")[[i]]^2 /
                    lavInspect(mgcfa_fit, "est")[[i]]$psi)
            )
        })
    )
}
# Test: Compute factor scores (M)
# GetFscoresM(test_dat, lambdam = 1, dlambdam = -0.5)

GetFscoresY <- function(dat, num_items = 20,
                        pars = NULL) {
    twopl_fit <- RunIrtY(dat,
        num_items = num_items,
        pars = pars
    )
    fscore_y <- mirt::fscores(twopl_fit, full.scores.SE = TRUE)
    psiy <- vapply(coef(twopl_fit), function(x) {
        sqrt(x$GroupPars["par", "COV_11"])
    }, numeric(1))
    out <- cbind(fscore_y, 1 - fscore_y[, 2]^2 / psiy[dat$group])
    colnames(out) <- c("fy_fs", "fy_fs_se", "fy_fs_rel")
    out
}
# Test: Compute factor scores (mirt)
# GetFscoresY(test_dat, pars = FIXED$irt_pars)

AddFscores <- function(condition, fixed_objects = NULL) {
    dat <- GenData(condition, fixed_objects = fixed_objects)
    fsm <- GetFscoresM(dat,
        lambdam = fixed_objects$lambdam,
        dlambdam = fixed_objects$dlambdam,
        num_items = fixed_objects$num_items_m
    )
    fsy <- GetFscoresY(dat,
        num_items = fixed_objects$num_items_y,
        pars = fixed_objects$irt_pars
    )
    cbind(dat, fsm, fsy)
}
# Test: Compute factor scores (mirt)
# test_dat <- AddFscores(DESIGNFACTOR[1,], fixed_objects = FIXED)

# Analysis function and subfunctions --------------------------------------

outer_dot <- function(x, y) {
    outer(x, y, paste, sep = ".")
}

ExtractMx <- function(model,
                      par = c("beta", "std_ind"),
                      parname = c("beta1", "beta2", "beta3", "std_ind")) {
    ci <- model@output$confidenceIntervals
    out_names <- outer_dot(parname, y = c("est", "ll", "ul"))
    out <- cbind(
        est = ci[, "estimate"],
        ll = ci[, "lbound"], ul = ci[, "ubound"]
    )
    setNames(as.vector(out), as.vector(out_names))
}

ExtractLavaan <-
    function(object,
             par = c("beta1", "beta2", "beta3", "std_ind")) {
        pe <- parameterEstimates(object)
        pe <- pe[which(pe$label %in% par), ]
        out_names <- outer_dot(
            par,
            c("est", "se", "ll", "ul")
        )
        out <- cbind(
            est = pe$est, se = pe$se,
            ll = pe$ci.lower,
            ul = pe$ci.upper
        )
        setNames(as.vector(out), as.vector(out_names))
    }

AnalyseSem <- function(condition, dat, fixed_objects = NULL) {
    num_items_m <- fixed_objects$num_items_m
    item_seq_m <- seq_len(num_items_m)
    lam1 <- seq_lam(fixed_objects$lambdam, num_items = num_items_m)
    lam2 <- add_vec(lam1, c(2, 5), fixed_objects$dlambdam[[1]])
    start_lds_m <- paste0(
        "start(c(", lam1 * sqrt(.9275), ", ",
        lam2 * sqrt(.9275), ")) * m", item_seq_m,
        collapse = " + "
    )
    num_items_y <- fixed_objects$num_items_y
    item_seq_y <- seq_len(num_items_y)
    discrim <- subset(fixed_objects$irt_pars, name == "a1",
        select = value, drop = TRUE
    ) / 1.7 * sqrt(.7666)
    start_lds_y <- paste0(
        "start(c(", discrim[item_seq_y], ", ",
        discrim[item_seq_y + num_items_y], ")) * y", item_seq_y,
        collapse = " + "
    )
    sem_mod <- paste0(
        " fm =~ NA * m1 + ", start_lds_m, " \n ",
        " fy =~ NA * y1 + ", start_lds_y, " \n ",
        paste0("y", item_seq_y, " ~~ 1 * y", item_seq_y,
               collapse = " \n "),
        " \n fm ~ c(0, NA) * 1 + c(alpham1, alpham2) * 1
          fm ~~ c(1, NA) * fm + c(psim1, psim2) * fm
          fm ~ c(b1, b1) * treat
          fy ~ c(b3, b3) * fm + c(b2, b2) * treat
          fy ~ c(0, NA) * 1 + c(alphay1, alphay2) * 1
          fy ~~ c(1, NA) * fy + c(psiy1, psiy2) * fy
          # Var
          varfm := (1 + psim2) / 2 + 0.25 * b1^2 + alpham2^2 / 4
          beta1 := b1 / sqrt(varfm)
          varfy := (1 + psiy2) / 2 + 0.25 * b3^2 +
                   ((1 + psim2) / 2 + 0.25 * b1^2) * b2^2 +
                   2 * (1 + psim2) / 2 * b1 * b2 * b3 * 0.25 +
                   (alphay2 + alpham2 * b2)^2 / 4
          beta2 := b2 / sqrt(varfy)
          beta3 := b3 * sqrt(varfm / varfy)
          std_ind := beta1 * beta3 "
    )
    sem_fit <- lavaan::sem(sem_mod,
        data = dat,
        group = "group",
        parameterization = "Theta",
        ordered = paste0("y", item_seq_y),
        group.equal = c("loadings", "intercepts"),
        group.partial = c(
            "fm=~m2", "fm=~m5",
            "xm~1", "m4~1", "m5~1",
            "fy=~y1", "fy=~y5", "fy=~y9",
            "y1|t1", "y2|t1", "y5|t1", "y8|t1", "y9|t1"
        ),
        fixed.x = TRUE,
        bounds = TRUE
    )
    if (!lavInspect(sem_fit, "converged")) {
        stop("SEM did not converge")
    }
    ExtractLavaan(sem_fit)
}
# Test: Full SEM
# AnalyseSem(DESIGNFACTOR[1,], test_dat, fixed_objects = FIXED)

AnalyseFspaMat <- function(condition, dat, fixed_objects = NULL) {
    freeA <- matrix(FALSE, nrow = 4, ncol = 4)
    valA <- 0 * freeA
    freeApos <- matrix(
        byrow = TRUE, ncol = 2,
        c(
            1, 2,
            1, 3,
            1, 4,
            2, 3,
            2, 4
        )
    )
    freeA[freeApos] <- TRUE
    valA[freeApos] <- c(condition$beta3, condition$beta2,
                        fixed_objects$dalphay,
                        condition$beta1, fixed_objects$dalpham)
    freeS <- diag(c(TRUE, TRUE, FALSE, FALSE))
    valS <- 0 * freeS
    diag(valS) <- c(condition$psiy, condition$psim,
                    0.25, 0.25)
    fspa_mx <- mxModel("FSPA",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fy_fs", "fm_fs", "treat", "group"),
        # A matrix
        mxMatrix("Full", nrow = 4, ncol = 4,
                 name = "A",
                 values = valA,
                 free = freeA),
        # S matrix
        mxMatrix("Symm", nrow = 4, ncol = 4,
                 name = "S",
                 values = valS[lower.tri(valS, diag = TRUE)],
                 free = freeS[lower.tri(freeS, diag = TRUE)]
        ),
        # F matrix
        mxMatrix("Iden", nrow = 4, ncol = 4, name = "F"),
        # M matrix (vector)
        mxMatrix("Full", nrow = 1, ncol = 4, name = "M",
                 free = c(TRUE, TRUE, FALSE, FALSE),
                 values = c(0, 0, 0.5, 0.5)),
        # SD
        mxMatrix(type = "Iden", nrow = 4, name = "I"),
        mxAlgebra(solve(I - A)[1:2,], name = "invIminusA"),
        mxAlgebra(rbind(diag2vec(invIminusA %&% S) ^ 0.5,
                        1, 1),
                  name = "Sv"),
        mxAlgebra(vec2diag(Sv ^ -1) %*% A %*% vec2diag(Sv),
                  name = "SA"),
        mxAlgebra(cbind(SA[2,3], SA[1,cbind(3,2)]), name = "beta"),
        mxAlgebra(beta[1,1] * beta[1,3], name = "std_ind"),
        mxExpectationRAM("A", "S", "F", "M",
            dimnames = c("fy_fs", "fm_fs", "treat", "group")
        ),
        mxFitFunctionML(),
        mxCI(c("beta", "std_ind"))
    )
    fspa_fit <- mxRun(fspa_mx, intervals = TRUE, silent = TRUE)
    if (!fspa_fit@output$status$status == 0) {
        stop("FSPA did not converge")
    }
    ExtractMx(fspa_fit)
}
# Test: FS-PA (matrix)
# AnalyseFspaMat(DESIGNFACTOR[1,], test_dat, fixed_objects = FIXED)

Analyse2spaMat <- function(condition, dat, fixed_objects = NULL) {
    freeA <- matrix(FALSE, nrow = 6, ncol = 6)
    valA <- 0 * freeA
    freeApos <- matrix(
        byrow = TRUE, ncol = 2,
        c(
            5, 6,
            5, 3,
            5, 4,
            6, 3,
            6, 4
        )
    )
    freeA[freeApos] <- TRUE
    valA[freeApos] <- c(condition$beta3, condition$beta2,
                        fixed_objects$dalphay,
                        condition$beta1, fixed_objects$dalpham)
    valA[1:2, 5:6] <- diag(2)
    valS <- matrix(0, nrow = 6, ncol = 6)
    diag(valS) <- c(
        mean(dat$fy_fs_rel * (1 - dat$fy_fs_rel)),
        mean(dat$fm_fs_rel * (1 - dat$fm_fs_rel)),
        0.25, 0.25,
        condition$psiy * mean(dat$fy_fs_rel) ^ 2,
        condition$psim * mean(dat$fm_fs_rel) ^ 2
    )
    labS <- NA * valS
    diag(labS)[1:2] <- c("ev_fs[1,1]", "ev_fs[1,2]")
    diag(labS)[5:6] <- c("ev_fy[1,1]", "ev_fm[1,1]")
    valF <- cbind(diag(4), 0, 0)
    tspa_mx <- mxModel("2SPA",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fy_fs", "fm_fs", "treat", "group"),
        latentVars = c("fy", "fm"),
        # Total variances
        mxMatrix(type = "Full", nrow = 1, ncol = 2,
                 values = .7,
                 free = TRUE, name = "v_fs"),
        # Definition variables
        mxMatrix(type = "Full", nrow = 1, ncol = 2, free = FALSE,
                 labels = c("data.fy_fs_rel", "data.fm_fs_rel"),
                 name = "fs_rel"),
        # True and error variances
        mxAlgebra(v_fs * (1 - fs_rel), name = "ev_fs"),
        mxAlgebra(v_fs * fs_rel, name = "vf"),
        # Sub-matrix of S
        mxMatrix("Symm", nrow = 2, ncol = 2,
                 name = "subS1",
                 values = valS[3:4, 3:4][lower.tri(valS[3:4, 3:4], diag = TRUE)],
                 free = FALSE,
                 labels = labS[3:4, 3:4]
        ),
        mxMatrix("Symm", nrow = 3, ncol = 3,
                 name = "subS2",
                 values = valS[c(3:4, 6), c(3:4, 6)][lower.tri(valS[c(3:4, 6), c(3:4, 6)], diag = TRUE)],
                 free = FALSE,
                 labels = labS[c(3:4, 6), c(3:4, 6)]
        ),
        mxAlgebra(vf[1,2] - A[6,3:4] %&% subS1, name = "ev_fm"),
        mxAlgebra(vf[1,1] - A[5,cbind(3,4,6)] %&% subS2, name = "ev_fy"),
        # A matrix
        mxMatrix("Full", nrow = 6, ncol = 6,
                 name = "A",
                 values = valA,
                 free = freeA),
        # S matrix
        mxMatrix("Symm", nrow = 6, ncol = 6,
                 name = "S",
                 values = valS[lower.tri(valS, diag = TRUE)],
                 free = FALSE,
                 labels = labS
        ),
        # F matrix
        mxMatrix("Full", nrow = 4, ncol = 6, name = "F",
                 values = valF, free = FALSE),
        # M matrix (vector)
        mxMatrix("Full", nrow = 1, ncol = 6, name = "M",
                 free = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                 values = c(0, 0, 0.5, 1.5, 0, 0)),
        # SD
        mxAlgebra(cbind(1, 1, vf ^ 0.5), name = "Sv"),
        mxAlgebra(A[-1:-2,-1:-2], name = "subA"),
        mxAlgebra(vec2diag(Sv ^ -1) %*% subA %*% vec2diag(Sv),
                  name = "SA"),
        mxAlgebra(cbind(SA[4,1], SA[3,cbind(1,4)]), name = "beta"),
        mxAlgebra(beta[1,1] * beta[1,3], name = "std_ind"),
        mxExpectationRAM("A", "S", "F", "M",
            dimnames = c("fy_fs", "fm_fs", "treat", "group", "fy", "fm")
        ),
        mxFitFunctionML(),
        mxCI(c("beta", "std_ind"))
    )
    tspa_fit <- mxRun(tspa_mx, intervals = TRUE, silent = TRUE)
    if (!tspa_fit@output$status$status == 0) {
        stop("2SPA did not converge")
    }
    ExtractMx(tspa_fit)
}
# Test: 2S-PA (matrix)
# Analyse2spaMat(DESIGNFACTOR[1,], test_dat, fixed_objects = FIXED)

Analyse2spaMat2 <- function(condition, dat, fixed_objects = NULL) {
    freeA <- matrix(FALSE, nrow = 4, ncol = 4)
    valA <- 0 * freeA
    freeApos <- matrix(
        byrow = TRUE, ncol = 2,
        c(
            1, 2,
            1, 3,
            1, 4,
            2, 3,
            2, 4
        )
    )
    freeA[freeApos] <- TRUE
    valA[freeApos] <- c(condition$beta3, condition$beta2,
                        fixed_objects$dalphay,
                        condition$beta1, fixed_objects$dalpham)
    valS <- diag(c(condition$psiy, condition$psim, 0.25, 0.25))
    labS <- NA * valS
    diag(labS)[1:2] <- c("ev_fy[1,1]", "ev_fm[1,1]")
    labE <- matrix(NA, nrow = 4, ncol = 4)
    diag(labE)[1:2] <- c("ev_fs[1,1]", "ev_fs[1,2]")
    tspa_mx <- mxModel("2SPA",
        mxData(observed = dat, type = "raw"),
        manifestVars = c("fy_fs", "fm_fs", "treat", "group"),
        # Total variances
        mxMatrix(type = "Full", nrow = 1, ncol = 2,
                 values = c(mean(dat$fy_fs_rel), mean(dat$fm_fs_rel)),
                 free = TRUE, name = "v_fs"),
        # Definition variables
        mxMatrix(type = "Full", nrow = 1, ncol = 2, free = FALSE,
                 labels = c("data.fy_fs_rel", "data.fm_fs_rel"),
                 name = "fs_rel"),
        # True and error variances
        mxAlgebra(v_fs * (1 - fs_rel), name = "ev_fs"),
        mxAlgebra(v_fs * fs_rel, name = "vf"),
        # Sub-matrix of S
        mxMatrix("Symm", nrow = 2, ncol = 2,
                 name = "subS1",
                 values = valS[3:4, 3:4][lower.tri(valS[3:4, 3:4], diag = TRUE)],
                 free = FALSE,
                 labels = labS[3:4, 3:4]
        ),
        mxMatrix("Symm", nrow = 3, ncol = 3,
                 name = "subS2",
                 values = valS[2:4, 2:4][lower.tri(valS[2:4, 2:4], diag = TRUE)],
                 free = FALSE,
                 labels = labS[2:4, 2:4]
        ),
        mxAlgebra(vf[1,2] - A[2,3:4] %&% subS1, name = "ev_fm"),
        mxAlgebra(vf[1,1] - A[1,2:4] %&% subS2, name = "ev_fy"),
        # A matrix
        mxMatrix("Full", nrow = 4, ncol = 4,
                 name = "A",
                 values = valA,
                 free = freeA),
        # S matrix
        mxMatrix("Symm", nrow = 4, ncol = 4,
                 name = "S",
                 values = valS[lower.tri(valS, diag = TRUE)],
                 free = FALSE,
                 labels = labS
        ),
        # E matrix (measurement error variances)
        mxMatrix("Diag", nrow = 4, ncol = 4, name = "E",
                 values = diag(c(0.1, 0.1, 0, 0)),
                 free = FALSE,
                 labels = labE),
        # M matrix (vector)
        mxMatrix("Full", nrow = 1, ncol = 4, name = "M",
                 free = c(TRUE, TRUE, FALSE, FALSE),
                 values = c(0, 0, 0.5, 1.5)),
        # SD
        mxMatrix(type = "Iden", nrow = 4, name = "I"),
        mxAlgebra(cbind(vf ^ 0.5, 1, 1), name = "Sv"),
        mxAlgebra(vec2diag(Sv ^ -1) %*% A %*% vec2diag(Sv), name = "SA"),
        mxAlgebra(cbind(SA[2,3], SA[1,cbind(3,2)]), name = "beta"),
        mxAlgebra(beta[1,1] * beta[1,3], name = "std_ind"),
        mxAlgebra(solve(I - A) %&% S + E, name = "Sigma"),
        mxExpectationNormal("Sigma", means = "M",
            dimnames = c("fy_fs", "fm_fs", "treat", "group")
        ),
        mxFitFunctionML(),
        mxCI(c("beta", "std_ind"))
    )
    tspa_fit <- mxRun(tspa_mx, intervals = TRUE, silent = TRUE)
    if (!tspa_fit@output$status$status == 0) {
        stop("2SPA did not converge")
    }
    ExtractMx(tspa_fit)
}
# Test: 2S-PA (matrix with normal expectation)
# Analyse2spaMat2(DESIGNFACTOR[1,], test_dat, fixed_objects = FIXED)

# Evaluate ----------------------------------------------------------------

Evaluate <- function(condition, results, fixed_objects = NULL) {
    # ind_pop <- condition$beta1 * condition$beta3
    meth_name <- c("fspa", "2spa", "sem")
    pars <- c("beta1", "beta2", "beta3", "std_ind")
    pop <- c(
        condition$beta1,
        condition$beta2,
        condition$beta3,
        condition$beta1 * condition$beta3
    )
    # Repeat by methods
    pop_vals <- rep(pop, each = length(meth_name))
    meth_pars <- c(outer(meth_name, pars, FUN = paste, sep = "."))
    est <- results[, paste0(meth_pars, ".est")]
    ci_names <- c(t(outer(meth_pars, c("ll", "ul"), paste, sep = ".")))
    cis <- results[, ci_names]
    centered_cis <- sweep(cis,
        MARGIN = 2, STATS = rep(pop_vals, each = 2)
    )
    c(
        bias = bias(est, pop_vals),
        rmse = RMSE(est, pop_vals),
        coverage = ECR(centered_cis,
            parameter = 0
        ),
        power = 1 - ECR(cis,
            parameter = 0
        )
    )
}

res <- runSimulation(
    design = DESIGNFACTOR,
    replications = 2500,
    generate = AddFscores,
    analyse = list(
        sem = AnalyseSem,
        fspa = AnalyseFspaMat,
        `2spa` = Analyse2spaMat
    ),
    summarise = Evaluate,
    fixed_objects = FIXED,
    packages = c("OpenMx", "psych", "lavaan", "mirt"),
    filename = "simulation/mediate-results",
    seed = rep(21390926, nrow(DESIGNFACTOR)),
    parallel = TRUE,
    ncores = 10L,
    save = TRUE,
    max_errors = 100,
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "mediate-results_"
    )
)
