
######## Tidy script for Empirical Example ########

#### Load libraries ####
library(tidyverse)
library(lavaan)
library(pinsearch)
library(kableExtra)
library(OpenMx)


#### Load data ####
orig_dat <- haven::read_sav(here::here("replication", "models_dat", 
                                       "LUI2018_class4mplus.sav"))
dat <- orig_dat %>%
  filter_at(vars(class14, audit1:audit3), ~ !is.na(.)) %>% 
  filter(eth %in% 1:4) %>%
  mutate_at(vars(class1:class15, audit1:audit10), as.numeric) %>%
  mutate(campuslive = factor(campuslive), 
         eth = factor(eth, 1:4)) %>%
  arrange(eth)


#### CLASS Invariance Testing ####

# baseline model
c_mod <- "
f1 =~ class1 + class2 + class3 + class4 + class5 + class6 + 
      class7 + class8 + class9 + class10 + class11 + class12 + 
      class13 + class14 + class15
"

# configural invariance
c_confit <- cfa(model = c_mod, 
                data = dat, 
                group = "eth", 
                group.label = 1:4, 
                std.lv = TRUE)
# metric invariance
c_metfit <- cfa(model = c_mod, 
                data = dat, 
                group = "eth", 
                group.label = 1:4, 
                group.equal = c("loadings"), 
                std.lv = TRUE)
# scalar invariance
c_scalfit <- cfa(model = c_mod, 
                 data = dat, 
                 group = "eth", 
                 group.label = 1:4, 
                 group.equal = c("loadings", "intercepts"), 
                 std.lv = TRUE)

# search for noninvariant intercepts 
# pinSearch(c_mod, 
#           data = filter(dat, eth %in% 1:4), 
#           group = "eth", 
#           group.label = 1:4, 
#           type = "intercepts")

# partial scalar model
c_ps_mod <- "
f1 =~ class1 + class2 + class3 + class4 + class5 + class6 + 
      class7 + class8 + class9 + class10 + class11 + class12 + 
      class13 + class14 + class15

class1 ~ c(i1, i1, i1, i1.h)*1
class2 ~ c(i2, i2, i2, i2.h)*1
class4 ~ c(i4.w, i4, i4, i4)*1
class6 ~ c(i6, i6, i6.b, i6)*1
class10 ~ c(i10, i10.a, i10, i10)*1
class11 ~ c(i11, i11, i11, i11.h)*1
class12 ~ c(i12, i12.a, i12, i12)*1
class13 ~ c(i13, i13, i13.b, i13)*1
class14 ~ c(i14.w, i14, i14, i14)*1
class15 ~ c(i15, i15.a, i15, i15)*1
"
# partial scalar invariance
c_psfit <- cfa(model = c_ps_mod, 
               data = dat, 
               group = "eth", 
               group.label = 1:4, 
               group.equal = c("loadings", "intercepts"), 
               std.lv = TRUE)

#### AUDIT Invariance Testing ####

# baseline model
a_mod <- "
f1 =~ audit4 + audit5 + audit6 + audit7 + audit8 + audit9 + audit10
"

# configural invariance
a_confit <- cfa(model = a_mod, 
                data = dat, 
                group = "eth", 
                group.label = 1:4, 
                std.lv = TRUE)
# metric invariance
a_metfit <- cfa(model = a_mod, 
                data = dat, 
                group = "eth", 
                group.label = 1:4, 
                group.equal = c("loadings"), 
                std.lv = TRUE)

# search for noninvariant loadings 
# pinSearch(a_mod, 
#           data = dat, 
#           group = "eth", 
#           group.label = 1:4, 
#           type = "loadings")

# partial metric model
a_pm_mod <- "
f1 =~ c(l4, l4, l4, l4.h)*audit4 + audit5 + c(l6, l6.a, l6.b, l6)*audit6 + audit7 + 
      c(l8.w, l8, l8, l8)*audit8 + audit9 + c(l10, l10, l10, l10.h)*audit10
"
# partial metric invariance
a_pmfit <- cfa(model = a_pm_mod, 
               data = dat, 
               group = "eth", 
               group.label = 1:4, 
               group.equal = c("loadings"), 
               std.lv = TRUE)

# scalar invariance
a_scalfit <- cfa(model = a_pm_mod, 
                 data = dat, 
                 group = "eth", 
                 group.label = 1:4, 
                 group.equal = c("loadings", "intercepts"), 
                 std.lv = TRUE)

# search for noninvariant intercepts 
# pinSearch(a_pm_mod, 
#           data = dat, 
#           group = "eth", 
#           group.label = 1:4, 
#           type = "intercepts")

# partial scalar model
a_ps_mod <- "
f1 =~ c(l4, l4, l4, l4.h)*audit4 + audit5 + c(l6, l6.a, l6.b, l6)*audit6 + audit7 + 
      c(l8.w, l8, l8, l8)*audit8 + audit9 + c(l10, l10, l10, l10.h)*audit10
      
audit5 ~ c(i5, i5, i5, i5.h)*1
audit6 ~ c(i6, i6, i6.b, i6)*1
audit8 ~ c(i8.w, i8, i8, i8)*1
audit9 ~ c(i9.w, i9, i9, i9)*1
"
# partial scalar invariance
a_psfit <- cfa(model = a_ps_mod, 
               data = dat, 
               group = "eth", 
               group.label = 1:4, 
               group.equal = c("loadings", "intercepts"), 
               std.lv = TRUE)


fit <- lapply(list(c_confit, c_metfit, c_psfit, c_scalfit, 
                   a_confit, a_pmfit, a_metfit, a_psfit, a_scalfit), 
              function(x) { fitmeasures(x)[c("chisq", "df", "rmsea", 
                                             "cfi", "tli", "srmr")] }) %>%
  sapply(c)
lrt <- t(rbind(NA, lavTestLRT(c_confit, c_metfit, c_psfit)[2:3, 5:7], 
               lavTestLRT(c_metfit, c_scalfit)[2, 5:7], 
               NA, lavTestLRT(a_confit, a_pmfit)[2, 5:7], 
               lavTestLRT(a_confit, a_metfit)[2, 5:7], 
               lavTestLRT(a_pmfit, a_psfit)[2, 5:7], 
               lavTestLRT(a_pmfit, a_scalfit)[2, 5:7]))
delta_fit <- cbind(NA, t(diff(t(fit[c("cfi", "rmsea"), 1:4]))), 
                   NA, t(diff(t(fit[c("cfi", "rmsea"), 5:9]))))

#### Summary Table ####
fit_tab <- rbind(fit, lrt, delta_fit)
rownames(fit_tab) <- c("chisq", "df", "rmsea", "cfi", "tli", "srmr", 
                       "delta_chisq", "delta_df", "p", "delta_cfi", "delta_rmsea")
colnames(fit_tab) <- c(paste0("c_", c("conf", "met", "ps", "scal")), 
                       paste0("a_", c("conf", "pm", "met", "ps", "scal")))


#### Extract Factor Scores, SEs, and Reliability ####

# CLASS (X, predictor)
fx_fs <- lavPredict(c_psfit, se = "standard")
fx_se <- rep(unlist(attributes(fx_fs)$se), 
             unlist(lapply(fx_fs, length)))
fx_psi <- rep(subset(parameterestimates(c_psfit), 
                     lhs == "f1" & op == "~~")$est, 
              unlist(lapply(fx_fs, length)))
fx_rel <- 1 - fx_se^2 / fx_psi

# AUDIT (Y, outcome)
fy_fs <- lavPredict(a_psfit, se = "standard")
fy_se <- rep(unlist(attributes(fy_fs)$se), 
             unlist(lapply(fy_fs, length)))
fy_psi <- rep(subset(parameterestimates(a_psfit), 
                     lhs == "f1" & op == "~~")$est, 
              unlist(lapply(fx_fs, length)))
fy_rel <- 1 - fy_se^2 / fy_psi

# Append to dat
dat <- dat %>%
  arrange(eth) %>%
  mutate(class_sum = rowSums(across(starts_with("class"))),
         audit_sum = rowSums(across(starts_with("audit"))),         eth2 = ifelse(eth == 2, 1, 0),
         eth3 = ifelse(eth == 3, 1, 0),
         eth4 = ifelse(eth == 4, 1, 0)) %>%
  add_column(fx_fs = unlist(fx_fs),
             fx_fs_se = fx_se,
             fx_fs_rel = fx_rel,
             fy_fs = unlist(fy_fs),
             fy_fs_se = fy_se,
             fy_fs_rel = fy_rel) %>%
  as.data.frame()


# #### Linear Regression with Composite Scores ####
# 
# lm_ac <- lm(scale(audit_sum) ~ scale(class_sum) + eth, data = dat)
# lm_b <- as.numeric(coef(lm_ac)[2])
# lm_b_ci <- as.numeric(confint(lm_ac)[2, ])
# 
# # standardized regression coefficient
# lm_b <- as.numeric(coef(lm_ac)[2])
# # confidence interval
# lm_b_ci <- as.numeric(confint(lm_ac)[2, ])

#### Linear Regression with lavaan ####

lm_mod <- "
audit_sum ~ b * class_sum + eth2 + eth3 + eth4
"

lmfit <- sem(model = lm_mod, 
             data = dat)

# standardized path coefficient
lm_b <- (standardizedSolution(lmfit) %>%
           filter(label == "b"))$est
# confidence interval
lm_b_ci <- (standardizedSolution(lmfit) %>%
              filter(label == "b"))[, c("ci.lower", "ci.upper")] %>%
  as.numeric()


#### FS-PA ####

fspa_mod <- "
fy_fs ~ b * fx_fs + eth2 + eth3 + eth4
"

fspafit <- sem(model = fspa_mod, 
               data = dat)

# standardized path coefficient
fspa_b <- (standardizedSolution(fspafit) %>%
           filter(label == "b"))$est
# confidence interval
fspa_b_ci <- (standardizedSolution(fspafit) %>%
              filter(label == "b"))[, c("ci.lower", "ci.upper")] %>%
  as.numeric()

#### Full SEM ####

n_g <- table(dat$eth)

# full SEM (poor fit)
sem_mod <- "
# measurement model of X

fx =~ class1 + class2 + class3 + class4 + class5 + class6 + 
      class7 + class8 + class9 + class10 + class11 + class12 + 
      class13 + class14 + class15

class1 ~ c(i1, i1, i1, i1.h)*1
class2 ~ c(i2, i2, i2, i2.h)*1
class4 ~ c(i4.w, i4, i4, i4)*1
class6 ~ c(i6, i6, i6.b, i6)*1
class10 ~ c(i10, i10.a, i10, i10)*1
class11 ~ c(i11, i11, i11, i11.h)*1
class12 ~ c(i12, i12.a, i12, i12)*1
class13 ~ c(i13, i13, i13.b, i13)*1
class14 ~ c(i14.w, i14, i14, i14)*1
class15 ~ c(i15, i15.a, i15, i15)*1

# measurement model of Y

fy =~ c(l4, l4, l4, l4.h)*audit4 + audit5 + c(l6, l6.a, l6.b, l6)*audit6 + audit7 + 
      c(l8.w, l8, l8, l8)*audit8 + audit9 + c(l10, l10, l10, l10.h)*audit10

audit5 ~ c(i5, i5, i5, i5.h)*1
audit6 ~ c(i6, i6, i6.b, i6)*1
audit8 ~ c(i8.w, i8, i8, i8)*1
audit9 ~ c(i9.w, i9, i9, i9)*1

# factor means
fx ~ c(0*mx.w, mx.a, mx.b, mx.h)*1
fy ~ c(0*my.w, my.a, my.b, my.h)*1

# factor variances
fx ~~ c(1*vx.w, vx.a, vx.b, vx.h)*fx
fy ~~ c(1*evy.w, evy.a, evy.b, evy.h)*fy

# structural model
fy ~ c(b, b, b, b)*fx

# grand factor means
gm_x := (392 * mx.w + 186 * mx.a + 101 * mx.b + 163 * mx.h) / 842
gm_y := (392 * my.w + 186 * my.a + 101 * my.b + 163 * my.h) / 842

# grand factor variances
gv_x := (392 * (vx.w + (mx.w - gm_x)^2) + 186 * (vx.a + (mx.a - gm_x)^2) +
          101 * (vx.b + (mx.b - gm_x)^2) + 163 * (vx.h + (mx.h - gm_x)^2)) / 842
gv_y := (392 * ((b^2 * vx.w + evy.w) + (my.w - gm_y)^2) +
          186 * ((b^2 * vx.a + evy.a) + (my.a - gm_y)^2) +
          101 * ((b^2 * vx.b + evy.b) + (my.b - gm_y)^2) +
          163 * ((b^2 * vx.h + evy.h) + (my.h - gm_y)^2)) / 842

# standardized beta
beta := b * sqrt(gv_x / gv_y)
"

# run SEM
semfit <- sem(model = sem_mod, 
              data = dat, 
              group = "eth",
              group.label = 1:4,
              group.equal = c("loadings", "intercepts", "regressions"))

# standardized path coefficient
sem_b <- (parameterestimates(semfit) %>%
            filter(label == "beta"))$est
# confidence interval
sem_b_ci <- (parameterestimates(semfit, ci = TRUE) %>%
               filter(label == "beta"))[, c("ci.lower", "ci.upper")] %>%
  as.numeric()


#### 2S-PA ####

sem1_mx <- mxModel("SEM1",
                   type = "RAM",
                   mxData(observed = dat[dat$eth == 1, ], type = "raw"),
                   manifestVars = c("fy_fs", "fx_fs"),
                   latentVars = c("fy", "fx"),
                   # Total variances
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fy_fs1", free = TRUE),
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fx_fs1", free = TRUE),
                   mxAlgebra(v_fy_fs1 * (1 - data.fy_fs_rel), name = "ev_fy_fs1"),
                   mxAlgebra(v_fx_fs1 * (1 - data.fx_fs_rel), name = "ev_fx_fs1"),
                   mxAlgebra(v_fy_fs1 * data.fy_fs_rel, name = "v_fy1"),
                   mxAlgebra(v_fx_fs1 * data.fx_fs_rel, name = "v_fx1"),
                   mxAlgebra(v_fy1 - v_fx1 * b^2, name = "ev_fy1"), 
                   # mxAlgebra(b * sqrt(v_fx1 / v_fy1), name = "beta1"),
                   # Loadings
                   mxPath(from = c("fy", "fx"), to = c("fy_fs", "fx_fs"),
                          values = c(1, 1), free = FALSE),
                   # Path
                   mxPath(from = "fx", to = "fy",
                          free = TRUE, labels = "b"),
                   # Variances
                   mxPath(from = c("fy", "fx"), arrows = 2,
                          free = FALSE, values = c(0.5, 0.5),
                          labels = c("ev_fy1[1,1]", "v_fx1[1,1]")),
                   # Unique variances
                   mxPath(from = c("fy_fs", "fx_fs"), arrows = 2,
                          free = c(FALSE, FALSE), values = c(0.2, 0.2),
                          labels = c("ev_fy_fs1[1,1]", "ev_fx_fs1[1,1]")),
                   # Mean
                   mxPath(from = "one", to = c("fy", "fx"),
                          free = c(TRUE, TRUE), values = c(0, 0), 
                          labels = c("fy_m1", "fx_m1")),
                   # Intercepts
                   mxPath(from = "one", to = c("fy_fs", "fx_fs"),
                          free = FALSE, values = c(0, 0)), 
                   mxCI(c("b"))
)
sem2_mx <- mxModel("SEM2",
                   type = "RAM",
                   mxData(observed = dat[dat$eth == 2, ], type = "raw"),
                   manifestVars = c("fy_fs", "fx_fs"),
                   latentVars = c("fy", "fx"),
                   # Total variances
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fy_fs2", free = TRUE),
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fx_fs2", free = TRUE),
                   mxAlgebra(v_fy_fs2 * (1 - data.fy_fs_rel), name = "ev_fy_fs2"),
                   mxAlgebra(v_fx_fs2 * (1 - data.fx_fs_rel), name = "ev_fx_fs2"),
                   mxAlgebra(v_fy_fs2 * data.fy_fs_rel, name = "v_fy2"),
                   mxAlgebra(v_fx_fs2 * data.fx_fs_rel, name = "v_fx2"),
                   mxAlgebra(v_fy2 - v_fx2 * b^2, name = "ev_fy2"), 
                   # mxAlgebra(b * sqrt(v_fx2 / v_fy2), name = "beta2"),
                   # Loadings
                   mxPath(from = c("fy", "fx"), to = c("fy_fs", "fx_fs"),
                          values = c(1, 1), free = FALSE),
                   # Path
                   mxPath(from = "fx", to = "fy",
                          free = TRUE, labels = "b"),
                   # Variances
                   mxPath(from = c("fy", "fx"), arrows = 2,
                          free = FALSE, values = c(0.5, 0.5),
                          labels = c("ev_fy2[1,1]", "v_fx2[1,1]")),
                   # Unique variances
                   mxPath(from = c("fy_fs", "fx_fs"), arrows = 2,
                          free = c(FALSE, FALSE), values = c(0.2, 0.2),
                          labels = c("ev_fy_fs2[1,1]", "ev_fx_fs2[1,1]")),
                   # Mean
                   mxPath(from = "one", to = c("fy", "fx"),
                          free = c(TRUE, TRUE), values = c(0, 0), 
                          labels = c("fy_m2", "fx_m2")),
                   # Intercepts
                   mxPath(from = "one", to = c("fy_fs", "fx_fs"),
                          free = FALSE, values = c(0, 0)), 
                   mxCI(c("b"))
)
sem3_mx <- mxModel("SEM3",
                   type = "RAM",
                   mxData(observed = dat[dat$eth == 3, ], type = "raw"),
                   manifestVars = c("fy_fs", "fx_fs"),
                   latentVars = c("fy", "fx"),
                   # Total variances
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fy_fs3", free = TRUE),
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fx_fs3", free = TRUE),
                   mxAlgebra(v_fy_fs3 * (1 - data.fy_fs_rel), name = "ev_fy_fs3"),
                   mxAlgebra(v_fx_fs3 * (1 - data.fx_fs_rel), name = "ev_fx_fs3"),
                   mxAlgebra(v_fy_fs3 * data.fy_fs_rel, name = "v_fy3"),
                   mxAlgebra(v_fx_fs3 * data.fx_fs_rel, name = "v_fx3"),
                   mxAlgebra(v_fy3 - v_fx3 * b^2, name = "ev_fy3"), 
                   # mxAlgebra(b * sqrt(v_fx3 / v_fy3), name = "beta3"),
                   # Loadings
                   mxPath(from = c("fy", "fx"), to = c("fy_fs", "fx_fs"),
                          values = c(1, 1), free = FALSE),
                   # Path
                   mxPath(from = "fx", to = "fy",
                          free = TRUE, labels = "b"),
                   # Variances
                   mxPath(from = c("fy", "fx"), arrows = 2,
                          free = FALSE, values = c(0.5, 0.5),
                          labels = c("ev_fy3[1,1]", "v_fx3[1,1]")),
                   # Unique variances
                   mxPath(from = c("fy_fs", "fx_fs"), arrows = 2,
                          free = c(FALSE, FALSE), values = c(0.2, 0.2),
                          labels = c("ev_fy_fs3[1,1]", "ev_fx_fs3[1,1]")),
                   # Mean
                   mxPath(from = "one", to = c("fy", "fx"),
                          free = c(TRUE, TRUE), values = c(0, 0), 
                          labels = c("fy_m3", "fx_m3")),
                   # Intercepts
                   mxPath(from = "one", to = c("fy_fs", "fx_fs"),
                          free = FALSE, values = c(0, 0)), 
                   mxCI(c("b"))
)
sem4_mx <- mxModel("SEM4",
                   type = "RAM",
                   mxData(observed = dat[dat$eth == 4, ], type = "raw"),
                   manifestVars = c("fy_fs", "fx_fs"),
                   latentVars = c("fy", "fx"),
                   # Total variances
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fy_fs4", free = TRUE),
                   mxMatrix("Full", nrow = 1, ncol = 1,
                            name = "v_fx_fs4", free = TRUE),
                   mxAlgebra(v_fy_fs4 * (1 - data.fy_fs_rel), name = "ev_fy_fs4"),
                   mxAlgebra(v_fx_fs4 * (1 - data.fx_fs_rel), name = "ev_fx_fs4"),
                   mxAlgebra(v_fy_fs4 * data.fy_fs_rel, name = "v_fy4"),
                   mxAlgebra(v_fx_fs4 * data.fx_fs_rel, name = "v_fx4"),
                   mxAlgebra(v_fy4 - v_fx4 * b^2, name = "ev_fy4"), 
                   # mxAlgebra(b * sqrt(v_fx4 / v_fy4), name = "beta4"),
                   # Loadings
                   mxPath(from = c("fy", "fx"), to = c("fy_fs", "fx_fs"),
                          values = c(1, 1), free = FALSE),
                   # Path
                   mxPath(from = "fx", to = "fy",
                          free = TRUE, labels = "b"),
                   # Variances
                   mxPath(from = c("fy", "fx"), arrows = 2,
                          free = FALSE, values = c(0.5, 0.5),
                          labels = c("ev_fy4[1,1]", "v_fx4[1,1]")),
                   # Unique variances
                   mxPath(from = c("fy_fs", "fx_fs"), arrows = 2,
                          free = c(FALSE, FALSE), values = c(0.2, 0.2),
                          labels = c("ev_fy_fs4[1,1]", "ev_fx_fs4[1,1]")),
                   # Mean
                   mxPath(from = "one", to = c("fy", "fx"),
                          free = c(TRUE, TRUE), values = c(0, 0), 
                          labels = c("fy_m4", "fx_m4")),
                   # Intercepts
                   mxPath(from = "one", to = c("fy_fs", "fx_fs"),
                          free = FALSE, values = c(0, 0)), 
                   mxCI(c("b"))
)
mgsem_mx <- mxModel(
  "MultipleGroupCFA",
  sem1_mx, sem2_mx, sem3_mx, sem4_mx, 
  mxAlgebra((392 * fx_m1 + 186 * fx_m2 + 101 * fx_m3 + 163 * fx_m4) / 842, 
            name = "gm_x"), 
  mxAlgebra((392 * fy_m1 + 186 * fy_m2 + 101 * fy_m3 + 163 * fy_m4) / 842, 
            name = "gm_y"),
  mxAlgebra((392 * (SEM1.v_fx1 + (fx_m1 - gm_x)^2) + 186 * (SEM2.v_fx2 + (fx_m2 - gm_x)^2) +
               101 * (SEM3.v_fx3 + (fx_m3 - gm_x)^2) + 163 * (SEM4.v_fx4 + (fx_m4 - gm_x)^2)) / 842,
            name = "gv_x"),
  mxAlgebra((392 * ((b^2 * SEM1.v_fx1 + SEM1.ev_fy1) + (fy_m1 - gm_y)^2) +
               186 * ((b^2 * SEM2.v_fx2 + SEM2.ev_fy2) + (fy_m2 - gm_y)^2) +
               101 * ((b^2 * SEM3.v_fx3 + SEM3.ev_fy3) + (fy_m3 - gm_y)^2) +
               163 * ((b^2 * SEM4.v_fx4 + SEM4.ev_fy4) + (fy_m4 - gm_y)^2)) / 842,
            name = "gv_y"),
  mxAlgebra(b * sqrt(gv_x / gv_y), name = "beta"),
  mxCI("beta"),
  mxFitFunctionMultigroup(c("SEM1", "SEM2", "SEM3", "SEM4"))
)
tspa_g_fit <- mxRun(mgsem_mx, intervals = TRUE, silent = TRUE)
(tspa_sum <- summary(tspa_g_fit))

tspa_b <- tspa_sum$CI[1, "estimate"]
tspa_b_ci <- as.numeric(tspa_sum$CI[1, c("lbound", "ubound")])


#### Comparison ####
lm_b
fspa_b
sem_b
tspa_b

lm_b_ci
fspa_b_ci
sem_b_ci
tspa_b_ci



