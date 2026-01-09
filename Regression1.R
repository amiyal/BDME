# Load required packages
library(Ecdat)
library(MASS)
library(COUNT)

# Load data
data(azdrg112)
df <- azdrg112
y <- df$los  

# Final formula with only significant covariates
final_formula <- los ~ type1 + gender + age75

# ---------------- Poisson Model ----------------
poisson_model <- glm(final_formula, data = df, family = poisson(link = "log"))
summary(poisson_model)

# ---------------- Dispersion Index ----------------
disp_index <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
cat("Dispersion Index:", round(disp_index, 2), "\n")

# ---------------- Negative Binomial Model ----------------
nb_model <- glm.nb(final_formula, data = df)

# ---------------- Geometric Model ----------------
ll_geom <- function(beta, data) {
  X <- model.matrix(final_formula, data = data)
  eta <- X %*% beta
  mu <- exp(eta)
  p <- 1 / (1 + mu)
  ll <- sum(log(p) + y * log(1 - p))  # ✅ y is used from global scope
  return(-ll)
}
init_beta <- coef(poisson_model)
fit_geom <- optim(init_beta, ll_geom, data = df, method = "BFGS", hessian = TRUE)
beta_geom <- fit_geom$par
log_like_geom <- -ll_geom(beta_geom, df)

# ---------------- BDME Model ----------------
ll_bdme <- function(beta, data) {
  X <- model.matrix(final_formula, data = data)
  eta <- X %*% beta
  mu <- exp(eta)
  A <- sqrt(mu^2 + 6 * mu + 1)
  ll <- sum(
    2 * log((3 + mu - A) / 2) +
      (2 * y - 1) * log((A - mu - 1) / 2) +
      log(y + ((A - mu - 1) * (3 + mu - A)) / 4) -
      (y + 2) * log(1 - ((A - mu - 1) * (3 + mu - A)) / 4)
  )
  return(-ll)
}
fit_bdme <- optim(init_beta, ll_bdme, data = df, method = "BFGS", hessian = TRUE)
beta_bdme <- fit_bdme$par
log_like_bdme <- -ll_bdme(beta_bdme, df)

# ---------------- Model Comparison ----------------
log_like_pois <- logLik(poisson_model)[1]
log_like_nb <- logLik(nb_model)[1]
n <- nrow(df)
k <- length(init_beta)
k_nb <- length(coef(nb_model)) + 1
loglog_n <- log(log(n))

model_compare <- data.frame(
  Model = c("Poisson", "Negative Binomial", "Geometric", "BDME"),
  LogLikelihood = round(c(log_like_pois, log_like_nb, log_like_geom, log_like_bdme), 3),
  AIC = round(c(
    -2 * log_like_pois + 2 * k,
    -2 * log_like_nb + 2 * k_nb,
    -2 * log_like_geom + 2 * k,
    -2 * log_like_bdme + 2 * k
  ), 2),
  BIC = round(c(
    -2 * log_like_pois + log(n) * k,
    -2 * log_like_nb + log(n) * k_nb,
    -2 * log_like_geom + log(n) * k,
    -2 * log_like_bdme + log(n) * k
  ), 2),
  HAIC = round(c(
    -2 * log_like_pois + 2 * k * loglog_n,
    -2 * log_like_nb + 2 * k_nb * loglog_n,
    -2 * log_like_geom + 2 * k * loglog_n,
    -2 * log_like_bdme + 2 * k * loglog_n
  ), 2)
)

cat("\n================ Model Comparison =================\n")
print(model_compare)

# ---------------- Coefficient Extraction ----------------
extract_results <- function(est, se) {
  z <- est / se
  pval <- 2 * (1 - pnorm(abs(z)))
  data.frame(
    Covariate = names(est),
    Estimate = round(est, 5),
    StdError = round(se, 5),
    z_value = round(z, 3),
    p_value = format.pval(pval, digits = 5, eps = .00001),
    row.names = NULL
  )
}

# Poisson
summary_pois <- summary(poisson_model)
pois_result <- extract_results(coef(poisson_model), summary_pois$coefficients[, "Std. Error"])

# Negative Binomial
summary_nb <- summary(nb_model)
nb_result <- extract_results(coef(nb_model), summary_nb$coefficients[, "Std. Error"])

# Geometric
vcov_geom <- MASS::ginv(fit_geom$hessian)  # Use ginv to avoid singularity error
se_geom <- sqrt(diag(vcov_geom))
geom_result <- extract_results(beta_geom, se_geom)

# BDME
vcov_bdme <- MASS::ginv(fit_bdme$hessian)
se_bdme <- sqrt(diag(vcov_bdme))
bdme_result <- extract_results(beta_bdme, se_bdme)

# ---------------- Print Results ----------------
cat("\n===== Poisson Regression Results =====\n")
print(pois_result, row.names = FALSE)

cat("\n===== Negative Binomial Regression Results =====\n")
print(nb_result, row.names = FALSE)

cat("\n===== Geometric Regression Results =====\n")
print(geom_result, row.names = FALSE)

cat("\n===== BDME Regression Results =====\n")
print(bdme_result, row.names = FALSE)
