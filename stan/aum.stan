
data {
  int<lower=1> n_person;
  int<lower=2> n_wave;
  int<lower=0> n_cov;                                       // initial-state predictors

  array[n_person, n_wave] int<lower=0, upper=1> observed;   // aspiration reported
  array[n_person, n_wave] real y;                           // aspiration, 0 where unobserved

  array[n_person, n_wave] real<lower=0> v_obs;

  array[n_wave] int<lower=0, upper=1> fielded_a;
  array[n_wave] int<lower=0, upper=1> fielded_sc;
  array[n_person, n_wave] int<lower=0, upper=1> obs_a;
  array[n_person, n_wave] int<lower=0, upper=1> obs_sc;
  array[n_person, n_wave] real s_a;                         // 0 where unobserved
  array[n_person, n_wave] real s_sc;

  array[n_wave] int<lower=0, upper=1> fielded_x;
  array[n_person, n_wave] int<lower=0, upper=1> obs_x;
  array[n_person, n_wave] real s_x;

  // The change in the mark since the previous fielded wave, and where both
  // marks were observed so that the change exists. Read only by model 5.
  array[n_person, n_wave] real news_a;
  array[n_person, n_wave] int<lower=0, upper=1> obs_news;

  matrix[n_person, n_cov] z;

  int<lower=0> n_mod;
  matrix[n_person, n_mod] w;

  int<lower=0, upper=7> model_id;
}

transformed data {
  // Model 3 is model 2 with one term added, so it takes model 2's blocks too.
  int n_beta  = (model_id == 1 || model_id >= 4) ? n_wave : 0;
  int n_phi   = (model_id >= 4) ? n_wave : 0;
  int n_news  = (model_id >= 5) ? 1 : 0;
  int n_bneg  = (model_id == 5 || model_id == 6) ? 1 : 0;
  int n_blvl  = (model_id == 6) ? 1 : 0;
  int n_kappa = (model_id == 2 || model_id == 3) ? n_wave : 0;
  int n_delta = (model_id == 2 || model_id == 3) ? n_mod : 0;
  int n_bx    = (model_id == 3) ? n_wave : 0;
  int n_cx    = (model_id == 3) ? n_cov : 0;
  int n_sx    = (model_id == 3) ? 1 : 0;
}

parameters {
  real gamma_0;
  vector[n_cov] gamma;

  real<lower=0.01> sigma_nu;                           // initial-state dispersion
  real<lower=0.01> sigma_eps;                          // measurement error
  real<lower=0.01> sigma_eta;                          // innovation

  // Auxiliary signal models. Identified by the observed signal values, and used
  // to integrate out the unobserved ones.
  vector[n_wave] mu_a;
  vector[n_wave] mu_sc;
  vector[n_cov] alpha_a;
  vector[n_cov] alpha_sc;
  real<lower=0.01> sigma_a;
  real<lower=0.01> sigma_sc;

  // M1: unconstrained response to signal levels
  vector[n_beta] beta_a;
  vector[n_beta] beta_sc;

  // M4: the autoregressive coefficient, free rather than fixed at one or tied
  // to the signal coefficients.
  vector[n_phi] phi;

  // M5: the response to good and to bad news about the mark. Under M7 the two
  // are tied, so b_neg is not allocated and b_pos carries both directions.
  vector[n_news] b_pos;
  vector[n_bneg] b_neg;

  vector[n_blvl] b_lvl;

  // M2: total permeability and its split between the two signals, on the logit
  // scale so the (0,1) and sum-to-under-one constraints hold by construction.
  vector[n_kappa] kappa_logit;
  vector[n_kappa] pi_logit;
  vector[n_delta] delta;

  // M3: the cross-domain push, plus its own auxiliary signal model so that an
  // unobserved cross-domain mark integrates out exactly as the others do.
  vector[n_bx] beta_x;
  vector[n_bx] mu_x;
  vector[n_cx] alpha_x;
  vector<lower=0.01>[n_sx] sigma_x;
}

model {
  gamma_0 ~ normal(0, 2);
  gamma ~ normal(0, 1);
  sigma_nu ~ normal(0, 1);
  sigma_eps ~ normal(0, 1);
  sigma_eta ~ normal(0, 1);

  mu_a ~ normal(0, 1);
  mu_sc ~ normal(0, 1);
  alpha_a ~ normal(0, 1);
  alpha_sc ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_sc ~ normal(0, 1);

  beta_a ~ normal(0, 0.5);
  beta_sc ~ normal(0, 0.5);
  // Centred on one, which is M1's assumption, so the prior does not lean
  // towards the decay that M2 implies.
  phi ~ normal(1, 0.5);
  b_pos ~ normal(0, 0.5);
  b_neg ~ normal(0, 0.5);
  b_lvl ~ normal(0, 0.5);
  kappa_logit ~ normal(-1, 1);                         // centred near a 27% gap closure
  pi_logit ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  beta_x ~ normal(0, 0.5);
  mu_x ~ normal(0, 1);
  alpha_x ~ normal(0, 1);
  sigma_x ~ normal(0, 1);

  for (i in 1:n_person) {
    real zg = (n_cov > 0) ? dot_product(z[i], gamma) : 0;
    real za = (n_cov > 0) ? dot_product(z[i], alpha_a) : 0;
    real zs = (n_cov > 0) ? dot_product(z[i], alpha_sc) : 0;
    real zx = (n_cx > 0) ? dot_product(z[i], alpha_x) : 0;

    real m = gamma_0 + zg;                             // Kalman prior at wave 1
    real p = square(sigma_nu);

    for (t in 1:n_wave) {
      // Auxiliary signal likelihood, which is what identifies mu, alpha, sigma.
      // Excluded from log_lik below: model comparison scores aspirations only.
      if (obs_a[i, t] == 1)  s_a[i, t] ~ normal(mu_a[t] + za, sigma_a);
      if (obs_sc[i, t] == 1) s_sc[i, t] ~ normal(mu_sc[t] + zs, sigma_sc);
      if (n_bx > 0 && obs_x[i, t] == 1)
        s_x[i, t] ~ normal(mu_x[t] + zx, sigma_x[1]);

      if (t > 1) {
        real c = 1;
        real a = 0;
        real q = 0;

        // Value entering the transition: the observation if there is one,
        // otherwise its conditional mean, whose variance is then carried in q.
        real val_a = (obs_a[i, t] == 1) ? s_a[i, t] : mu_a[t] + za;
        real val_sc = (obs_sc[i, t] == 1) ? s_sc[i, t] : mu_sc[t] + zs;

        if (model_id == 1 || model_id >= 4) {
          if (model_id >= 4) c = phi[t];
          if (fielded_a[t] == 1) {
            a += beta_a[t] * val_a;
            if (obs_a[i, t] == 0) q += square(beta_a[t] * sigma_a);
            if (model_id >= 5 && obs_news[i, t] == 1) {
              a += ((model_id == 7 || news_a[i, t] > 0) ? b_pos[1] : b_neg[1])
                   * news_a[i, t];
              if (model_id == 6) {
                a += b_lvl[1] * news_a[i, t] * (val_a - news_a[i, t]);
              }
            }
          }
          if (fielded_sc[t] == 1) {
            a += beta_sc[t] * val_sc;
            if (obs_sc[i, t] == 0) q += square(beta_sc[t] * sigma_sc);
          }
        } else if (model_id == 2 || model_id == 3) {
          real kappa = inv_logit(kappa_logit[t] +
                                 ((n_mod > 0) ? dot_product(w[i], delta) : 0));
          real lam_a = 0;
          real lam_sc = 0;

          if (fielded_a[t] == 1 && fielded_sc[t] == 1) {
            real share = inv_logit(pi_logit[t]);
            lam_a = kappa * share;
            lam_sc = kappa * (1 - share);
          } else if (fielded_a[t] == 1) {
            lam_a = kappa;
          } else if (fielded_sc[t] == 1) {
            lam_sc = kappa;
          }

          c = 1 - lam_a - lam_sc;
          a = lam_a * val_a + lam_sc * val_sc;

          if (fielded_a[t] == 1 && obs_a[i, t] == 0)   q += square(lam_a * sigma_a);
          if (fielded_sc[t] == 1 && obs_sc[i, t] == 0) q += square(lam_sc * sigma_sc);

          // Added after the assignments above, which use `=`, not `+=`.
          if (model_id == 3 && fielded_x[t] == 1) {
            real val_x = (obs_x[i, t] == 1) ? s_x[i, t] : mu_x[t] + zx;
            a += beta_x[t] * val_x;
            if (obs_x[i, t] == 0) q += square(beta_x[t] * sigma_x[1]);
          }
        }

        m = c * m + a;
        p = square(c) * p + square(sigma_eta) + q;
      }

      if (observed[i, t] == 1) {
        real v = y[i, t] - m;                          // innovation
        real f = p + square(sigma_eps) + v_obs[i, t];  // innovation variance
        target += normal_lpdf(v | 0, sqrt(f));
        m += p * v / f;
        p -= p * p / f;
      }
    }
  }
}

generated quantities {
  vector[n_person] log_lik;

  real bayes_r2;
  real s_m = 0;
  real s_mm = 0;
  real s_vv = 0;
  int n_obs = 0;

  // Updating rates on the natural scale at the reference category, so the
  // developmental trajectory can be read straight off the draws.
  vector[n_kappa] lambda_a;
  vector[n_kappa] lambda_sc;

  for (t in 1:n_kappa) {
    real kappa = inv_logit(kappa_logit[t]);
    real share = inv_logit(pi_logit[t]);
    // No transition enters the first wave, so no rate is identified there and
    // kappa_logit[1] merely samples its prior. Do not report it.
    if (t == 1) {
      lambda_a[t] = 0;
      lambda_sc[t] = 0;
    } else if (fielded_a[t] == 1 && fielded_sc[t] == 1) {
      lambda_a[t] = kappa * share;
      lambda_sc[t] = kappa * (1 - share);
    } else if (fielded_a[t] == 1) {
      lambda_a[t] = kappa;
      lambda_sc[t] = 0;
    } else if (fielded_sc[t] == 1) {
      lambda_a[t] = 0;
      lambda_sc[t] = kappa;
    } else {
      lambda_a[t] = 0;
      lambda_sc[t] = 0;
    }
  }

  for (i in 1:n_person) {
    real zg = (n_cov > 0) ? dot_product(z[i], gamma) : 0;
    real za = (n_cov > 0) ? dot_product(z[i], alpha_a) : 0;
    real zs = (n_cov > 0) ? dot_product(z[i], alpha_sc) : 0;
    real zx = (n_cx > 0) ? dot_product(z[i], alpha_x) : 0;

    real m = gamma_0 + zg;
    real p = square(sigma_nu);
    log_lik[i] = 0;

    for (t in 1:n_wave) {
      if (t > 1) {
        real c = 1;
        real a = 0;
        real q = 0;
        real val_a = (obs_a[i, t] == 1) ? s_a[i, t] : mu_a[t] + za;
        real val_sc = (obs_sc[i, t] == 1) ? s_sc[i, t] : mu_sc[t] + zs;

        if (model_id == 1 || model_id >= 4) {
          if (model_id >= 4) c = phi[t];
          if (fielded_a[t] == 1) {
            a += beta_a[t] * val_a;
            if (obs_a[i, t] == 0) q += square(beta_a[t] * sigma_a);
            if (model_id >= 5 && obs_news[i, t] == 1) {
              a += ((model_id == 7 || news_a[i, t] > 0) ? b_pos[1] : b_neg[1])
                   * news_a[i, t];
              if (model_id == 6) {
                a += b_lvl[1] * news_a[i, t] * (val_a - news_a[i, t]);
              }
            }
          }
          if (fielded_sc[t] == 1) {
            a += beta_sc[t] * val_sc;
            if (obs_sc[i, t] == 0) q += square(beta_sc[t] * sigma_sc);
          }
        } else if (model_id == 2 || model_id == 3) {
          real kappa = inv_logit(kappa_logit[t] +
                                 ((n_mod > 0) ? dot_product(w[i], delta) : 0));
          real lam_a = 0;
          real lam_sc = 0;

          if (fielded_a[t] == 1 && fielded_sc[t] == 1) {
            real share = inv_logit(pi_logit[t]);
            lam_a = kappa * share;
            lam_sc = kappa * (1 - share);
          } else if (fielded_a[t] == 1) {
            lam_a = kappa;
          } else if (fielded_sc[t] == 1) {
            lam_sc = kappa;
          }

          c = 1 - lam_a - lam_sc;
          a = lam_a * val_a + lam_sc * val_sc;

          if (fielded_a[t] == 1 && obs_a[i, t] == 0)   q += square(lam_a * sigma_a);
          if (fielded_sc[t] == 1 && obs_sc[i, t] == 0) q += square(lam_sc * sigma_sc);

          // Added after the assignments above, which use `=`, not `+=`.
          if (model_id == 3 && fielded_x[t] == 1) {
            real val_x = (obs_x[i, t] == 1) ? s_x[i, t] : mu_x[t] + zx;
            a += beta_x[t] * val_x;
            if (obs_x[i, t] == 0) q += square(beta_x[t] * sigma_x[1]);
          }
        }

        m = c * m + a;
        p = square(c) * p + square(sigma_eta) + q;
      }

      if (observed[i, t] == 1) {
        real v = y[i, t] - m;
        real f = p + square(sigma_eps) + v_obs[i, t];
        log_lik[i] += normal_lpdf(v | 0, sqrt(f));
        s_m += m;
        s_mm += square(m);
        s_vv += square(v);
        n_obs += 1;
        m += p * v / f;
        p -= p * p / f;
      }
    }
  }

  {
    real var_fit = s_mm / n_obs - square(s_m / n_obs);
    real var_res = s_vv / n_obs;
    bayes_r2 = var_fit / (var_fit + var_res);
  }
}

