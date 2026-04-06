# PII

This repository contains MATLAB code for computing the **Practical Identifiability Index (PII)** for a range of epidemic and growth models under different observation-error structures and calibration-window lengths.

## Repository structure

- `Final_code_PII/`  
  Main MATLAB code for simulation, estimation, bootstrap confidence intervals, and plotting.

- `options/`  
  Option files specifying the model, estimated parameters, true values, bounds, and other configuration details.

## How to run the code

To run an experiment, first make sure that the required **options file** is available in the working folder used by the main code. Then edit the user settings in the main script according to the experiment you want to run.

The main user settings are:

```matlab
% ============================================================================
% ============================================================================
% USER SETTINGS — Change these as needed
% ============================================================================
% ============================================================================

% Options file handle (must be located in the current folder).
%   @options_forecast_PII_EXPO_r_dist1_1        (EXP, Poisson)
%   @options_forecast_PII_EXPO_r_dist1_3        (EXP, Negbin)
%   @options_forecast_PII_GGM_r_p_dist1_1       (GGM, Poisson)
%   @options_forecast_PII_GGM_r_p_dist1_3       (GGM, Negbin)
%   @options_forecast_PII_GLM_r_p_K_dist1_1     (GLM, Poisson)
%   @options_forecast_PII_GLM_r_p_K_dist1_3     (GLM, Negbin)
%   @options_forecast_PII_SIR_beta_dist1_1      (SIR S1, Poisson)
%   @options_forecast_PII_SIR_beta_dist1_3      (SIR S1, Negbin)
%   @options_forecast_PII_SIR_beta_gamma_dist1_3          (SIR S2)
%   @options_forecast_PII_SEIR_beta_dist1_3               (SEIR S1)
%   @options_forecast_PII_SEIR_beta_gamma_dist1_3         (SEIR S2)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3   (SEIR S3)
%   @options_forecast_PII_SEIRD_beta_dist1_3              (SEIRD S1)
%   @options_forecast_PII_SEIRD_beta_rho_dist1_3          (SEIRD S2)
%   @options_forecast_PII_SEIRD_beta_rho_gamma_dist1_3    (SEIRD S3)
%   @options_forecast_PII_SEIRunrep_beta_dist1_3          (SEIR_unrep S1)
%   @options_forecast_PII_SEIRunrep_beta_rho_dist1_3      (SEIR_unrep S2)
%   @options_forecast_PII_SEIRunrep_beta_rho_gamma_dist1_3 (SEIR_unrep S3)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_dist1_3            (SEIR_asymp S1)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_rho_dist1_3        (SEIR_asymp S2)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_rho_gamma_dist1_3  (SEIR_asymp S3)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_2vars  (SEIR_vars 2vars)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_3vars  (SEIR_vars 3vars)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_dCdt   (SEIR_vars dCdt)
options_handle = @options_forecast_PII_EXPO_r_dist1_3;

% Error structure label (for plot titles / PDF filenames):
%   'Poisson', 'Negbin5', or 'Negbin10'
error_type = 'Negbin5';

% Calibration window size — scalar or array
%   scalar:  windowsize1 = 20;
%   array:   windowsize1 = 20:10:100;
windowsize1 = 20:10:50;

% Number of replicates
num_replicates = 500;

% Whether to run the computation (true) or skip to plotting (false)
run_flag = true;

% Plot type (after computation or from existing output):
%   'PII'     — PII median + 95% CI vs calibration period
%   'CI_grid' — 95% CI caterpillar plot per parameter
%   'both'    — generate both PII and CI_grid plots
%   'none'    — no plotting (run only)
plot_type = 'both';
