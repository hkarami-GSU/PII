function main(taskID)
% ============================================================================
% CLUSTER VERSION — each SLURM task calls main(taskID)
% ============================================================================
% The taskID (integer) maps to a unique (windowsize, replicate) pair.
% Submit as a job array with taskID = 1 : num_windowsizes * num_replicates.
%
% Example SLURM submission: EXP: windowsize:20,30,40,50-num_replicates=500
% array=1-2000 -- meaning 4 (length windowsize) * 500 (num_replicates).
% ============================================================================

% ============================================================================
% USER SETTINGS — Change these as needed
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

% Error structure overdispersion parameter for Run_PracticalIndentifiability_ODEModel:
%   1  for Poisson
%   5  for Negbin5
%   10 for Negbin10
factor1 = 5;

% Calibration window sizes — array of values to sweep
windowsize_array = 20:10:50;

% Number of replicates
num_replicates = 500;

% ============================================================================
% END USER SETTINGS
% ============================================================================

rng(taskID);

num_windowsizes = length(windowsize_array);

windowsize_idx = ceil(taskID / num_replicates);
replicate_idx = mod(taskID - 1, num_replicates) + 1;

current_windowsize = windowsize_array(windowsize_idx);

[~,~,~,~,~,~,model_INP,~,~,~,~,~,~,~,~] = options_handle();

outfile = sprintf('./output/results-replicate-%d-model_name-%s-calibrationperiod-%d.mat', ...
    replicate_idx, model_INP.name, current_windowsize);
if isfile(outfile)
    fprintf('Output exists: %s - skipping.\n', outfile);
    return
end

Run_PracticalIndentifiability_ODEModel(...
    options_handle,...
    current_windowsize, factor1, 1, replicate_idx);
end
