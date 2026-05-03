# Practical Identifiability Index (PII)

This repository accompanies the manuscript **"A Practical Identifiability Index for Dynamical Models"** by Hamed Karami, Alexandra Smirnova, Sunmi Lee, and Gerardo Chowell.

The project introduces and evaluates a Practical Identifiability Index (PII) for ordinary differential equation (ODE) models. PII summarizes the width of a parameter confidence interval on a base-10 logarithmic scale, making it easier to compare practical identifiability across parameters, models, observation windows, and error structures.

For a parameter `theta_j` with 95% confidence interval `[theta_L, theta_U]`, the manuscript defines

```text
PII_j = log10(theta_U / (theta_L + epsilon)),   epsilon = 1e-3.
```

Small PII values indicate tightly constrained parameters, while larger values indicate weak practical identifiability.

## Current Identifiability Scale

The current manuscript uses the following interpretation:

| PII range | Classification |
| --- | --- |
| `PII < 0.1` | Identifiable |
| `0.1 <= PII < 1` | Weakly identifiable |
| `PII >= 1` | Non-identifiable |
| `PII < 0.1` and empirical 95% CI coverage `< 90%` | Identifiable with low coverage |

Across the 135 model, scenario, error-structure, and parameter combinations summarized in the manuscript, the current results classify:

- 77 combinations as identifiable.
- 27 combinations as identifiable with low coverage.
- 30 combinations as weakly identifiable.
- 1 combination as non-identifiable: SEIAR Scenario 3, parameter `beta_1`, under Poisson error.

## Study Design

The manuscript evaluates PII using synthetic-data experiments where the true parameter values are known. For each model, scenario, error structure, and calibration window:

1. A noiseless ODE trajectory is generated.
2. Observation noise is added.
3. The model is fitted to the first `T` observations.
4. Parametric bootstrap confidence intervals are computed.
5. PII is calculated from the bootstrap confidence interval span.
6. Results are summarized across 500 independent replicates.

The manuscript considers three observation-error structures:

- `Poisson`
- `Negbin5`: negative binomial with dispersion parameter `alpha = 5`
- `Negbin10`: negative binomial with dispersion parameter `alpha = 10`

The model set includes growth models and compartmental epidemic models:

- `EXP`: exponential growth
- `GGM`: generalized growth model
- `GLM`: generalized logistic model
- `SIR`
- `SEIR`
- `SEIR-UR`: SEIR with under-reporting
- `SEIAR`: SEIR with asymptomatic transmission
- `SEIRD`: SEIR with disease-induced mortality
- `SEIRMO`: multi-observable SEIR analysis using `dC/dt`, `I`, and `R`


|-- Final_code_PII/
|-- Final_code_PII_cluster/
|-- Options
```

## Requirements

MATLAB requirements:

- MATLAB with ODE solvers such as `ode45` and `ode15s`.
- Optimization Toolbox, used through `fmincon` and optimization problem setup.
- Global Optimization Toolbox, used through `MultiStart` and `CustomStartPointSet`.
- Statistics and Machine Learning Toolbox, used for functions such as `poissrnd`, `nbinrnd`, `lhsdesign`, and `quantile`.
- Parallel Computing Toolbox is recommended for the local `parfor` workflow.

## Local MATLAB Workflow

Use `Final_code_PII/` for local or workstation runs.

1. Open `Final_code_PII/plotPracticalIdentifiabilityResults.m`.
2. Edit the user settings at the top of the file:

```matlab
options_handle = @options_forecast_PII_EXPO_r_dist1_3;
error_type = 'Negbin5';
windowsize1 = 30:10:50;
num_replicates = 20;
run_flag = true;
plot_type = 'both';
```

3. Make sure the selected model function and option file are on the MATLAB path.
4. Run the script from inside `Final_code_PII/`.

```matlab
plotPracticalIdentifiabilityResults
```

The script writes output under `Final_code_PII/output/`. For full manuscript-scale runs, increase `num_replicates` to 500 and use the calibration windows listed in the manuscript.

When switching negative-binomial settings, keep the displayed `error_type` and the overdispersion factor passed to `Run_PracticalIndentifiability_ODEModel` synchronized:

```text
Poisson  -> factor1 = 1
Negbin5  -> factor1 = 5
Negbin10 -> factor1 = 10
```

## Cluster MATLAB Workflow

Use `Final_code_PII_cluster/` for HPC runs.

The cluster workflow is controlled by `Final_code_PII_cluster/main.m`. Each SLURM array task calls

```matlab
main(taskID)
```

where `taskID` maps to one `(calibration window, replicate)` pair. The total number of tasks is

```text
length(windowsize_array) * num_replicates
```

For example, if

```matlab
windowsize_array = 20:10:50;
num_replicates = 500;
```

then the job array should contain `4 * 500 = 2000` tasks.

Minimal SLURM-style call:

```bash
matlab -batch "main(${SLURM_ARRAY_TASK_ID})"
```

Before submitting jobs, edit the user settings at the top of `Final_code_PII_cluster/main.m`:

```matlab
options_handle = @options_forecast_PII_EXPO_r_dist1_3;
factor1 = 5;
windowsize_array = 20:10:50;
num_replicates = 500;
```

The cluster workflow writes one `.mat` file per replicate and calibration window under `Final_code_PII_cluster/output/`, unless the script is copied into a model-specific run directory with its own `output/` folder.

## Option Files and Scenarios

Option files encode the model, observation operator, estimated parameters, true values, bounds, bootstrap size, and error model. The naming convention is:

```text
options_forecast_PII_<MODEL>_<estimated_parameters>_dist1_<error_id>.m
```

where `dist1_1` denotes Poisson and `dist1_3` denotes negative binomial.

Examples:

| Option handle | Meaning |
| --- | --- |
| `@options_forecast_PII_EXPO_r_dist1_1` | EXP, estimate `r`, Poisson |
| `@options_forecast_PII_EXPO_r_dist1_3` | EXP, estimate `r`, negative binomial |
| `@options_forecast_PII_GGM_r_p_dist1_3` | GGM, estimate `r` and `p` |
| `@options_forecast_PII_GLM_r_p_K_dist1_3` | GLM, estimate `r`, `p`, and `K` |
| `@options_forecast_PII_SIR_beta_gamma_dist1_3` | SIR Scenario 2, estimate `beta` and `gamma` |
| `@options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3` | SEIR Scenario 3, estimate `beta`, `kappa`, and `gamma` |
| `@options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_dCdt` | SEIRMO Scenario 1, observe `dC/dt` |
| `@options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_2vars` | SEIRMO Scenario 2, observe `dC/dt` and `I` |
| `@options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_3vars` | SEIRMO Scenario 3, observe `dC/dt`, `I`, and `R` |



## Code Availability

```text
https://github.com/hkarami-GSU/PII/tree/main
```
