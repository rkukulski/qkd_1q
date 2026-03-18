# Validation Harness Outputs

This folder contains the Sprint 1 validation pipeline for manuscript-ready artifacts.

## Run Validation

```bash
julia --project=. src/analysis/run_validation.jl \
  --output-dir results/validation \
  --eps-min 0.0 \
  --eps-max 0.25 \
  --eps-step 0.01 \
  --restarts 3 \
  --tol 1e-7 \
  --seed 20260318 \
  --suppress-warnings true
```

Generated files:

- `validation_points.csv`
- `validation_summary.csv`
- `validation_meta.json`

## Generate Plots

```bash
julia --project=. src/analysis/plot_validation.jl \
  results/validation/validation_points.csv \
  results/validation
```

Generated files:

- `validation_curves.pdf`, `validation_curves.png`
- `validation_gaps.pdf`, `validation_gaps.png`

## CSV/JSON Contract

`validation_points.csv` columns:

`epsilon,case1_R,case2_R,case3_R,upper_bound,gap_case1,gap_case2,gap_case3,best_case,case1_x_opt,case2_x_opt,case2_y_opt,case1_best_seed,case2_best_seed`

`validation_summary.csv` columns:

`row_type,case,epsilon,mean,std,min,max,best_seed,best_R,crossover_from,crossover_to,crossover_est_eps,bracket_eps_left,bracket_eps_right`

`validation_meta.json` keys:

`created_at,script,eps_min,eps_max,eps_step,eps_count,restarts,tol,seed,suppress_warnings,git_commit,runtime_seconds,output_dir,points_csv,summary_csv,meta_file`
