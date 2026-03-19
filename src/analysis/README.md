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

## Sprint 2A: N=1 Theorem Draft Validation

```bash
julia --project=. src/analysis/run_n1_theorem_validation.jl \
  --output-dir results/sprint2_n1 \
  --eps-min 0.0 \
  --eps-max 0.25 \
  --eps-step 0.005 \
  --restarts 5 \
  --tol 1e-7 \
  --seed 20260319 \
  --suppress-warnings true \
  --u-min -3.0 \
  --u-max 1.0 \
  --u-step 0.01
```

Generated files:

- `n1_points.csv`
- `n1_summary.csv`
- `n1_meta.json`
- `upper_bound_audit.csv`

Optional plots:

```bash
julia --project=. src/analysis/plot_n1_theorem_validation.jl \
  results/sprint2_n1/n1_points.csv \
  results/sprint2_n1
```

Generated files:

- `n1_theorem_curves.pdf`, `n1_theorem_curves.png`
- `n1_theorem_mismatch.pdf`, `n1_theorem_mismatch.png`

## Sprint 2B: N=1 Proof Diagnostics

```bash
julia --project=. src/analysis/run_n1_proof_diagnostics.jl \
  --input-points results/sprint2_n1/n1_points.csv \
  --output-dir results/sprint2_n1 \
  --local-u-step 0.01 \
  --x-sensitivity-step 0.001 \
  --x-sensitivity-tol 1e-8 \
  --x-sensitivity-strict-drop 1e-5 \
  --lambda-tol 1e-6 \
  --tol 1e-7 \
  --suppress-warnings true
```

Generated files:

- `n1_proof_diagnostics.csv`
- `n1_proof_diagnostics_summary.csv`

## Script-to-Manuscript Mapping

Use this section as the single source of truth for where each analysis script
lands in `txt/darlings/main.tex`.

| Script | Primary outputs | Manuscript targets |
|---|---|---|
| `src/analysis/run_validation.jl` | `results/validation/validation_points.csv`, `validation_summary.csv`, `validation_meta.json` | Subsection `Validation summary (production grid)`; tables `tab:validation-case-summary`, `tab:validation-crossovers`; appendix table `tab:validation-grid-full-appendix` |
| `src/analysis/plot_validation.jl` | `results/validation/validation_curves.pdf/png`, `validation_gaps.pdf/png` | Figure `fig:validation-curves` (main body) |
| `src/analysis/run_n1_theorem_validation.jl` | `results/sprint2_n1/n1_points.csv`, `n1_summary.csv`, `n1_meta.json` | Subsection `N=1 optimality theorem (draft)`; table `tab:n1-theorem-summary`; appendix table `tab:n1-theorem-grid-full-appendix` |
| `src/analysis/run_upper_bound_audit.jl` | `results/sprint2_n1/upper_bound_audit.csv` | Referenced in N=1 theorem subsection as bound-dominance evidence |
| `src/analysis/plot_n1_theorem_validation.jl` | `results/sprint2_n1/n1_theorem_curves.pdf/png`, `n1_theorem_mismatch.pdf/png` | Figure `fig:n1-theorem-curves` |
| `src/analysis/run_n1_proof_diagnostics.jl` | `results/sprint2_n1/n1_proof_diagnostics.csv`, `n1_proof_diagnostics_summary.csv` | KKT/Lemma 2 evidence table `tab:n1-proof-diagnostics` in N=1 theorem subsection |

## Label Index in `main.tex`

- `Numerical artifact contract`
- `tab:n1-theorem-summary`
- `tab:n1-proof-diagnostics`
- `fig:n1-theorem-curves`
- `tab:validation-case-summary`
- `tab:validation-crossovers`
- `fig:validation-curves`
- `tab:validation-grid-full-appendix`
- `tab:n1-theorem-grid-full-appendix`

## Recommended Regeneration Order

1. `run_validation.jl`
2. `plot_validation.jl`
3. `run_n1_theorem_validation.jl`
4. `plot_n1_theorem_validation.jl`
5. `run_n1_proof_diagnostics.jl`
