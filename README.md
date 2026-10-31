# psoriasis — mechanism-aware variant effect prediction for the psoriasis spectrum

Pilot experiments and a paper proposal for disease-specific, **mechanism-aware** variant effect
prediction (VEP) across the psoriasis spectrum: plaque psoriasis, generalised pustular psoriasis
(GPP), palmoplantar pustulosis, acrodermatitis continua of Hallopeau, pityriasis rubra pilaris and
`CARD14`-associated papulosquamous eruption.

## The one-paragraph version

Generic missense predictors look excellent on the psoriasis spectrum by AUROC (0.90–0.98) and fail
where it counts. At AlphaMissense's own calibrated threshold, it recovers **85% of
loss-of-function pathogenic variants but only 42% of gain-of-function ones** — and four of the five
founding `CARD14` psoriasis alleles, including p.Gly117Ser (score 0.236), are called *likely
benign*. Simply conditioning the decision threshold on disease mechanism nearly doubles
gain-of-function sensitivity at matched false-positive rate (0.258 → 0.495). Separately, **one in
six** psoriasis-spectrum variants flips its clinical class call depending on which transcript is
scored. And the single most common GPP allele in Vietnam (`IL36RN` c.115+6T>C, 44% of patients) is
a splice-region variant that no missense predictor scores at all.

Full numbers, methods, limitations and the system design: **[`docs/technical-report.md`](docs/technical-report.md)**.
System diagram also as [`docs/system-design.svg`](docs/system-design.svg).

## Layout

```
src/        experiment code (stdlib + numpy + scikit-learn only)
results/    metrics as JSON and the raw run logs
docs/       technical report and system design diagram
```

| Script | What it does |
|---|---|
| `src/extract_clinvar.py` | pulls ClinVar records for the 30-gene panel in `src/mech.json` |
| `src/join_am.py` | joins them to canonical AlphaMissense scores |
| `src/join_am_all.py` | same, keeping every transcript (used to confirm the canonical table is one row per variant) |
| `src/analyze.py` | Experiment 1: mechanism-stratified performance, VUS burden, named literature alleles |
| `src/exp2_calib_cluster.py` | Experiment 2: mechanism-conditional calibration. Experiment 3: positional clustering |
| `src/exp4_isoform.py` | Experiment 4: cross-isoform score instability |

Data is not committed (ClinVar 147 MB, AlphaMissense 2 GB). Download commands are in §7 of the
technical report.

## Status

Exploratory. **No model has been trained.** Every number comes from precomputed predictor scores —
the session that produced this had no GPU and no route to model weights. Experiments 1, 3 and 4
measure how the current baseline behaves; Experiment 2 is a recalibration, not a new model. The
limitations are listed in full in §6 of the report, and they should be read before any of the
numbers are quoted.

## Licensing of inputs

ClinVar is public domain. AlphaMissense predictions are CC BY-NC-SA 4.0 (non-commercial use only);
no AlphaMissense data is redistributed in this repository.
