# β2-Adrenergic Receptor Ensemble Virtual Screening

This repository contains a reproducible computational screening workflow for the β2-adrenergic receptor (β2-AR, ADRB2). The workflow uses multiple receptor conformations to account for GPCR conformational flexibility and applies an EnOpt-style weighted consensus score to prioritize candidate ligands from a 30,000-compound ChEMBL library.

## Project Overview

Classical high-throughput virtual screening often docks each ligand against a single receptor structure. For GPCR targets, this can be limiting because transmembrane helix rearrangements, especially around TM6 activation-associated motion, change the binding pocket geometry. This project therefore uses an ensemble of inactive and active β2-AR conformations, then combines docking results into a single ranking matrix.

## Dataset Summary

- Ligand source: ChEMBL37 subset generated for β2-AR screening
- Input ligand records: 30,000
- Prepared ligand structures: 29,865
- Receptor conformations: 5
- Completed docking score rows: 149,325
- Ranking output: top 200 compounds by EnOpt-style weighted score

Ligand preparation status:

- `ready`: 29,865
- `pdbqt_failed`: 109
- `embed_failed`: 20
- `prep_timeout`: 6

## Receptor Ensemble

The receptor set includes inactive and active-state β2-AR structures so that the screen is not tied to a single rigid receptor geometry.

| conformation_id         |   weight |
|:------------------------|---------:|
| active_bi167107_4LDE    |   0.198  |
| active_gs_complex_3SN6  |   0.1943 |
| active_hbi_4LDL         |   0.2018 |
| inactive_carazolol_2RH1 |   0.2031 |
| inactive_carazolol_5D5A |   0.2028 |

## Methods

1. **Ligand library construction**: A 30,000-compound ChEMBL37 ligand table was prepared with compound identifiers and SMILES strings.
2. **Protein structure preparation**: Five β2-AR crystal structures were selected to represent active and inactive conformational states.
3. **Ligand preparation**: Ligands were desalted, embedded into 3D conformers, minimized, and converted into docking-ready PDBQT files.
4. **Batch docking**: AutoDock Vina was used to dock each prepared ligand against each receptor conformation.
5. **Score matrix construction**: Docking scores were reshaped into a ligand-by-conformation matrix.
6. **EnOpt-style ranking**: Conformation weights were estimated from the score ensemble and used to calculate weighted consensus scores.
7. **Result review**: Top-ranked molecules were summarized in tables and visualized with score distribution and top-hit figures.

## Key Results

Top-ranked compounds from the ensemble weighted score:

|   rank | ligand_id     |   weighted_score |   best_score |   mean_score |   score_sd |   n_conf |
|-------:|:--------------|-----------------:|-------------:|-------------:|-----------:|---------:|
|      1 | CHEMBL3311247 |          -13.164 |       -14.48 |      -13.146 |      1.249 |        5 |
|      2 | CHEMBL16965   |          -12.896 |       -13.35 |      -12.888 |      0.553 |        5 |
|      3 | CHEMBL36113   |          -12.787 |       -13.71 |      -12.772 |      1.06  |        5 |
|      4 | CHEMBL4436402 |          -12.754 |       -13.65 |      -12.738 |      1.065 |        5 |
|      5 | CHEMBL33607   |          -12.721 |       -13.48 |      -12.706 |      1.007 |        5 |
|      6 | CHEMBL12018   |          -12.641 |       -13.91 |      -12.628 |      1.023 |        5 |
|      7 | CHEMBL12143   |          -12.63  |       -13.25 |      -12.62  |      0.687 |        5 |
|      8 | CHEMBL26449   |          -12.616 |       -13.43 |      -12.604 |      0.951 |        5 |
|      9 | CHEMBL11778   |          -12.571 |       -13.33 |      -12.558 |      1.055 |        5 |
|     10 | CHEMBL21788   |          -12.564 |       -13.67 |      -12.55  |      1.013 |        5 |

Main figures:

- `results/figures/enopt_weighted_top_hits.png`
- `results/figures/score_distributions.png`

## Repository Structure

```text
configs/                    Workflow configuration
data/raw/                   Input ligand library
scripts/                    Reproducible preparation, docking, and analysis scripts
notebooks/                  Analysis notebook for result review
docs/                       Method notes and final run summary
results/tables/             Docking scores, score matrix, weights, and top hits
results/figures/            Publication-style result figures
examples/top_hit_poses/     Representative top-hit docking outputs
```

## Reproducing the Workflow

Create the environment and run a smoke test:

```bash
bash scripts/00_setup_env.sh
bash scripts/run_smoke_test.sh
```

Run the full workflow:

```bash
bash scripts/run_pipeline.sh
```

The full run is CPU-intensive because it performs docking across 30,000 ligands and five receptor conformations. The supplied results were generated with parallel CPU docking.

## Notes and Limitations

- Docking scores are computational prioritization signals, not experimental binding affinities.
- The ranking should be interpreted as a shortlist for follow-up analysis, not as confirmed biological activity.
- EnOpt-style weighting improves ensemble-level interpretation, but additional validation such as redocking controls, molecular dynamics, or experimental assays would be required before biological claims.
