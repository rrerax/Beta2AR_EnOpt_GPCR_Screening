# Final Run Summary

## Run Scope

- Target: β2-adrenergic receptor (β2-AR / ADRB2)
- Ligand library: 30,000 ChEMBL37 records
- Receptor conformations: 5
- Prepared ligands used for docking: 29,865
- Docking score rows generated: 149,325
- Final ranked output: top 200 compounds

## Ligand Preparation Status

- `ready`: 29,865
- `pdbqt_failed`: 109
- `embed_failed`: 20
- `prep_timeout`: 6

## Ensemble Weights

| conformation_id         |   weight |
|:------------------------|---------:|
| active_bi167107_4LDE    |   0.198  |
| active_gs_complex_3SN6  |   0.1943 |
| active_hbi_4LDL         |   0.2018 |
| inactive_carazolol_2RH1 |   0.2031 |
| inactive_carazolol_5D5A |   0.2028 |

## Top 10 Ranked Compounds

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

## Output Files

- `results/tables/docking_scores.csv`
- `results/tables/enopt_style_score_matrix.csv`
- `results/tables/enopt_style_top_hits.csv`
- `results/tables/conformation_weights.csv`
- `results/tables/ligand_manifest.csv`
- `results/tables/receptor_manifest.csv`
- `results/figures/enopt_weighted_top_hits.png`
- `results/figures/score_distributions.png`

## Interpretation

The complete run satisfies the planned second-stage workflow: β2-AR receptor preparation, 30,000-compound ligand processing, five-conformation docking, ensemble score integration, EnOpt-style weighted ranking, result tables, figures, and reproducible scripts. The results are suitable for a computational portfolio project and should be described as virtual screening prioritization rather than experimentally confirmed activity.
