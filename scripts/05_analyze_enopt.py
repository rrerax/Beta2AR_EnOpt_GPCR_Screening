from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import ensure_dir, load_config, project_path


def conformation_weights(pivot: pd.DataFrame) -> pd.Series:
    consensus = pivot.mean(axis=1, skipna=True)
    weights: dict[str, float] = {}
    for column in pivot.columns:
        paired = pd.concat([pivot[column], consensus], axis=1).dropna()
        if len(paired) < 10:
            weights[column] = 0.0
            continue
        rho, _ = spearmanr(paired.iloc[:, 0], paired.iloc[:, 1])
        weights[column] = max(float(rho) if np.isfinite(rho) else 0.0, 0.0)
    weight_series = pd.Series(weights, dtype=float)
    if weight_series.sum() <= 0:
        weight_series[:] = 1.0 / max(len(weight_series), 1)
    else:
        weight_series = weight_series / weight_series.sum()
    return weight_series


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/beta2ar_screen.yml")
    parser.add_argument("--score-file", default="results/tables/docking_scores.csv")
    parser.add_argument("--output-prefix", default="")
    args = parser.parse_args()

    config = load_config(args.config)
    top_n = int(config["analysis"].get("top_n", 200))
    docking = pd.read_csv(project_path(args.score_file))
    docking = docking[docking["status"] == "ready"].copy()
    docking["best_score_kcal_mol"] = pd.to_numeric(docking["best_score_kcal_mol"], errors="coerce")
    docking = docking.dropna(subset=["best_score_kcal_mol"])
    if docking.empty:
        raise ValueError("No completed docking scores found.")

    pivot = docking.pivot_table(
        index="ligand_id",
        columns="conformation_id",
        values="best_score_kcal_mol",
        aggfunc="min",
    )
    weights = conformation_weights(pivot)
    weighted = pivot.mul(weights, axis=1).sum(axis=1, skipna=True)
    matrix = pd.DataFrame(
        {
            "ligand_id": pivot.index,
            "docked_conformation_count": pivot.notna().sum(axis=1).values,
            "ensemble_best_score_kcal_mol": pivot.min(axis=1, skipna=True).values,
            "ensemble_mean_score_kcal_mol": pivot.mean(axis=1, skipna=True).values,
            "ensemble_median_score_kcal_mol": pivot.median(axis=1, skipna=True).values,
            "ensemble_score_std_kcal_mol": pivot.std(axis=1, skipna=True).fillna(0.0).values,
            "enopt_style_weighted_score_kcal_mol": weighted.values,
        }
    )
    ligands = pd.read_csv(project_path(config["ligands"]["input_csv"]))
    matrix = matrix.merge(ligands, on="ligand_id", how="left")
    matrix = matrix.sort_values("enopt_style_weighted_score_kcal_mol", ascending=True)

    table_dir = ensure_dir("results/tables")
    figure_dir = ensure_dir("results/figures")
    prefix = args.output_prefix
    matrix_path = table_dir / f"{prefix}enopt_style_score_matrix.csv"
    top_path = table_dir / f"{prefix}enopt_style_top_hits.csv"
    weights_path = table_dir / f"{prefix}conformation_weights.csv"
    matrix.to_csv(matrix_path, index=False)
    matrix.head(top_n).to_csv(top_path, index=False)
    weights.rename("weight").reset_index().rename(columns={"index": "conformation_id"}).to_csv(weights_path, index=False)

    plt.figure(figsize=(10, 5))
    for conformation_id, group in docking.groupby("conformation_id"):
        group["best_score_kcal_mol"].plot(kind="kde", label=conformation_id)
    plt.xlabel("Vina score (kcal/mol)")
    plt.ylabel("Density")
    plt.title("β2-AR docking score distributions across receptor conformations")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_dir / f"{prefix}score_distributions.png", dpi=300)
    plt.close()

    top30 = matrix.head(30).sort_values("enopt_style_weighted_score_kcal_mol", ascending=True)
    plt.figure(figsize=(10, 8))
    plt.barh(top30["ligand_id"], top30["enopt_style_weighted_score_kcal_mol"], color="#4C78A8")
    plt.xlabel("EnOpt-style weighted score (kcal/mol)")
    plt.ylabel("Ligand")
    plt.title("Top β2-AR candidates after ensemble reranking")
    plt.tight_layout()
    plt.savefig(figure_dir / f"{prefix}enopt_weighted_top_hits.png", dpi=300)
    plt.close()

    summary_path = project_path(f"docs/{prefix}analysis_summary.md")
    top_hit = matrix.iloc[0]
    summary_path.write_text(
        "# β2-AR Ensemble Screening Summary\n\n"
        f"- Completed docking score rows: {len(docking):,}\n"
        f"- Ligands with at least one conformation score: {len(matrix):,}\n"
        f"- Receptor conformations used: {pivot.shape[1]}\n"
        f"- Top ranked ligand: {top_hit['ligand_id']}\n"
        f"- Weighted docking score: {top_hit['enopt_style_weighted_score_kcal_mol']:.3f} kcal/mol\n\n"
        "The weighted score is used as a computational prioritization metric. It is not an experimental affinity value.\n",
        encoding="utf-8",
    )
    print(f"Wrote {matrix_path}")
    print(f"Wrote {top_path}")
    print(f"Wrote {weights_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
