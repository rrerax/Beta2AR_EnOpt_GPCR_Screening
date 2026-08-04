from __future__ import annotations

import argparse
import math
import urllib.request
from pathlib import Path

from common import PROJECT_ROOT, ensure_dir, load_config, receptor_entries, run_command, write_csv


def download_pdb(pdb_id: str, output_path: Path) -> None:
    if output_path.exists() and output_path.stat().st_size > 0:
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    with urllib.request.urlopen(url, timeout=90) as response:
        output_path.write_bytes(response.read())


def keep_receptor_chain(raw_pdb: Path, output_pdb: Path, chain_id: str | None) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for line in raw_pdb.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain_id and line[21].strip() != chain_id:
            continue
        alt_loc = line[16].strip()
        if alt_loc not in {"", "A"}:
            continue
        lines.append(line[:16] + " " + line[17:] if alt_loc == "A" else line)
    if not lines:
        raise ValueError(f"No ATOM records kept for {raw_pdb.name}, chain={chain_id}")
    output_pdb.write_text("\n".join(lines) + "\nTER\nEND\n", encoding="utf-8")


def ligand_coordinates(raw_pdb: Path, chain_id: str | None, resname: str | None) -> list[tuple[float, float, float]]:
    if not resname:
        return []
    coords: list[tuple[float, float, float]] = []
    for line in raw_pdb.read_text(errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        if chain_id and line[21].strip() != chain_id:
            continue
        if line[17:20].strip() != resname:
            continue
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def box_from_native_ligand(coords: list[tuple[float, float, float]], default_size: list[float], padding: float) -> tuple[list[float], list[float]]:
    if not coords:
        return [0.0, 0.0, 0.0], list(default_size)
    mins = [min(point[i] for point in coords) for i in range(3)]
    maxs = [max(point[i] for point in coords) for i in range(3)]
    center = [(mins[i] + maxs[i]) / 2 for i in range(3)]
    ligand_span = [maxs[i] - mins[i] + 2 * padding for i in range(3)]
    size = [max(float(default_size[i]), ligand_span[i]) for i in range(3)]
    return center, size


def prepare_pdbqt(input_pdb: Path, output_pdbqt: Path, force: bool) -> str:
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 0 and not force:
        return "ready"
    log_path = PROJECT_ROOT / "logs" / "receptors" / f"{output_pdbqt.stem}.log"
    command = ["obabel", str(input_pdb), "-O", str(output_pdbqt), "-xr", "-h"]
    return_code = run_command(command, log_path=log_path)
    if return_code != 0 or not output_pdbqt.exists() or output_pdbqt.stat().st_size == 0:
        return "pdb_ready_pdbqt_failed"
    return "ready"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/beta2ar_screen.yml")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    raw_dir = ensure_dir("data/raw/receptors")
    prepared_dir = ensure_dir("data/processed/receptors")
    ensure_dir("logs/receptors")

    default_size = [float(value) for value in config["receptors"]["default_box_size_angstrom"]]
    padding = float(config["receptors"].get("box_padding_angstrom", 8.0))

    rows: list[dict[str, object]] = []
    for entry in receptor_entries(config):
        pdb_id = str(entry["pdb_id"]).upper()
        conformation_id = str(entry.get("conformation_id") or pdb_id)
        chain_id = entry.get("receptor_chain")
        native_ligand = entry.get("native_ligand_resname")
        raw_pdb = raw_dir / f"{pdb_id}.pdb"
        download_pdb(pdb_id, raw_pdb)

        prepared_pdb = prepared_dir / f"{conformation_id}.pdb"
        prepared_pdbqt = prepared_dir / f"{conformation_id}.pdbqt"
        if args.force or not prepared_pdb.exists():
            keep_receptor_chain(raw_pdb, prepared_pdb, chain_id)
        coords = ligand_coordinates(raw_pdb, chain_id, native_ligand)
        center, size = box_from_native_ligand(coords, default_size, padding)
        status = prepare_pdbqt(prepared_pdb, prepared_pdbqt, force=args.force)
        rows.append(
            {
                "conformation_id": conformation_id,
                "pdb_id": pdb_id,
                "state_label": entry.get("state_label", "unspecified"),
                "receptor_chain": chain_id or "",
                "native_ligand_resname": native_ligand or "",
                "native_ligand_atoms": len(coords),
                "center_x": round(center[0], 3),
                "center_y": round(center[1], 3),
                "center_z": round(center[2], 3),
                "size_x": round(size[0], 3),
                "size_y": round(size[1], 3),
                "size_z": round(size[2], 3),
                "prepared_pdb": str(prepared_pdb.relative_to(PROJECT_ROOT)),
                "receptor_pdbqt": str(prepared_pdbqt.relative_to(PROJECT_ROOT)),
                "status": status,
            }
        )

    out = write_csv(
        "results/tables/receptor_manifest.csv",
        rows,
        [
            "conformation_id",
            "pdb_id",
            "state_label",
            "receptor_chain",
            "native_ligand_resname",
            "native_ligand_atoms",
            "center_x",
            "center_y",
            "center_z",
            "size_x",
            "size_y",
            "size_z",
            "prepared_pdb",
            "receptor_pdbqt",
            "status",
        ],
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

