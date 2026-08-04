from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from common import PROJECT_ROOT, load_config, project_path, receptor_entries


def main() -> None:
    config = load_config()
    ligand_csv = project_path(config["ligands"]["input_csv"])
    if not ligand_csv.exists():
        raise FileNotFoundError(f"Missing ligand input: {ligand_csv}")

    ligands = pd.read_csv(ligand_csv)
    ligand_count = len(ligands)
    expected = int(config["ligands"].get("full_limit", 30000))
    if ligand_count < expected:
        raise ValueError(f"Ligand library has {ligand_count} rows, expected at least {expected}")

    missing_receptors: list[str] = []
    for entry in receptor_entries(config):
        pdb_path = PROJECT_ROOT / "data" / "raw" / "receptors" / f"{entry['pdb_id']}.pdb"
        if not pdb_path.exists():
            missing_receptors.append(str(pdb_path))

    if missing_receptors:
        raise FileNotFoundError("Missing receptor PDB files:\n" + "\n".join(missing_receptors))

    tool_status = {
        "vina": shutil.which("vina") or "not found",
        "obabel": shutil.which("obabel") or "not found",
    }
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Ligand rows: {ligand_count}")
    print(f"Receptor conformations: {len(receptor_entries(config))}")
    for tool_name, tool_path in tool_status.items():
        print(f"{tool_name}: {tool_path}")


if __name__ == "__main__":
    main()

