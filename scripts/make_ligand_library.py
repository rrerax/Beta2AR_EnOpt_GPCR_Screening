from __future__ import annotations

import argparse
import csv
import gzip
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")


CHEMBL_CHEMREPS_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz"


def largest_fragment(mol: Chem.Mol) -> Chem.Mol | None:
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return None
    return max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())


def clean_smiles(smiles: str) -> tuple[str, Chem.Mol] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fragment = largest_fragment(mol)
    if fragment is None:
        return None
    canonical = Chem.MolToSmiles(fragment, isomericSmiles=True, canonical=True)
    clean_mol = Chem.MolFromSmiles(canonical)
    if clean_mol is None:
        return None
    return canonical, clean_mol


def passes_filters(mol: Chem.Mol) -> bool:
    molecular_weight = Descriptors.MolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    return (
        150.0 <= molecular_weight <= 550.0
        and 10 <= heavy_atoms <= 45
        and logp <= 6.0
        and hbd <= 8
        and hba <= 12
    )


def descriptors(mol: Chem.Mol) -> dict[str, float | int]:
    return {
        "mol_weight": round(float(Descriptors.MolWt(mol)), 3),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=30000)
    parser.add_argument("--output", default="data/raw/chembl37_ligand_library_30000.csv")
    parser.add_argument("--url", default=CHEMBL_CHEMREPS_URL)
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ligand_id",
        "smiles",
        "source",
        "mol_weight",
        "heavy_atoms",
        "logp",
        "hbd",
        "hba",
        "tpsa",
    ]
    seen_smiles: set[str] = set()
    kept = 0
    scanned = 0

    with urllib.request.urlopen(args.url, timeout=120) as response, output_path.open("w", newline="", encoding="utf-8") as handle:
        gz_handle = gzip.GzipFile(fileobj=response)
        text_iter = (line.decode("utf-8", "replace") for line in gz_handle)
        reader = csv.DictReader(text_iter, delimiter="\t")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            scanned += 1
            chembl_id = row.get("chembl_id", "").strip()
            smiles = row.get("canonical_smiles", "").strip()
            if not chembl_id or not smiles:
                continue
            cleaned = clean_smiles(smiles)
            if cleaned is None:
                continue
            canonical_smiles, mol = cleaned
            if canonical_smiles in seen_smiles or not passes_filters(mol):
                continue
            seen_smiles.add(canonical_smiles)
            output_row = {
                "ligand_id": chembl_id,
                "smiles": canonical_smiles,
                "source": "ChEMBL37_chemreps",
                **descriptors(mol),
            }
            writer.writerow(output_row)
            kept += 1
            if kept % 5000 == 0:
                print(f"kept={kept} scanned={scanned}", flush=True)
            if kept >= args.limit:
                break

    print(f"Wrote {kept} ligands to {output_path}")


if __name__ == "__main__":
    main()

