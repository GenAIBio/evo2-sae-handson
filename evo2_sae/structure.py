"""立体構造に残基ごとの活性を載せる."""
import numpy as np
import py3Dmol

from . import data

HELIX_COLOUR = "#2563eb"
SHEET_COLOUR = "#f59e0b"
OTHER_COLOUR = "#d4d4d8"
RNA_COLOUR = "#7a4fa3"

# 表示名: (PDB, {チェーン: 遺伝子}, RNA のチェーン)
SAMPLES = {
    "tufB": ("eftu_trna.pdb", {"A": "tufB"}, "B"),
    "rpoB-rpoC": ("rpoBC.pdb", {"A": "rpoB", "B": "rpoC"}, None),
    "pspA": ("pspA.pdb", {"A": "pspA"}, None),
    "ompG": ("ompG.pdb", {"A": "ompG"}, None),
    "crp": ("crp.pdb", {"A": "crp"}, None),
    "groL": ("groL.pdb", {"A": "groL"}, None),
}


def set_bfactors(pdb_path, values):
    """values は {(チェーン, 残基番号): 値}. PDB の B-factor 欄 (61-66 桁) を差し替える."""
    out = []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            value = values.get((line[21], int(line[22:26])))
            if value is not None:
                line = f"{line[:60]}{min(value, 999):6.2f}{line[66:]}"
        out.append(line)
    return "".join(out)


def classify(acts, gene, helix_feature, sheet_feature):
    """発火した残基を, 値の大きいほうの特徴に割り当てる. 番号は 1 始まり."""
    helix, sheet = acts(gene, helix_feature), acts(gene, sheet_feature)
    is_helix = (helix > 0) & (helix >= sheet)
    is_sheet = (sheet > 0) & (sheet > helix)
    return (
        (np.flatnonzero(is_helix) + 1).tolist(),
        (np.flatnonzero(is_sheet) + 1).tolist(),
        np.maximum(helix, sheet),
    )


def show_complex(sample, acts, helix_feature, sheet_feature):
    """acts(gene, feature) が返す残基ごとの活性で塗り分ける."""
    pdb_name, chains, rna_chain = SAMPLES[sample]
    bfactors, styles = {}, []
    for chain, gene in chains.items():
        helix, sheet, strength = classify(acts, gene, helix_feature, sheet_feature)
        styles += [(chain, helix, HELIX_COLOUR), (chain, sheet, SHEET_COLOUR)]
        for i, value in enumerate(strength):
            bfactors[(chain, i + 1)] = float(value)

    view = py3Dmol.view(width=760, height=460)
    view.addModel(set_bfactors(data.download(pdb_name), bfactors), "pdb")
    view.setStyle({"cartoon": {"color": OTHER_COLOUR}})
    for chain, residues, colour in styles:
        view.addStyle({"chain": chain, "resi": residues}, {"cartoon": {"color": colour}})
    if rna_chain:
        view.addStyle({"chain": rna_chain}, {"cartoon": {"color": RNA_COLOUR}})
    view.addSurface("VDW", {"opacity": 0.35, "color": "#e5e7eb"}, {"chain": list(chains)})
    view.zoomTo()
    return view.show()
