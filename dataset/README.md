## Evo 2 sparse autoencoder activations for the *Escherichia coli* K-12 MG1655 genome

SAE activations for every base of the chromosome (RefSeq GCF_000005845.2, NC_000913.3,
4,641,652 bp), with the embeddings, annotation and structures used by the
[hands-on notebook](https://github.com/GenAIBio/evo2-sae-handson).

### How the activations were made

Evo 2 7B (`arcinstitute/evo2_7b`) was run over the circular chromosome in 16,384 bp
windows at a stride of 15,360 bp, discarding the first 1,024 bp of each window as context.
Each base's `blocks.26` embedding, 4,096 values, was encoded with the Goodfire BatchTopK
SAE (`Goodfire/Evo-2-Layer-26-Mixed`, `sae-layer26-mixed-expansion_8-k_64.pt`, sha256
prefix `598f131a2e71f1c4`) the way the
[Evo 2 notebook](https://github.com/arcinstitute/evo2/blob/main/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb)
does it:

    a = relu(x @ W + b_enc)
    z = BatchTopK_64(a)

BatchTopK keeps the largest 64 x L activations over the L bases it is given, so the result
depends on the window. Rows here hold the largest 128 activations per base, twice what
BatchTopK keeps on average, so the cut-off can be reapplied for any window. A base that
exceeds the cut-off more than 128 times loses the surplus — 0.16% of the kept activations
over the 30,720 bp of `region_emb.pt`, 0.57% over a 5 kb window, always the smallest values
just above the cut-off.

### Files

| File | Bytes | Description |
|---|---:|---|
| `ecoli_values.f16` | 1,188,262,912 | the largest 128 activations per base, float16 |
| `ecoli_indices.u16` | 1,188,262,912 | their feature ids, uint16, 0–32767 |
| `ecoli_indptr.i64` | 37,133,224 | row pointers, int64, 4,641,653 entries, 128 apart |
| `ecoli_meta.json` | 1,000 | provenance of the three files above |
| `region_emb.pt` | 251,659,966 | `blocks.26` embeddings for NC_000913.3:4,162,560-4,193,280; 30,720 × 4,096 bfloat16 |
| `context.pt` | 83,888,037 | five window overlaps, each 1,024 bases embedded twice, with 0–1,023 and with ≥ 15,360 bp of left context |
| `scramble.pt` | 183,506,633 | the CRISPR array embedded sixteen ways: natural, and five draws of each of three edits; 16 × 1,400 × 4,096 bfloat16 |
| `class_sums.npz` | 656,598 | per-class totals of the per-region mean activation, over 7,279 regions and 32,768 features, with the genome-wide BatchTopK cut-off applied |
| `ecoli.bed` | 182,381 | 4,340 CDS, 86 tRNA, 22 rRNA, 20 CRISPR repeats, 18 spacers, 6 prophages |
| `proteins.json` | 5,468 | CDS coordinates, UniProt ids and DSSP strings for seven proteins; residue *i* is codon *i* |
| `eftu_trna.pdb` | 379,071 | EF-Tu (AlphaFold P0CE48) on 1B23, keeping the 1B23 tRNA as chain B |
| `rpoBC.pdb` | 1,701,413 | RNA polymerase β and β′ chains from 6C9Y |
| `pspA.pdb` | 149,849 | AlphaFold P0AFM6 |
| `ompG.pdb` | 205,982 | AlphaFold P76045 |
| `crp.pdb` | 139,643 | AlphaFold P0ACJ8 |
| `groL.pdb` | 331,451 | AlphaFold P0A6F5 |
| `sha256sums.txt` | — | checksums for all of the above |

### Reading one region

Rows `lo:hi` occupy `[2*indptr[lo], 2*indptr[hi])` of the values and indices files, and the
host serves HTTP Range requests, so a few-kb region costs a MB or two. Apply BatchTopK over
the window before reading a feature.

```python
import numpy as np, requests

BASE = "https://huggingface.co/datasets/suzuki-2001/evo2-sae-handson/resolve/main"

def part(name, first, last):
    r = requests.get(f"{BASE}/{name}",
                     headers={"Range": f"bytes={first}-{last - 1}"})
    r.raise_for_status()
    return r.content

lo, hi = 4_175_250, 4_177_250  # tRNA array and tufB
ptr = np.frombuffer(part("ecoli_indptr.i64", 8 * lo, 8 * (hi + 1)), np.int64)
first, last = int(ptr[0]), int(ptr[-1])
values = np.frombuffer(part("ecoli_values.f16", 2 * first, 2 * last), np.float16)
indices = np.frombuffer(part("ecoli_indices.u16", 2 * first, 2 * last), np.uint16)
rows = np.repeat(np.arange(hi - lo), np.diff(ptr))

value = values.astype(np.float32)                       # BatchTopK over this window
want = 64 * (hi - lo)
floor = np.partition(value, len(value) - want)[len(value) - want]
keep = value >= floor

hit = keep & (indices == 30262)  # the paper's tRNA feature
track = np.zeros(hi - lo, np.float32)
track[rows[hit]] = value[hit]
```

### Sources and licences

- Genome NC_000913.3 — NCBI imposes no restrictions on use or redistribution; rights
  remain with the original submitters
- Evo 2 7B weights, `arcinstitute/evo2_7b` — Apache 2.0
- SAE checkpoint, `Goodfire/Evo-2-Layer-26-Mixed` — MIT
- AlphaFold models P0CE48, P0A8V2, P0A8T7, P0AFM6, P76045, P0ACJ8, P0A6F5 — CC-BY 4.0,
  EMBL-EBI
- Experimental structures 1B23 and 6C9Y — CC0 1.0, RCSB PDB
- DSSP strings — mkdssp 4.6.1 (BSD-2-Clause), run on the AlphaFold models above
- Prophage intervals — geNomad (Lawrence Berkeley National Laboratory academic /
  non-commercial licence); CPZ-55, doi:10.3389/fgene.2019.00065
- CRISPR repeats and spacers — derived from the genome sequence

Released under CC-BY 4.0.

### Citation

Brixi G, Durrant MG, Ku J, Naghipourfar M, Poli M, Sun G, et al. Genome modelling and
design across all domains of life with Evo 2. *Nature* 2026.
doi:10.1038/s41586-026-10176-5

Archived at https://doi.org/10.5281/zenodo.21856795
