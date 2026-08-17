"""1 つのゲノム領域と, その上の列 (埋め込みの次元, または SAE の特徴)."""
import json

import numpy as np

from . import data, genome, structure

K = 64  # チェックポイントの k. ファイル名の k_64 にあたる
LENGTH = 4_641_652  # NC_000913.3 の全長
N_FEATURE = 32_768

_proteins = None


def proteins():
    """proteins.json. ORF の座標, UniProt id, DSSP 文字列."""
    global _proteins
    if _proteins is None:
        _proteins = json.load(open(data.download("proteins.json")))
    return _proteins


def batch_topk_floor(value, k, n_rows):
    """n_rows x k 個を残すときの足切り値. 配布行列は 1 塩基あたり 128 個を持つ."""
    want = k * n_rows
    if want >= len(value):
        return 0.0
    return float(np.partition(value, len(value) - want)[len(value) - want])


class Region:
    """from_matrix は手元の行列から, from_genome は全ゲノムの活性行列から作る."""

    def __init__(self, chrom, start, end, column):
        self.chrom, self.start, self.end = chrom, start, end
        self.column = column
        self._elsewhere = {}

    def __len__(self):
        return self.end - self.start

    def track(self, index):
        """列 1 本を, この領域ぶんの 1 次元配列で返す."""
        return self.column(index)

    def show(self, columns, lo=None, hi=None, floor=True, top=None):
        """columns は (列番号, 色, 表示名) の並び. top は全トラック共通の縦軸上限."""
        tracks = [
            genome.Signal(
                self.track(i),
                self.chrom,
                self.start,
                color=colour,
                title=label,
                min_value=0 if floor else "auto",
                max_value="auto" if top is None else top,
            )
            for i, colour, label in columns
        ]
        return genome.browse(
            tracks,
            self.chrom,
            self.start if lo is None else lo,
            self.end if hi is None else hi,
        )

    def residue_acts(self, gene, feature):
        """ORF を 3 塩基ずつに区切り, コドンごとの平均を残基の活性とする."""
        p = proteins()[gene]
        inside = self.start <= p["start"] and p["end"] <= self.end
        source = self if inside else self._outside(gene, p)
        offset = p["start"] - source.start
        codons = source.track(feature)[offset:offset + 3 * p["n_res"]]
        return codons.reshape(p["n_res"], 3).mean(axis=1)

    def show_complex(self, sample, helix, sheet):
        """立体構造に活性を載せる. helix と sheet は二次構造に対応する特徴番号."""
        return structure.show_complex(
            sample,
            lambda gene, f: self.residue_acts(gene, f),
            helix,
            sheet,
        )

    def _outside(self, gene, p):
        # 配布した領域の外にある遺伝子は, 全ゲノムの活性行列から引く
        if gene not in self._elsewhere:
            self._elsewhere[gene] = Region.from_genome(self.chrom, p["start"], p["end"])
        return self._elsewhere[gene]

    @classmethod
    def from_matrix(cls, chrom, start, matrix):
        width = matrix.shape[1]

        def column(index):
            if not 0 <= index < width:
                raise ValueError(f"列は 0 から {width - 1} です. 受け取ったのは {index}")
            return np.asarray(matrix[:, index])

        return cls(chrom, start, start + len(matrix), column)

    @classmethod
    def from_genome(cls, chrom, lo, hi, k=K):
        """取り出した範囲を 1 つの入力として BatchTopK をかける. Evo 2 の公開実装と同じ."""
        lo, hi = int(lo), int(hi)
        if not 0 <= lo < hi <= LENGTH:
            raise ValueError(
                f"座標は 0 から {LENGTH:,} で, 開始 < 終了 です. "
                f"受け取ったのは {lo:,}-{hi:,}"
            )
        ptr = np.frombuffer(data.fetch("ecoli_indptr.i64", 8 * lo, 8 * (hi + 1)), np.int64)
        first, last = int(ptr[0]), int(ptr[-1])
        values = np.frombuffer(data.fetch("ecoli_values.f16", 2 * first, 2 * last), np.float16)
        indices = np.frombuffer(data.fetch("ecoli_indices.u16", 2 * first, 2 * last), np.uint16)
        rows = np.repeat(np.arange(hi - lo), np.diff(ptr))
        value = values.astype(np.float32)
        floor = batch_topk_floor(value, k, hi - lo)
        keep = value >= floor
        rows, indices, value = rows[keep], indices[keep], value[keep]

        def column(index):
            if not 0 <= index < N_FEATURE:
                raise ValueError(
                    f"特徴のインデックスは 0 から {N_FEATURE - 1} です. 受け取ったのは {index}"
                )
            out = np.zeros(hi - lo, np.float32)
            hit = indices == index
            out[rows[hit]] = value[hit]
            return out

        return cls(chrom, lo, hi, column)
