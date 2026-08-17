"""ゲノムブラウザのトラックと, 領域オブジェクト."""
import logging

import numpy as np
import pandas as pd
from coolbox.api import BED, Spacer, XAxis
from coolbox.core.track.hist.base import HistBase

from . import data

# bgzip/tabix が無いと coolbox は毎回警告を出しつつメモリ上で読む. この規模なら問題ない
for _name in list(logging.root.manager.loggerDict):
    if _name.startswith("coolbox"):
        logging.getLogger(_name).setLevel(logging.CRITICAL)

ANNO_STYLE = {
    "ORF": "#64748b",
    "tRNA": "#7c3aed",
    "rRNA": "#059669",
    "repeat": "#334155",
    "spacer": "#db2777",
    "prophage": "#b91c1c",
}

_annotation = None


def annotation():
    """ecoli.bed の全区間. 初回に取得する."""
    global _annotation
    if _annotation is None:
        _annotation = []
        for line in open(data.download("ecoli.bed")):
            field = line.rstrip("\n").split("\t")
            kind, name = field[3].split(":", 1)
            if kind == "CDS":  # RefSeq の CDS を論文にならって ORF と呼ぶ
                kind = "ORF"
            _annotation.append({
                "start": int(field[1]),
                "end": int(field[2]),
                "name": name,
                "kind": kind,
                "strand": field[5],
            })
    return _annotation


def overlapping(lo, hi, kind=None):
    return [r for r in annotation()
            if r["end"] > lo and r["start"] < hi and (kind is None or r["kind"] == kind)]


class Signal(HistBase):
    """1 次元配列をそのまま表示するトラック.

    coolbox は表示範囲が変わるたびトラックに問い合わせるので, 配列を保持して
    範囲ぶんを返す形にする.
    """

    def __init__(self, values, chrom, start, **kwargs):
        super().__init__(style=HistBase.STYLE_FILL, **kwargs)
        self.values, self.chrom, self.start = np.asarray(values, np.float32), chrom, start

    def fetch_data(self, gr, **kwargs):
        lo = max(gr.start, self.start)
        hi = min(gr.end, self.start + len(self.values))
        return pd.DataFrame({
            "pos": np.arange(lo, hi) + 0.5,
            "score": self.values[lo - self.start:hi - self.start],
        })

    def fetch_plot_data(self, gr, **kwargs):
        return self.fetch_data(gr, **kwargs)


class Axis(XAxis):
    """coolbox 既定の目盛は kb 丸めなので, 狭い範囲だと同じ数字が並ぶ."""

    def plot(self, ax, gr, **kwargs):
        self.ax = ax
        ax.set_xlim(gr.start, gr.end)
        ticks, span = ax.get_xticks(), gr.end - gr.start
        bp = span <= 20_000
        labels = [f"{t:,.0f}" if bp else f"{t / 1e3:,.1f}" for t in ticks]
        labels[-2] += " bp" if bp else " kb"
        ax.axis["x"] = ax.new_floating_axis(0, 0.5)
        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis["x"].major_ticklabels.set(size=int(self.properties["fontsize"]))


def write_bed(hits, chrom, lo, hi, kind):
    """表示範囲で切り詰めた区間を書き出し, そのパスを返す.

    bgzip が無いと coolbox は BED を pandas でそのまま読む. 1 行目がヘッダとして
    消費されるのでダミーを先頭に置き, 範囲に完全に含まれる区間しか返らないので
    表示範囲ごとに切り詰めて書き直す.
    """
    path = f"{data.TRACKS}/view_{kind}.bed"
    with open(path, "w") as fh:
        fh.write(f"{chrom}\t0\t1\theader\t0\t+\n")
        for r in hits:
            fh.write(
                f"{chrom}\t{max(r['start'], lo)}\t{min(r['end'], hi)}"
                f"\t{r['name']}\t0\t{r['strand']}\n"
            )
    return path


def browse(tracks, chrom, lo, hi):
    """数値トラックを上に, アノテーションを下に積む."""
    frame = Axis()
    for track in tracks:
        frame += track
    frame += Spacer(height=0.8)
    for kind, colour in ANNO_STYLE.items():
        hits = overlapping(lo, hi, kind)
        if not hits:  # 範囲に入らない種類は表示しない
            continue
        # 一番短い区間が表示幅の 1% 未満なら, 区間名は表示しない
        named = min(r["end"] - r["start"] for r in hits) > (hi - lo) * 0.01
        frame += BED(
            write_bed(hits, chrom, lo, hi, kind),
            color=colour,
            title=kind,
            display="stacked",
            labels="auto" if named else False,
            row_height=0.42,
        )
        frame += Spacer(height=0.6)
    return frame.plot(f"{chrom}:{lo}-{hi}")


def intervals(tracks, chrom, lo, hi):
    """区間の並びだけを表示する. tracks は (トラック名, 色, [(開始, 終了, ラベル), ...]) の並び."""
    frame = Axis()
    for name, colour, spans in tracks:
        path = f"{data.TRACKS}/diagram_{name}.bed"
        with open(path, "w") as fh:
            fh.write(f"{chrom}\t0\t1\theader\t0\t+\n")
            for start, end, label in spans:
                fh.write(f"{chrom}\t{start}\t{end}\t{label}\t0\t+\n")
        frame += BED(
            path,
            color=colour,
            title=name,
            display="stacked",
            labels="auto",
            row_height=0.5,
        )
        frame += Spacer(height=0.4)
    return frame.plot(f"{chrom}:{lo}-{hi}")
