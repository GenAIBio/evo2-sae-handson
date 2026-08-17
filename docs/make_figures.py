"""ノートブックが画像として貼る図を描画する."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))


def context_windows(path):
    """縮尺は合っておらず, 重なりは見えるように広げてある."""
    width, shift = 10.0, 6.0
    base = shift + 1.5

    fig, ax = plt.subplots(figsize=(8.6, 2.5))
    rows = [
        (1.6, 0.0, "window 1", "rich context"),
        (0.6, shift, "window 2", "shallow context"),
    ]
    for row, begin, name, label in rows:
        ax.add_patch(Rectangle((begin, row - 0.16), width, 0.32, color="#334155"))
        ax.text(begin - 0.35, row, name, ha="right", va="center", fontsize=10)
        ax.add_patch(FancyArrowPatch(
            (begin, row - 0.33),
            (base, row - 0.33),
            arrowstyle="->",
            color="#dc2626",
            lw=1.3,
            mutation_scale=12,
        ))
        ax.text((begin + base) / 2, row - 0.53, label, ha="center", fontsize=10, color="#dc2626")

    ax.plot([base, base], [0.26, 1.9], color="#dc2626", lw=1.3)
    ax.text(base, 2.02, "same position", ha="center", fontsize=10, color="#dc2626")
    ax.set_xlim(-3.0, shift + width + 0.4)
    ax.set_ylim(0.0, 2.25)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150, transparent=True)
    print(f"{path}  {os.path.getsize(path) / 1e3:.0f} KB")


if __name__ == "__main__":
    context_windows(f"{HERE}/context_windows.png")
