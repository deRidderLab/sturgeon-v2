"""
Plots for ``sturgeon-v2 single`` and ``sturgeon-v2 multi``.

One plot per sample: a vertical bar for every class, in the order of the
model's classification_system.yaml, with the score averaged over all
samplings. Class colors come from the ``colors: type:`` section of that file;
classes without a color are drawn in grey. Plots are saved as .png.
"""

from typing import Dict, List

import numpy as np

DEFAULT_COLOR = "#999999"


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def get_class_colors(classification_file: dict) -> Dict[str, str]:
    """Return {class name: hex color} from a loaded classification_system.yaml."""
    colors = (classification_file.get("colors") or {}).get("type") or {}
    return {str(k): str(v).strip() for k, v in colors.items()}


def plot_scores(
    sample_name: str,
    class_names: List[str],
    scores: np.ndarray,
    colors: Dict[str, str],
    n_cpgs: int,
    output_file: str,
) -> None:
    """Vertical bar chart of the score of every class for one sample.

    ``class_names`` is in the order of the classification system yaml and
    ``scores`` holds the score per class, averaged over all samplings.
    """
    plt = _pyplot()

    scores = np.asarray(scores, dtype=float)
    best = int(np.argmax(scores))
    x = np.arange(len(class_names))

    fig, ax = plt.subplots(figsize=(max(8.0, 0.17 * len(class_names) + 2), 5))
    ax.bar(
        x, scores,
        color=[colors.get(c, DEFAULT_COLOR) for c in class_names],
        edgecolor="black", linewidth=0.3, width=0.8,
    )
    ax.text(
        best, scores[best] + 0.015, f"{scores[best]:.3f}",
        ha="center", va="bottom", fontsize=8,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(class_names, rotation=90, fontsize=6)
    ax.get_xticklabels()[best].set_fontweight("bold")
    ax.set_xlim(-0.7, len(class_names) - 0.3)
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("Score")
    ax.set_title(
        f"{sample_name}: {class_names[best]} ({scores[best]:.3f}), {n_cpgs} CpG sites",
        fontsize=10,
    )
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_file, dpi=200, format="png")
    plt.close(fig)
