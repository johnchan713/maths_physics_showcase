#!/usr/bin/env python3
"""Plot recorded observations only, with no smoothing or extrapolation."""
import json

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "late_growth_t014"
import matplotlib.pyplot as plt

from analyze import ROOT, history


def main():
    report = json.loads((ROOT / "assessment.json").read_text())
    if report["status"] == "analysis-in-progress":
        raise ValueError("Wait for the completed paired and sampling assessment")
    protocol = json.loads((ROOT / "protocol.json").read_text())
    data = {grid: history(grid) for grid in protocol["evolution_grids"]}
    fig, panels = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    axes = panels.ravel()
    titles = [r"Critical $\dot H^{1/2}$ norm", r"Critical $L^3$ norm",
              r"Relative $\dot H^{1/2}$ growth rate", "Sampled maximum vorticity"]
    for axis, title in zip(axes, titles):
        axis.set_title(title, fontsize=12)
        axis.set_xlabel("Simulation time")
        axis.axvspan(protocol["start_time"], protocol["endpoints"][-1],
                     color="#e7e7e7", zorder=0)
        axis.grid(alpha=.22)
        axis.spines[["top", "right"]].set_visible(False)
    for grid, color, style, marker in ((64, "#b85012", "--", "o"),
                                      (128, "#1267a4", "-", ".")):
        joined, rows, _ = data[grid]
        time = [row["time"] for row in rows]
        options = dict(color=color, linestyle=style, marker=marker,
                       markersize=4, linewidth=1.6, label=f"{grid}³ evolution")
        axes[0].plot(time, [row["h_half_ratio"] for row in rows], **options)
        axes[1].plot(time, [row["l3_dense_ratio"] for row in rows], **options)
        axes[2].plot([row["time"] for row in joined["budgets"]],
                     [row["h_half_logarithmic_rate"] for row in joined["budgets"]],
                     **options)
        axes[3].plot(time, [row["vorticity_dense_ratio"] for row in rows], **options)
    for axis in axes[:2]:
        axis.set_ylabel("Norm / original initial norm")
        axis.axhline(1, color="#666666", linewidth=.7)
    axes[2].set_ylabel(r"$d\log\|u\|_{\dot H^{1/2}}/dt$")
    axes[2].axhline(0, color="#666666", linewidth=.7)
    axes[2].set_ylim(bottom=0)
    axes[3].set_ylabel("Maximum / original initial maximum")
    axes[3].text(.02, .96, "Both fields sampled on a 128³ physical grid",
                 transform=axes[3].transAxes, va="top", fontsize=9)
    final = str(protocol["endpoints"][-1])
    gap = report["pairs"][final]["maximum_trajectory_relative_gaps"]["vorticity_dense_ratio"]
    axes[3].text(.98, .08, f"Largest gap: {100 * gap:.2f}%\nDeclared limit: 10%",
                 transform=axes[3].transAxes, ha="right", fontsize=10)
    axes[0].legend(frameon=False, fontsize=9, loc="upper left")
    fig.suptitle("Frozen candidate: observed finite growth, no extrapolation\n"
                 "Shading marks the new .12–.14 interval; lines connect saved samples.", fontsize=12)
    fig.savefig(ROOT / "growth.svg", metadata={"Date": None})
    vector = ROOT / "growth.svg"
    vector.write_text("\n".join(line.rstrip() for line in vector.read_text().splitlines()) + "\n")
    fig.savefig(ROOT / "growth.png", dpi=170)
    plt.close(fig)
    print("Saved growth.svg and growth.png from the verified joined histories")


if __name__ == "__main__":
    main()
