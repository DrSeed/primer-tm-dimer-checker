# Self-contained demo: score 300 random primers and compare Tm methods.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from primer_check import check_primer, check_pair


def random_primers(n=300, seed=1):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        length = int(rng.integers(16, 29))
        gc = rng.uniform(0.3, 0.7)
        probs = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        out.append("".join(rng.choice(list("AGCT"), size=length, p=probs)))
    return out


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    df = pd.DataFrame([check_primer(p) for p in random_primers()])
    df.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    sc = ax.scatter(df["tm_wallace"], df["tm_nn"], c=df["length"], cmap="viridis", s=18)
    lims = [df[["tm_wallace", "tm_nn"]].min().min() - 2, df[["tm_wallace", "tm_nn"]].max().max() + 2]
    ax.plot(lims, lims, "k--", lw=1, label="perfect agreement")
    ax.set_xlabel("Wallace rule Tm (C)")
    ax.set_ylabel("Nearest-neighbour Tm (C)")
    ax.set_title("Rule of thumb vs thermodynamics")
    ax.legend()
    fig.colorbar(sc, ax=ax, label="primer length")

    ax = axes[1]
    err = df["tm_wallace"] - df["tm_nn"]
    ax.scatter(df["length"] + np.random.default_rng(0).uniform(-0.2, 0.2, len(df)), err, s=14, alpha=0.7)
    ax.axhline(0, color="k", lw=1)
    ax.set_xlabel("Primer length (nt)")
    ax.set_ylabel("Wallace minus NN (C)")
    ax.set_title("The error grows with length")

    ax = axes[2]
    flags = df["flags"].str.split(";").explode()
    counts = flags[flags != ""].value_counts()
    ax.barh(counts.index, counts.values, color="#c0504d")
    ax.set_title(f"Why primers fail ({(~df['pass']).sum()} of {len(df)} flagged)")
    ax.set_xlabel("Primers flagged")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)

    good = df[df["pass"]]["primer"].tolist()
    if len(good) >= 2:
        print("Example pair check:", check_pair(good[0], good[1]))
    print(df[["length", "gc", "tm_wallace", "tm_nn"]].describe().round(1))


if __name__ == "__main__":
    main()
