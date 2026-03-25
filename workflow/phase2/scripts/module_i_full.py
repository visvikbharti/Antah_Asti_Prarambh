#!/usr/bin/env python3
"""
Module I: Publication-quality figures for Phase 2 full-scale analysis.

6 figures:
  Fig 1: Domain architecture (CATH class, top superfamilies, domain count)
  Fig 2: N-vs-C stability (contact order violins, pLDDT violins, heatmap)
  Fig 3: GroEL class effects (CO by class, pLDDT by class)
  Fig 4: MTS targeting (classification bars, gap histogram, scatter)
  Fig 5: Orthology (overlap, conservation scatter)
  Fig 6: Summary (enriched folds, N-C asymmetry bars, MTS pie)
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

PROJECT_DIR = os.environ.get("PROJECT_DIR",
    "/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc")
RESULTS = f"{PROJECT_DIR}/results/phase2"
FIG_DIR = f"{RESULTS}/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Style
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 12, "axes.titlesize": 13,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
    "figure.dpi": 300, "savefig.dpi": 300,
})

# Colorblind-friendly palette
DATASET_COLORS = {"groel": "#4ECDC4", "hsp60": "#FF6B6B",
                  "mito_bg": "#45B7D1", "matrix_bg": "#96CEB4",
                  "proteome_bg": "#BDC3C7"}
NC_COLORS = {"N-domain": "#4ECDC4", "C-region": "#FF6B6B"}
CATH_CLASS_NAMES = {1: "Mainly Alpha", 2: "Mainly Beta",
                    3: "Alpha-Beta", 4: "Few SS"}

print("=" * 70)
print("Module I: Phase 2 Figure Generation")
print("=" * 70)


def safe_load(path, desc):
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        print(f"  {desc}: {len(df)} records")
        return df
    print(f"  [SKIP] {desc}: {path}")
    return None

def save_figure(fig, name):
    for fmt in ["pdf", "png"]:
        fig.savefig(f"{FIG_DIR}/{name}.{fmt}", bbox_inches="tight")
    plt.close(fig)
    sz = os.path.getsize(f"{FIG_DIR}/{name}.pdf")
    print(f"  Generated: {name} ({sz/1024:.0f} KB)")


# Load data
print("\nLoading data...")
domains = safe_load(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", "Domains")
dist = safe_load(f"{RESULTS}/domains/domain_distribution_full.tsv", "Distribution")
paired = safe_load(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", "N-vs-C paired")
contact_order = safe_load(f"{RESULTS}/stability/contact_order_full.tsv", "Contact order")
regions = safe_load(f"{RESULTS}/stability/region_boundaries_full.tsv", "Regions")
pvalues = safe_load(f"{RESULTS}/stats/corrected_pvalues_full.tsv", "P-values")
cath = safe_load(f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv", "CATH")
targeting = safe_load(f"{PROJECT_DIR}/results/mts/combined_targeting.tsv", "Targeting")
mts_domain = safe_load(f"{PROJECT_DIR}/results/mts/mts_domain_relationship.tsv", "MTS-domain")
homologs = safe_load(f"{PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv", "Homologs")

groel = safe_load(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", "GroEL")
hsp60 = safe_load(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", "HSP60")

# Accession sets (correct column names)
groel_acc = set()
hsp60_acc = set()
if groel is not None:
    col = "current_accession" if "current_accession" in groel.columns else "accession"
    groel_acc = set(groel[col].dropna().values)
if hsp60 is not None:
    col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
    hsp60_acc = set(hsp60[col].dropna().values)

figures_generated = 0


# ===========================================================================
# Figure 1: Domain Architecture Overview
# ===========================================================================
print("\n--- Figure 1: Domain Architecture ---")
try:
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 3, width_ratios=[1, 1.2, 0.8], wspace=0.3)

    # Panel A: CATH class distribution
    ax = fig.add_subplot(gs[0])
    if cath is not None and "cath_class" in cath.columns:
        cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
        plot_data = []
        for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
            sub = cath[cath[cath_acc_col].isin(ds_accs)]
            for cls in [1, 2, 3, 4]:
                count = (sub["cath_class"] == cls).sum()
                total = len(sub)
                plot_data.append({"Dataset": ds_name, "Class": CATH_CLASS_NAMES.get(cls, str(cls)),
                                  "Fraction": count / total if total > 0 else 0})
        if plot_data:
            pdf = pd.DataFrame(plot_data)
            pivot = pdf.pivot(index="Dataset", columns="Class", values="Fraction")
            pivot.plot(kind="bar", stacked=True, ax=ax, color=["#8B5CF6", "#D8B4FE", "#9CA3AF", "#3B82F6"])
            ax.set_ylabel("Fraction")
            ax.set_title("A. CATH Class Distribution")
            ax.legend(fontsize=8, loc="upper right")
            ax.tick_params(axis="x", rotation=0)
    ax.set_xlabel("")

    # Panel B: Top superfamilies
    ax = fig.add_subplot(gs[1])
    if cath is not None:
        cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
        sf_col = "cath_superfamily" if "cath_superfamily" in cath.columns else None
        if sf_col:
            groel_sf = cath[cath[cath_acc_col].isin(groel_acc)][sf_col].value_counts().head(8)
            hsp60_sf = cath[cath[cath_acc_col].isin(hsp60_acc)][sf_col].value_counts().head(8)
            all_sfs = list(dict.fromkeys(list(groel_sf.index) + list(hsp60_sf.index)))[:10]

            y = np.arange(len(all_sfs))
            w = 0.35
            g_counts = [groel_sf.get(sf, 0) for sf in all_sfs]
            h_counts = [hsp60_sf.get(sf, 0) for sf in all_sfs]

            ax.barh(y + w/2, g_counts, w, label="GroEL", color=DATASET_COLORS["groel"])
            ax.barh(y - w/2, h_counts, w, label="HSP60", color=DATASET_COLORS["hsp60"])
            ax.set_yticks(y)
            sf_labels = [sf[:15] for sf in all_sfs]
            ax.set_yticklabels(sf_labels, fontsize=8)
            ax.set_xlabel("Count")
            ax.set_title("B. Top CATH Superfamilies")
            ax.legend(fontsize=8)

    # Panel C: Domain count distribution
    ax = fig.add_subplot(gs[2])
    if dist is not None and len(dist) > 0:
        for ds in ["groel", "hsp60"]:
            ds_data = dist[dist["dataset"] == ds]
            if len(ds_data) > 0:
                ax.bar(ds_data["n_domains"] + (0.15 if ds == "groel" else -0.15),
                       ds_data["percent"], width=0.3,
                       label=ds.upper(), color=DATASET_COLORS[ds], alpha=0.8)
        ax.set_xlabel("Number of Domains")
        ax.set_ylabel("Percentage")
        ax.set_title("C. Domain Count")
        ax.legend(fontsize=8)
        ax.set_xlim(-0.5, 6.5)

    save_figure(fig, "fig1_domain_architecture")
    figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


# ===========================================================================
# Figure 2: N-vs-C Structural Comparison
# ===========================================================================
print("\n--- Figure 2: N-vs-C Stability ---")
try:
    if paired is not None and len(paired) > 0:
        fig = plt.figure(figsize=(16, 5))
        gs = GridSpec(1, 3, wspace=0.35)

        # Build long-form data for violins
        datasets_to_plot = ["groel", "hsp60", "matrix_bg", "mito_bg"]
        metrics_for_violin = [
            ("relative_contact_order", "Relative Contact Order"),
            ("mean_plddt", "Mean pLDDT"),
        ]

        for panel_idx, (metric, ylabel) in enumerate(metrics_for_violin):
            ax = fig.add_subplot(gs[panel_idx])
            plot_rows = []
            for ds in datasets_to_plot:
                sub = paired[paired["datasets"].str.contains(ds)]
                for _, row in sub.iterrows():
                    n_val = row.get(f"{metric}_n_domain")
                    c_val = row.get(f"{metric}_c_region")
                    if pd.notna(n_val):
                        plot_rows.append({"Dataset": ds.upper(), "Region": "N-domain", "Value": n_val})
                    if pd.notna(c_val):
                        plot_rows.append({"Dataset": ds.upper(), "Region": "C-region", "Value": c_val})

            if plot_rows:
                pdf = pd.DataFrame(plot_rows)
                sns.violinplot(data=pdf, x="Dataset", y="Value", hue="Region",
                               split=True, ax=ax, palette=NC_COLORS, cut=0, inner="quartile")
                ax.set_ylabel(ylabel)
                ax.set_title(f"{'AB'[panel_idx]}. {ylabel}")
                if panel_idx > 0:
                    ax.get_legend().remove()
                else:
                    ax.legend(fontsize=8, loc="upper right")
                ax.tick_params(axis="x", rotation=30)

        # Panel C: Heatmap of N-C differences
        ax = fig.add_subplot(gs[2])
        heatmap_metrics = ["relative_contact_order", "mean_plddt", "frac_helix",
                           "frac_strand", "frac_hydrophobic", "frac_plddt_gt70"]
        heatmap_data = []
        for ds in datasets_to_plot:
            sub = paired[paired["datasets"].str.contains(ds)]
            row_data = {}
            for m in heatmap_metrics:
                diff_col = f"{m}_diff"
                if diff_col in sub.columns:
                    row_data[m] = sub[diff_col].median()
                else:
                    row_data[m] = np.nan
            heatmap_data.append(row_data)

        if heatmap_data:
            hm_df = pd.DataFrame(heatmap_data, index=[ds.upper() for ds in datasets_to_plot])
            short_names = {"relative_contact_order": "RCO", "mean_plddt": "pLDDT",
                           "frac_helix": "Helix", "frac_strand": "Strand",
                           "frac_hydrophobic": "Hydro.", "frac_plddt_gt70": "pLDDT>70"}
            hm_df.columns = [short_names.get(c, c) for c in hm_df.columns]
            sns.heatmap(hm_df.T, ax=ax, cmap="RdBu_r", center=0, annot=True,
                        fmt=".3f", linewidths=0.5, cbar_kws={"label": "N - C"})
            ax.set_title("C. N-C Difference (median)")

        save_figure(fig, "fig2_n_vs_c_stability")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


# ===========================================================================
# Figure 3: GroEL Class-Specific Effects
# ===========================================================================
print("\n--- Figure 3: GroEL Class Effects ---")
try:
    if paired is not None and "groel_class" in paired.columns:
        groel_sub = paired[paired["datasets"].str.contains("groel")]
        groel_sub = groel_sub[groel_sub["groel_class"].isin(["I", "II", "III"])]

        if len(groel_sub) > 10:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            for idx, (metric, ylabel) in enumerate([
                ("relative_contact_order_diff", "RCO Difference (N - C)"),
                ("mean_plddt_diff", "pLDDT Difference (N - C)"),
            ]):
                ax = axes[idx]
                if metric in groel_sub.columns:
                    order = ["I", "II", "III"]
                    sns.boxplot(data=groel_sub, x="groel_class", y=metric,
                                order=order, ax=ax, palette="Set2", fliersize=2)
                    sns.stripplot(data=groel_sub, x="groel_class", y=metric,
                                  order=order, ax=ax, color="black", alpha=0.3,
                                  size=3, jitter=True)
                    ax.axhline(0, color="red", linestyle="--", alpha=0.5)
                    ax.set_xlabel("GroEL Substrate Class")
                    ax.set_ylabel(ylabel)
                    ax.set_title(f"{'AB'[idx]}. {ylabel}")

                    # Annotate sample sizes
                    for i, cls in enumerate(order):
                        n = len(groel_sub[groel_sub["groel_class"] == cls])
                        ax.text(i, ax.get_ylim()[1], f"n={n}", ha="center",
                                fontsize=8, style="italic")

            plt.tight_layout()
            save_figure(fig, "fig3_groel_class_comparison")
            figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


# ===========================================================================
# Figure 4: MTS Targeting
# ===========================================================================
print("\n--- Figure 4: MTS Targeting ---")
try:
    has_mts_data = (targeting is not None or mts_domain is not None)
    if has_mts_data:
        fig = plt.figure(figsize=(16, 5))
        gs = GridSpec(1, 3, wspace=0.35)

        # Panel A: Targeting classification
        ax = fig.add_subplot(gs[0])
        if targeting is not None:
            targ_acc_col = "uniprot_accession" if "uniprot_accession" in targeting.columns else "accession"
            targ_col = None
            for col in ["targeting_classification", "localization", "sub_localization"]:
                if col in targeting.columns:
                    targ_col = col
                    break
            if targ_col:
                hsp60_targ = targeting[targeting[targ_acc_col].isin(hsp60_acc)]
                if len(hsp60_targ) > 0:
                    counts = hsp60_targ[targ_col].value_counts().head(8)
                    colors = plt.cm.Set2(np.linspace(0, 1, len(counts)))
                    counts.plot(kind="barh", ax=ax, color=colors)
                    for i, v in enumerate(counts.values):
                        ax.text(v + 0.5, i, str(v), va="center", fontsize=9)
                    ax.set_xlabel("Count")
                    ax.set_title("A. HSP60 Substrate Targeting")

        # Panel B: MTS-to-domain gap
        ax = fig.add_subplot(gs[1])
        if mts_domain is not None and "gap_length" in mts_domain.columns:
            gaps = mts_domain["gap_length"].dropna()
            if len(gaps) > 0:
                ax.hist(gaps, bins=30, color="#45B7D1", edgecolor="white", alpha=0.8)
                ax.axvline(0, color="red", linestyle="--", linewidth=2, label="Gap = 0")
                ax.axvline(gaps.median(), color="orange", linestyle="--",
                           label=f"Median = {gaps.median():.0f}")
                ax.set_xlabel("Gap Length (residues)")
                ax.set_ylabel("Count")
                ax.set_title("B. MTS-to-Domain Gap")
                ax.legend(fontsize=9)

                pre_col = "mts_is_pre_domain" if "mts_is_pre_domain" in mts_domain.columns else None
                if pre_col:
                    n_pre = mts_domain[pre_col].sum()
                    n_total = len(mts_domain[mts_domain[pre_col].notna()])
                    ax.text(0.95, 0.95, f"Pre-domain: {n_pre}/{n_total}\n({100*n_pre/n_total:.1f}%)",
                            transform=ax.transAxes, ha="right", va="top", fontsize=9,
                            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        # Panel C: MTS cleavage vs first domain start
        ax = fig.add_subplot(gs[2])
        if mts_domain is not None:
            tp_col = "transit_peptide_end" if "transit_peptide_end" in mts_domain.columns else None
            fd_col = "first_domain_start" if "first_domain_start" in mts_domain.columns else None
            if tp_col and fd_col:
                valid = mts_domain[[tp_col, fd_col]].dropna()
                if len(valid) > 0:
                    ax.scatter(valid[tp_col], valid[fd_col], alpha=0.5,
                               s=20, color="#4ECDC4", edgecolors="none")
                    lim = max(valid[tp_col].max(), valid[fd_col].max()) * 1.1
                    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, label="y = x")
                    ax.set_xlabel("MTS Cleavage Site")
                    ax.set_ylabel("First Domain Start")
                    ax.set_title("C. MTS vs Domain Boundary")
                    ax.legend(fontsize=9)

        save_figure(fig, "fig4_mts_targeting")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


# ===========================================================================
# Figure 5: Cross-Species Orthology
# ===========================================================================
print("\n--- Figure 5: Orthology ---")
try:
    if homologs is not None and paired is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Panel A: Orthogroup overlap
        ax = axes[0]
        ortho_path = f"{PROJECT_DIR}/results/homology/substrate_orthogroups.tsv"
        if os.path.exists(ortho_path):
            ortho = pd.read_csv(ortho_path, sep="\t")
            cat_col = "category" if "category" in ortho.columns else None
            if cat_col:
                counts = ortho[cat_col].value_counts()
                counts.plot(kind="bar", ax=ax, color=["#4ECDC4", "#FF6B6B", "#96CEB4"])
                ax.set_ylabel("Count")
                ax.set_title("A. Substrate Orthogroup Categories")
                ax.tick_params(axis="x", rotation=30)
        else:
            ax.text(0.5, 0.5, "Orthogroup data\nnot available",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_title("A. Orthogroup Overlap")

        # Panel B: Conservation of N-domain contact order
        ax = axes[1]
        groel_col = "groel_accession" if "groel_accession" in homologs.columns else None
        hsp60_col = "hsp60_accession" if "hsp60_accession" in homologs.columns else None

        if groel_col and hsp60_col:
            # Build RCO lookup from paired data
            rco_lookup = {}
            if "relative_contact_order_n_domain" in paired.columns:
                for _, row in paired.iterrows():
                    rco_lookup[row["accession"]] = row["relative_contact_order_n_domain"]

            scatter_x, scatter_y = [], []
            for _, row in homologs.iterrows():
                g_acc = row[groel_col]
                h_acc = row[hsp60_col]
                if g_acc in rco_lookup and h_acc in rco_lookup:
                    g_rco = rco_lookup[g_acc]
                    h_rco = rco_lookup[h_acc]
                    if pd.notna(g_rco) and pd.notna(h_rco):
                        scatter_x.append(g_rco)
                        scatter_y.append(h_rco)

            if len(scatter_x) > 3:
                ax.scatter(scatter_x, scatter_y, alpha=0.6, s=30,
                           color="#4ECDC4", edgecolors="gray", linewidth=0.5)
                lim = max(max(scatter_x), max(scatter_y)) * 1.1
                ax.plot([0, lim], [0, lim], "k--", alpha=0.3)
                r, p = stats.pearsonr(scatter_x, scatter_y)
                ax.text(0.05, 0.95, f"r = {r:.2f}\np = {p:.1e}\nn = {len(scatter_x)}",
                        transform=ax.transAxes, va="top", fontsize=10,
                        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
                ax.set_xlabel("GroEL Homolog N-domain RCO")
                ax.set_ylabel("HSP60 Homolog N-domain RCO")
                ax.set_title("B. N-Domain Contact Order Conservation")
            else:
                ax.text(0.5, 0.5, f"Insufficient homolog pairs\nwith RCO data (n={len(scatter_x)})",
                        transform=ax.transAxes, ha="center", va="center")
                ax.set_title("B. RCO Conservation")

        plt.tight_layout()
        save_figure(fig, "fig5_orthology")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


# ===========================================================================
# Figure 6: Key Findings Summary
# ===========================================================================
print("\n--- Figure 6: Summary ---")
try:
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 3, wspace=0.35)

    # Panel A: Mean N-C RCO difference by dataset
    ax = fig.add_subplot(gs[0])
    if paired is not None and "relative_contact_order_diff" in paired.columns:
        ds_labels = []
        ds_means = []
        ds_sems = []
        for ds in ["groel", "hsp60", "matrix_bg", "mito_bg"]:
            sub = paired[paired["datasets"].str.contains(ds)]
            diffs = sub["relative_contact_order_diff"].dropna()
            if len(diffs) >= 5:
                ds_labels.append(ds.upper())
                ds_means.append(diffs.mean())
                ds_sems.append(diffs.sem())

        if ds_labels:
            colors = [DATASET_COLORS.get(dl.lower(), "#999") for dl in ds_labels]
            ax.bar(range(len(ds_labels)), ds_means, yerr=ds_sems,
                   color=colors, capsize=5, edgecolor="gray")
            ax.axhline(0, color="black", linestyle="-", linewidth=0.5)
            ax.set_xticks(range(len(ds_labels)))
            ax.set_xticklabels(ds_labels, rotation=30)
            ax.set_ylabel("Mean RCO Diff (N - C)")
            ax.set_title("A. N-C Contact Order Asymmetry")

    # Panel B: Statistical significance overview
    ax = fig.add_subplot(gs[1])
    if pvalues is not None and len(pvalues) > 0:
        sig_col = "significant_overall" if "significant_overall" in pvalues.columns else "significant"
        if sig_col in pvalues.columns:
            sig_counts = pvalues.groupby("family")[sig_col].agg(["sum", "count"])
            sig_counts.columns = ["Significant", "Total"]
            sig_counts["Not Significant"] = sig_counts["Total"] - sig_counts["Significant"]

            sig_counts[["Significant", "Not Significant"]].plot(
                kind="barh", stacked=True, ax=ax,
                color=["#E74C3C", "#BDC3C7"])
            ax.set_xlabel("Number of Tests")
            ax.set_title("B. Statistical Significance")
            ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "P-values\nnot available",
                transform=ax.transAxes, ha="center", va="center")
        ax.set_title("B. Statistical Tests")

    # Panel C: MTS-domain relationship pie
    ax = fig.add_subplot(gs[2])
    if mts_domain is not None:
        pre_col = "mts_is_pre_domain" if "mts_is_pre_domain" in mts_domain.columns else None
        if pre_col:
            valid = mts_domain[mts_domain[pre_col].notna()]
            n_pre = valid[pre_col].sum()
            n_overlap = len(valid) - n_pre
            ax.pie([n_pre, n_overlap],
                   labels=[f"Pre-domain\n(n={n_pre})", f"Overlaps domain\n(n={n_overlap})"],
                   colors=["#96CEB4", "#FF6B6B"],
                   autopct="%1.1f%%", startangle=90)
            ax.set_title("C. MTS-Domain Relationship")
    else:
        ax.text(0.5, 0.5, "MTS data\nnot available",
                transform=ax.transAxes, ha="center", va="center")
        ax.set_title("C. MTS-Domain Relationship")

    save_figure(fig, "fig6_summary")
    figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")


print(f"\n{'=' * 70}")
print(f"Module I complete. Generated {figures_generated}/6 figures in {FIG_DIR}")
print(f"{'=' * 70}")
