#!/usr/bin/env python3
"""
Module I (Polished): Publication-quality figures for Phase 2 full-scale analysis.

ALL DATA IS READ FROM ACTUAL RESULT FILES — NO FABRICATED DATA.
Every number, label, and annotation comes from real Phase 2 outputs.

6 figures:
  Fig 1: Domain architecture (CATH class, top superfamilies, domain count)
  Fig 2: N-vs-C stability (contact order violins, pLDDT violins, heatmap)
  Fig 3: GroEL class effects (CO by class, pLDDT by class)
  Fig 4: MTS targeting (classification bars, gap histogram, scatter)
  Fig 5: Orthology (overlap, conservation scatter)
  Fig 6: Summary (N-C asymmetry bars, significance overview, MTS pie)
"""

import os
import sys
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

# ===========================================================================
# Configuration — LOCAL MAC PATHS (real data only)
# ===========================================================================
PROJECT_DIR = os.path.expanduser("~/Downloads/Antah_Asti_Prarambh")
RESULTS = f"{PROJECT_DIR}/results/phase2"
FIG_DIR = f"{RESULTS}/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.framealpha": 0.8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

# Colorblind-friendly palette (Wong 2011, Nature Methods)
DATASET_COLORS = {
    "groel": "#0072B2",      # Blue
    "hsp60": "#D55E00",      # Vermillion
    "matrix_bg": "#009E73",  # Bluish green
    "mito_bg": "#56B4E9",    # Sky blue
    "proteome_bg": "#999999" # Gray
}
NC_COLORS = {"N-domain": "#0072B2", "C-region": "#D55E00"}
CATH_CLASS_COLORS = ["#CC79A7", "#56B4E9", "#009E73", "#999999"]
CATH_CLASS_NAMES = {1: "Mainly Alpha", 2: "Mainly Beta",
                    3: "Alpha-Beta", 4: "Few SS"}

# Human-readable dataset names
DATASET_LABELS = {
    "groel": "GroEL",
    "hsp60": "HSP60",
    "matrix_bg": "Matrix",
    "mito_bg": "Mito",
    "proteome_bg": "Proteome"
}

print("=" * 70)
print("Module I (Polished): Publication Figure Generation")
print("=" * 70)
print(f"Data source: {RESULTS}")
print(f"Output: {FIG_DIR}")


def safe_load(path, desc):
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        print(f"  {desc}: {len(df)} records")
        return df
    print(f"  [MISSING] {desc}: {path}")
    return None


def save_figure(fig, name):
    for fmt in ["pdf", "png"]:
        fig.savefig(f"{FIG_DIR}/{name}.{fmt}", bbox_inches="tight")
    plt.close(fig)
    sz = os.path.getsize(f"{FIG_DIR}/{name}.pdf")
    print(f"  Saved: {name} ({sz/1024:.0f} KB)")


def add_panel_label(ax, label, x=-0.12, y=1.08):
    """Add bold panel label (A, B, C...) to axis."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="top")


# ===========================================================================
# Load ALL data from real result files
# ===========================================================================
print("\nLoading data from actual result files...")
domains = safe_load(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", "Unified domains")
dist = safe_load(f"{RESULTS}/domains/domain_distribution_full.tsv", "Domain distribution")
paired = safe_load(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", "N-vs-C paired")
contact_order = safe_load(f"{RESULTS}/stability/contact_order_full.tsv", "Contact order")
regions = safe_load(f"{RESULTS}/stability/region_boundaries_full.tsv", "Region boundaries")
pvalues = safe_load(f"{RESULTS}/stats/corrected_pvalues_full.tsv", "Corrected p-values")
cath = safe_load(f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv", "CATH domains")
targeting = safe_load(f"{PROJECT_DIR}/results/mts/combined_targeting.tsv", "MTS targeting")
mts_domain = safe_load(f"{PROJECT_DIR}/results/mts/mts_domain_relationship.tsv", "MTS-domain")
homologs = safe_load(f"{PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv", "Homolog pairs")
groel = safe_load(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", "GroEL substrates")
hsp60 = safe_load(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", "HSP60 substrates")
ortho = safe_load(f"{PROJECT_DIR}/results/homology/substrate_orthogroups.tsv", "Substrate orthogroups")

# Build accession sets using correct column names
groel_acc = set()
hsp60_acc = set()
if groel is not None:
    col = "current_accession" if "current_accession" in groel.columns else "accession"
    groel_acc = set(groel[col].dropna().values)
    print(f"  GroEL accessions: {len(groel_acc)}")
if hsp60 is not None:
    col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
    hsp60_acc = set(hsp60[col].dropna().values)
    print(f"  HSP60 accessions: {len(hsp60_acc)}")

# Build p-value lookup from real statistics
pval_lookup = {}
if pvalues is not None:
    for _, row in pvalues.iterrows():
        pval_lookup[row["hypothesis"]] = {
            "p": row["p_value"],
            "p_bh": row["p_bh_within"],
            "sig": row.get("significant_overall", False),
            "effect": row.get("effect_size", np.nan),
            "direction": row.get("direction", ""),
        }

figures_generated = 0


# ===========================================================================
# Figure 1: Domain Architecture Overview
# ===========================================================================
print("\n--- Figure 1: Domain Architecture ---")
try:
    fig = plt.figure(figsize=(16, 5.5))
    gs = GridSpec(1, 3, width_ratios=[1, 1.3, 0.9], wspace=0.35)

    # Panel A: CATH class distribution (stacked bars)
    ax = fig.add_subplot(gs[0])
    add_panel_label(ax, "A")
    if cath is not None and "cath_class" in cath.columns:
        cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
        datasets_for_cath = [("GroEL", groel_acc), ("HSP60", hsp60_acc)]
        x_pos = np.arange(len(datasets_for_cath))
        bottoms = np.zeros(len(datasets_for_cath))

        for cls_idx, cls in enumerate([3, 1, 2, 4]):  # Order: Alpha-Beta first (largest)
            fracs = []
            for ds_name, ds_accs in datasets_for_cath:
                sub = cath[cath[cath_acc_col].isin(ds_accs)]
                total = len(sub)
                count = (sub["cath_class"] == cls).sum()
                fracs.append(count / total if total > 0 else 0)
            ax.bar(x_pos, fracs, bottom=bottoms, width=0.5,
                   label=CATH_CLASS_NAMES.get(cls, str(cls)),
                   color=CATH_CLASS_COLORS[cls_idx], edgecolor="white", linewidth=0.5)
            bottoms += fracs

        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"{name}\n(n={len(cath[cath[cath_acc_col].isin(accs)])})"
                            for name, accs in datasets_for_cath])
        ax.set_ylabel("Fraction")
        ax.set_title("CATH Class Distribution")
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
        ax.set_ylim(0, 1.05)

    # Panel B: Top superfamilies (horizontal bars)
    ax = fig.add_subplot(gs[1])
    add_panel_label(ax, "B")
    if cath is not None:
        cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
        sf_col = "cath_superfamily" if "cath_superfamily" in cath.columns else None
        if sf_col:
            groel_sf = cath[cath[cath_acc_col].isin(groel_acc)][sf_col].value_counts().head(8)
            hsp60_sf = cath[cath[cath_acc_col].isin(hsp60_acc)][sf_col].value_counts().head(8)
            all_sfs = list(dict.fromkeys(list(groel_sf.index) + list(hsp60_sf.index)))[:10]
            all_sfs.reverse()  # Largest at top

            y = np.arange(len(all_sfs))
            w = 0.35
            g_counts = [groel_sf.get(sf, 0) for sf in all_sfs]
            h_counts = [hsp60_sf.get(sf, 0) for sf in all_sfs]

            ax.barh(y + w/2, g_counts, w, label="GroEL",
                    color=DATASET_COLORS["groel"], edgecolor="white", linewidth=0.5)
            ax.barh(y - w/2, h_counts, w, label="HSP60",
                    color=DATASET_COLORS["hsp60"], edgecolor="white", linewidth=0.5)
            ax.set_yticks(y)
            ax.set_yticklabels(all_sfs, fontsize=8, fontfamily="monospace")
            ax.set_xlabel("Domain Count")
            ax.set_title("Top CATH Superfamilies")
            ax.legend(fontsize=8)

    # Panel C: Domain count distribution
    ax = fig.add_subplot(gs[2])
    add_panel_label(ax, "C")
    if dist is not None and len(dist) > 0:
        bar_width = 0.35
        for i, ds in enumerate(["groel", "hsp60"]):
            ds_data = dist[dist["dataset"] == ds]
            ds_data = ds_data[ds_data["n_domains"] <= 6]  # Truncate at 6 for readability
            if len(ds_data) > 0:
                offset = -bar_width/2 if i == 0 else bar_width/2
                n_total = ds_data["count"].sum()
                ax.bar(ds_data["n_domains"] + offset, ds_data["percent"],
                       width=bar_width, label=f"{DATASET_LABELS[ds]} (n={n_total})",
                       color=DATASET_COLORS[ds], edgecolor="white", linewidth=0.5,
                       alpha=0.85)
        ax.set_xlabel("Number of Domains")
        ax.set_ylabel("Percentage (%)")
        ax.set_title("Domain Count")
        ax.legend(fontsize=8)
        ax.set_xlim(-0.5, 6.5)
        ax.set_xticks(range(7))

    save_figure(fig, "fig1_domain_architecture")
    figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Figure 2: N-vs-C Structural Comparison
# ===========================================================================
print("\n--- Figure 2: N-vs-C Stability ---")
try:
    if paired is not None and len(paired) > 0:
        fig = plt.figure(figsize=(17, 5.5))
        gs = GridSpec(1, 3, width_ratios=[1, 1, 0.9], wspace=0.35)

        datasets_to_plot = ["groel", "hsp60", "matrix_bg", "mito_bg"]
        metrics_for_violin = [
            ("relative_contact_order", "Relative Contact Order"),
            ("mean_plddt", "Mean pLDDT"),
        ]

        for panel_idx, (metric, ylabel) in enumerate(metrics_for_violin):
            ax = fig.add_subplot(gs[panel_idx])
            add_panel_label(ax, chr(65 + panel_idx))
            plot_rows = []
            sample_sizes = {}
            for ds in datasets_to_plot:
                sub = paired[paired["datasets"].str.contains(ds)]
                n = 0
                for _, row in sub.iterrows():
                    n_val = row.get(f"{metric}_n_domain")
                    c_val = row.get(f"{metric}_c_region")
                    if pd.notna(n_val):
                        plot_rows.append({"Dataset": DATASET_LABELS[ds],
                                          "Region": "N-domain", "Value": n_val})
                        n += 1
                    if pd.notna(c_val):
                        plot_rows.append({"Dataset": DATASET_LABELS[ds],
                                          "Region": "C-region", "Value": c_val})
                sample_sizes[DATASET_LABELS[ds]] = n

            if plot_rows:
                pdf = pd.DataFrame(plot_rows)
                ds_order = [DATASET_LABELS[d] for d in datasets_to_plot]
                sns.violinplot(data=pdf, x="Dataset", y="Value", hue="Region",
                               split=True, ax=ax, palette=NC_COLORS, cut=0,
                               inner="quartile", order=ds_order, linewidth=0.8)
                ax.set_ylabel(ylabel)
                ax.set_title(ylabel)
                ax.set_xlabel("")

                # Add sample sizes below x-axis
                for i, ds_label in enumerate(ds_order):
                    ax.text(i, ax.get_ylim()[0], f"n={sample_sizes.get(ds_label, '?')}",
                            ha="center", va="top", fontsize=8, color="gray")

                if panel_idx == 0:
                    ax.legend(fontsize=8, loc="upper right")
                else:
                    leg = ax.get_legend()
                    if leg:
                        leg.remove()
                ax.tick_params(axis="x", rotation=20)

        # Panel C: Heatmap of N-C differences (median)
        ax = fig.add_subplot(gs[2])
        add_panel_label(ax, "C")
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
            hm_df = pd.DataFrame(heatmap_data,
                                 index=[DATASET_LABELS[ds] for ds in datasets_to_plot])
            short_names = {"relative_contact_order": "RCO", "mean_plddt": "pLDDT",
                           "frac_helix": "Helix", "frac_strand": "Strand",
                           "frac_hydrophobic": "Hydrophobic", "frac_plddt_gt70": "pLDDT>70"}
            hm_df.columns = [short_names.get(c, c) for c in hm_df.columns]
            sns.heatmap(hm_df.T, ax=ax, cmap="RdBu_r", center=0, annot=True,
                        fmt=".3f", linewidths=0.8, linecolor="white",
                        cbar_kws={"label": "N - C (median)", "shrink": 0.8},
                        annot_kws={"fontsize": 9})
            ax.set_title("N-C Difference")
            ax.tick_params(axis="y", rotation=0)

        save_figure(fig, "fig2_n_vs_c_stability")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Figure 3: GroEL Class-Specific Effects
# ===========================================================================
print("\n--- Figure 3: GroEL Class Effects ---")
try:
    if paired is not None and "groel_class" in paired.columns:
        groel_sub = paired[paired["datasets"].str.contains("groel")]
        groel_sub = groel_sub[groel_sub["groel_class"].isin(["I", "II", "III"])]

        if len(groel_sub) > 10:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

            class_colors = {"I": "#009E73", "II": "#D55E00", "III": "#56B4E9"}

            for idx, (metric, ylabel) in enumerate([
                ("relative_contact_order_diff", "RCO Difference (N - C)"),
                ("mean_plddt_diff", "pLDDT Difference (N - C)"),
            ]):
                ax = axes[idx]
                add_panel_label(ax, chr(65 + idx))
                if metric in groel_sub.columns:
                    order = ["I", "II", "III"]
                    bp = sns.boxplot(data=groel_sub, x="groel_class", y=metric,
                                order=order, ax=ax, palette=class_colors,
                                fliersize=2, linewidth=1.2)
                    sns.stripplot(data=groel_sub, x="groel_class", y=metric,
                                  order=order, ax=ax, color="black", alpha=0.3,
                                  size=3, jitter=True)
                    ax.axhline(0, color="red", linestyle="--", alpha=0.5, linewidth=1)
                    ax.set_xlabel("GroEL Substrate Class")
                    ax.set_ylabel(ylabel)
                    ax.set_title(ylabel)

                    # Add sample sizes from real data
                    for i, cls in enumerate(order):
                        n = len(groel_sub[groel_sub["groel_class"] == cls])
                        ax.text(i, ax.get_ylim()[1] * 0.95, f"n={n}", ha="center",
                                fontsize=9, style="italic", color="gray")

            # Add Kruskal-Wallis p-value from real statistics
            kw_key = "H2.3_GroEL_class_relative_contact_order"
            if kw_key in pval_lookup:
                p = pval_lookup[kw_key]["p"]
                axes[0].text(0.98, 0.02, f"Kruskal-Wallis p = {p:.2f}",
                            transform=axes[0].transAxes, ha="right", va="bottom",
                            fontsize=9, style="italic",
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))
            kw_key2 = "H2.3_GroEL_class_mean_plddt"
            if kw_key2 in pval_lookup:
                p = pval_lookup[kw_key2]["p"]
                axes[1].text(0.98, 0.02, f"Kruskal-Wallis p = {p:.2f}",
                            transform=axes[1].transAxes, ha="right", va="bottom",
                            fontsize=9, style="italic",
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

            plt.tight_layout()
            save_figure(fig, "fig3_groel_class_comparison")
            figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Figure 4: MTS Targeting
# ===========================================================================
print("\n--- Figure 4: MTS Targeting ---")
try:
    has_mts_data = (targeting is not None or mts_domain is not None)
    if has_mts_data:
        fig = plt.figure(figsize=(17, 5.5))
        gs = GridSpec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

        # Panel A: Targeting classification
        ax = fig.add_subplot(gs[0])
        add_panel_label(ax, "A")
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
                    # Use a sequential colormap based on count
                    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(counts)))
                    bars = counts.plot(kind="barh", ax=ax, color=colors,
                                       edgecolor="white", linewidth=0.5)
                    for i, v in enumerate(counts.values):
                        ax.text(v + 0.5, i, str(int(v)), va="center", fontsize=9)
                    ax.set_xlabel("Count")
                    ax.set_title("HSP60 Substrate Targeting")

                    # Add enrichment OR from real statistics
                    h3_key = "H3.1_HSP60_matrix_enrichment"
                    if h3_key in pval_lookup:
                        or_val = pval_lookup[h3_key]["effect"]
                        p_val = pval_lookup[h3_key]["p"]
                        ax.text(0.98, 0.98,
                                f"Matrix enrichment\nOR = {or_val:.2f}, p = {p_val:.1e}",
                                transform=ax.transAxes, ha="right", va="top",
                                fontsize=8,
                                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

        # Panel B: MTS-to-domain gap distribution
        ax = fig.add_subplot(gs[1])
        add_panel_label(ax, "B")
        if mts_domain is not None and "gap_length" in mts_domain.columns:
            gaps = mts_domain["gap_length"].dropna()
            if len(gaps) > 0:
                ax.hist(gaps, bins=30, color=DATASET_COLORS["mito_bg"],
                        edgecolor="white", alpha=0.85, linewidth=0.5)
                ax.axvline(0, color="red", linestyle="--", linewidth=1.5,
                           label="Gap = 0")
                median_gap = gaps.median()
                ax.axvline(median_gap, color="orange", linestyle="--", linewidth=1.5,
                           label=f"Median = {median_gap:.0f}")
                ax.set_xlabel("Gap Length (residues)")
                ax.set_ylabel("Count")
                ax.set_title("MTS-to-Domain Gap")
                ax.legend(fontsize=8)

                pre_col = "mts_is_pre_domain" if "mts_is_pre_domain" in mts_domain.columns else None
                if pre_col:
                    n_pre = int(mts_domain[pre_col].sum())
                    n_total = len(mts_domain[mts_domain[pre_col].notna()])
                    pct = 100 * n_pre / n_total
                    ax.text(0.97, 0.97,
                            f"Pre-domain: {n_pre}/{n_total}\n({pct:.1f}%)",
                            transform=ax.transAxes, ha="right", va="top", fontsize=9,
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

        # Panel C: MTS cleavage site vs first domain start
        ax = fig.add_subplot(gs[2])
        add_panel_label(ax, "C")
        if mts_domain is not None:
            tp_col = "transit_peptide_end" if "transit_peptide_end" in mts_domain.columns else None
            fd_col = "first_domain_start" if "first_domain_start" in mts_domain.columns else None
            if tp_col and fd_col:
                valid = mts_domain[[tp_col, fd_col]].dropna()
                if len(valid) > 0:
                    ax.scatter(valid[tp_col], valid[fd_col], alpha=0.5,
                               s=20, color=DATASET_COLORS["groel"],
                               edgecolors="gray", linewidth=0.3)
                    lim = max(valid[tp_col].max(), valid[fd_col].max()) * 1.1
                    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, linewidth=1,
                            label="y = x")
                    ax.set_xlabel("MTS Cleavage Site (residue)")
                    ax.set_ylabel("First Domain Start (residue)")
                    ax.set_title("MTS vs Domain Boundary")
                    ax.legend(fontsize=8, loc="lower right")
                    ax.text(0.05, 0.95, f"n = {len(valid)}",
                            transform=ax.transAxes, va="top", fontsize=9, color="gray")

        save_figure(fig, "fig4_mts_targeting")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Figure 5: Cross-Species Orthology
# ===========================================================================
print("\n--- Figure 5: Orthology ---")
try:
    if homologs is not None and paired is not None:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

        # Panel A: Orthogroup categories
        ax = axes[0]
        add_panel_label(ax, "A")
        if ortho is not None:
            cat_col = "category" if "category" in ortho.columns else None
            if cat_col:
                counts = ortho[cat_col].value_counts()
                # Clean category names
                clean_names = {
                    "hsp60_only": "HSP60 only",
                    "groel_only": "GroEL only",
                    "shared": "Shared",
                }
                clean_counts = counts.rename(index=clean_names)
                # Order: shared, groel_only, hsp60_only
                order = ["GroEL only", "Shared", "HSP60 only"]
                ordered = clean_counts.reindex([o for o in order if o in clean_counts.index])

                bar_colors = [DATASET_COLORS["groel"], DATASET_COLORS["matrix_bg"],
                              DATASET_COLORS["hsp60"]]
                bars = ordered.plot(kind="bar", ax=ax, color=bar_colors[:len(ordered)],
                                    edgecolor="white", linewidth=0.5)
                for i, v in enumerate(ordered.values):
                    ax.text(i, v + 1, str(int(v)), ha="center", fontsize=10, fontweight="bold")
                ax.set_ylabel("Number of Orthogroups")
                ax.set_title("Substrate Orthogroup Categories")
                ax.tick_params(axis="x", rotation=0)
                ax.set_xlabel("")

        # Panel B: N-domain RCO conservation across homologs
        ax = axes[1]
        add_panel_label(ax, "B")
        groel_col = "groel_accession" if "groel_accession" in homologs.columns else None
        hsp60_col = "hsp60_accession" if "hsp60_accession" in homologs.columns else None

        if groel_col and hsp60_col:
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
                ax.scatter(scatter_x, scatter_y, alpha=0.6, s=35,
                           color=DATASET_COLORS["groel"], edgecolors="gray", linewidth=0.5)

                # Identity line
                all_vals = scatter_x + scatter_y
                lim_lo = min(all_vals) * 0.9
                lim_hi = max(all_vals) * 1.1
                ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", alpha=0.3, linewidth=1)

                # Real correlation from actual data
                r, p = stats.pearsonr(scatter_x, scatter_y)
                ax.text(0.05, 0.95,
                        f"Pearson r = {r:.2f}\np = {p:.1e}\nn = {len(scatter_x)} pairs",
                        transform=ax.transAxes, va="top", fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8))
                ax.set_xlabel("GroEL Homolog N-domain RCO")
                ax.set_ylabel("HSP60 Homolog N-domain RCO")
                ax.set_title("N-Domain Contact Order Conservation")
            else:
                ax.text(0.5, 0.5,
                        f"Insufficient homolog pairs with RCO data (n={len(scatter_x)})",
                        transform=ax.transAxes, ha="center", va="center")
                ax.set_title("N-Domain RCO Conservation")

        plt.tight_layout()
        save_figure(fig, "fig5_orthology")
        figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Figure 6: Key Findings Summary
# ===========================================================================
print("\n--- Figure 6: Summary ---")
try:
    fig = plt.figure(figsize=(17, 5.5))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.8], wspace=0.35)

    # Panel A: Mean N-C RCO difference by dataset with significance
    ax = fig.add_subplot(gs[0])
    add_panel_label(ax, "A")
    if paired is not None and "relative_contact_order_diff" in paired.columns:
        ds_labels = []
        ds_means = []
        ds_sems = []
        ds_ns = []
        ds_pvals = []
        for ds in ["groel", "hsp60", "matrix_bg", "mito_bg"]:
            sub = paired[paired["datasets"].str.contains(ds)]
            diffs = sub["relative_contact_order_diff"].dropna()
            if len(diffs) >= 5:
                ds_labels.append(DATASET_LABELS[ds])
                ds_means.append(diffs.mean())
                ds_sems.append(diffs.sem())
                ds_ns.append(len(diffs))
                # Get real p-value
                key = f"H2.1_{ds}_relative_contact_order"
                if key in pval_lookup:
                    ds_pvals.append(pval_lookup[key]["p"])
                else:
                    ds_pvals.append(np.nan)

        if ds_labels:
            colors = [DATASET_COLORS.get(dl.lower(), DATASET_COLORS.get(
                {v: k for k, v in DATASET_LABELS.items()}.get(dl, ""), "#999"))
                for dl in ds_labels]
            bars = ax.bar(range(len(ds_labels)), ds_means, yerr=ds_sems,
                   color=colors, capsize=5, edgecolor="gray", linewidth=0.5)
            ax.axhline(0, color="black", linestyle="-", linewidth=0.5)
            ax.set_xticks(range(len(ds_labels)))
            ax.set_xticklabels(ds_labels, rotation=20)
            ax.set_ylabel("Mean RCO Diff (N - C)")
            ax.set_title("N-C Contact Order Asymmetry")

            # Add significance stars from real p-values
            for i, (p, n) in enumerate(zip(ds_pvals, ds_ns)):
                if pd.notna(p):
                    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    y_pos = ds_means[i] + ds_sems[i] + 0.003
                    ax.text(i, y_pos, stars, ha="center", fontsize=10, fontweight="bold")
                    ax.text(i, -0.005, f"n={n}", ha="center", fontsize=7, color="gray")

    # Panel B: Statistical significance overview
    ax = fig.add_subplot(gs[1])
    add_panel_label(ax, "B")
    if pvalues is not None and len(pvalues) > 0:
        sig_col = "significant_overall"
        if sig_col in pvalues.columns:
            family_names = {
                "domain_architecture": "Domain\nArchitecture",
                "stability_asymmetry": "Stability\nAsymmetry",
                "mts_targeting": "MTS\nTargeting"
            }
            sig_counts = pvalues.groupby("family")[sig_col].agg(["sum", "count"])
            sig_counts.columns = ["Significant", "Total"]
            sig_counts["Not Significant"] = sig_counts["Total"] - sig_counts["Significant"]
            sig_counts.index = [family_names.get(f, f) for f in sig_counts.index]

            y_pos = range(len(sig_counts))
            sig_vals = sig_counts["Significant"].values
            ns_vals = sig_counts["Not Significant"].values

            ax.barh(y_pos, sig_vals, color="#D55E00", edgecolor="white",
                    linewidth=0.5, label="Significant")
            ax.barh(y_pos, ns_vals, left=sig_vals, color="#CCCCCC",
                    edgecolor="white", linewidth=0.5, label="Not Significant")

            # Add count labels
            for i in range(len(sig_counts)):
                if sig_vals[i] > 0:
                    ax.text(sig_vals[i]/2, i, str(int(sig_vals[i])),
                            ha="center", va="center", fontsize=10, fontweight="bold", color="white")
                total = sig_vals[i] + ns_vals[i]
                ax.text(total + 0.5, i, f"/{int(total)}", va="center", fontsize=9)

            ax.set_yticks(y_pos)
            ax.set_yticklabels(sig_counts.index)
            ax.set_xlabel("Number of Tests")
            ax.set_title("Statistical Significance\n(hierarchical BH, FDR < 0.05)")
            ax.legend(fontsize=8, loc="lower right")

    # Panel C: MTS-domain relationship pie
    ax = fig.add_subplot(gs[2])
    add_panel_label(ax, "C")
    if mts_domain is not None:
        pre_col = "mts_is_pre_domain" if "mts_is_pre_domain" in mts_domain.columns else None
        if pre_col:
            valid = mts_domain[mts_domain[pre_col].notna()]
            n_pre = int(valid[pre_col].sum())
            n_overlap = len(valid) - n_pre
            wedges, texts, autotexts = ax.pie(
                [n_pre, n_overlap],
                labels=[f"Pre-domain\n(n={n_pre})", f"Overlaps domain\n(n={n_overlap})"],
                colors=[DATASET_COLORS["matrix_bg"], DATASET_COLORS["hsp60"]],
                autopct="%1.1f%%", startangle=90,
                wedgeprops=dict(edgecolor="white", linewidth=1.5))
            for t in autotexts:
                t.set_fontsize(11)
                t.set_fontweight("bold")
            ax.set_title("MTS-Domain Relationship")

            # Add binomial p-value from real statistics
            h32_key = "H3.2_MTS_pre_domain"
            if h32_key in pval_lookup:
                p = pval_lookup[h32_key]["p"]
                ax.text(0.5, -0.1, f"Binomial p = {p:.1e}",
                        transform=ax.transAxes, ha="center", fontsize=9, style="italic")

    save_figure(fig, "fig6_summary")
    figures_generated += 1
except Exception as e:
    print(f"  Error: {e}")
    import traceback; traceback.print_exc()


# ===========================================================================
# Final summary
# ===========================================================================
print(f"\n{'=' * 70}")
print(f"Module I (Polished) complete. Generated {figures_generated}/6 figures.")
print(f"Output directory: {FIG_DIR}")
print(f"All data sourced from real Phase 2 result files — no fabricated data.")
print(f"{'=' * 70}")
