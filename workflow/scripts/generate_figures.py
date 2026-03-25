#!/usr/bin/env python3
"""
Generate all publication-quality figures for the Antah Asti Prarambh project.
Comparative structural proteomics of chaperonin substrates (GroEL/HSP60).

Three goals:
  1. Domain architecture comparison
  2. N-domain vs C-region stability asymmetry
  3. Mitochondrial targeting and MTS-domain relationship

Output: PDF + PNG figures at 300 DPI in results/figures/
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"
OUT  = os.path.join(BASE, "results", "figures")
os.makedirs(OUT, exist_ok=True)

# Input files
CATH_PROTEIN  = os.path.join(BASE, "results/domains/cath_protein_summary.tsv")
CATH_DOMAINS  = os.path.join(BASE, "results/domains/cath_domain_assignments.tsv")
DOMAIN_METRICS = os.path.join(BASE, "results/domains/domain_structural_metrics.tsv")
FOLDSEEK      = os.path.join(BASE, "results/domains/foldseek_clusters.tsv")
NVC_PAIRED    = os.path.join(BASE, "results/termini/n_vs_c_paired.tsv")
CONTACT_ORDER = os.path.join(BASE, "results/termini/contact_order.tsv")
REGION_BOUNDS = os.path.join(BASE, "results/termini/region_boundaries.tsv")
SEQ_METRICS   = os.path.join(BASE, "results/termini/sequence_metrics.tsv")
STRUCT_METRICS = os.path.join(BASE, "results/termini/structure_metrics.tsv")
COMBINED_TARG = os.path.join(BASE, "results/mts/combined_targeting.tsv")
MTS_DOMAIN    = os.path.join(BASE, "results/mts/mts_domain_relationship.tsv")
RBH_ANNOT     = os.path.join(BASE, "results/homology/rbh_groel_hsp60_annotated.tsv")
ORTHOGROUPS   = os.path.join(BASE, "results/homology/substrate_orthogroups.tsv")
GROEL_SUBS    = os.path.join(BASE, "data/processed/groel_substrates_standardized.tsv")
HSP60_SUBS    = os.path.join(BASE, "data/processed/hsp60_tier1_substrates.tsv")
STRUCT_INDEX  = os.path.join(BASE, "results/structures/structure_index.tsv")

# ---------------------------------------------------------------------------
# Style setup
# ---------------------------------------------------------------------------
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Colorblind-friendly palette
CB_PALETTE = sns.color_palette("Set2", 8)
DATASET_COLORS = {
    "GroEL": CB_PALETTE[0],
    "HSP60": CB_PALETTE[1],
    "Mito": CB_PALETTE[2],
    "Matrix": CB_PALETTE[3],
}
CATH_CLASS_COLORS = {
    "Mainly Alpha": CB_PALETTE[4],
    "Mainly Beta": CB_PALETTE[5],
    "Alpha Beta": CB_PALETTE[6],
    "Few SS": CB_PALETTE[7],
}
NC_COLORS = {"N-domain": CB_PALETTE[0], "C-region": CB_PALETTE[1]}

# ---------------------------------------------------------------------------
# Helper: assign dataset label from source_dataset string
# ---------------------------------------------------------------------------
def assign_dataset(source_str):
    """Return a list of dataset labels for a source_dataset string."""
    s = str(source_str).lower()
    labels = []
    if "groel" in s:
        labels.append("GroEL")
    if "hsp60" in s:
        labels.append("HSP60")
    if "matrix" in s:
        labels.append("Matrix")
    if "mito" in s:
        labels.append("Mito")
    return labels if labels else ["Other"]


def assign_primary_dataset(source_str):
    """Assign single primary dataset for bar charts where we need unique assignment."""
    s = str(source_str).lower()
    if "groel" in s:
        return "GroEL"
    if "hsp60" in s:
        return "HSP60"
    if "matrix" in s:
        return "Matrix"
    if "mito" in s:
        return "Mito"
    return "Other"


def cath_class_name(code):
    """Map CATH class number to readable name."""
    mapping = {
        "1": "Mainly Alpha",
        "2": "Mainly Beta",
        "3": "Alpha Beta",
        "4": "Few SS",
        "6": "Few SS",
    }
    return mapping.get(str(code), "Other")


SUPERFAMILY_NAMES = {
    "3.20.20.70": "TIM Barrel",
    "3.40.50.720": "Rossmann Fold",
    "3.40.50.300": "P-loop NTPase",
    "3.40.50.150": "Thioredoxin-like",
    "1.10.510.10": "Transferase",
    "3.50.50.60": "FAD/NAD(P)-binding",
    "1.50.40.10": "Beta-Propeller",
    "3.40.30.10": "Aldolase Class I",
    "1.25.40.10": "Leucine-rich Repeat",
    "3.40.640.10": "Flavodoxin-like",
    "3.90.1150.10": "PLP-dependent Enz.",
    "1.10.238.10": "Helix-Hairpin-Helix",
    "3.90.226.10": "Alkaline Phosphatase",
    "3.40.50.620": "HUP-like",
    "3.30.830.10": "Aldolase Class II",
    "2.40.50.100": "OB Fold",
    "3.40.50.1820": "CoA-dependent Acyltransf.",
    "2.40.30.10": "Lipocalin",
    "1.10.10.10": "Arc Repressor",
    "3.40.47.10": "Periplasmic Binding Prot.",
    "3.30.930.10": "BirA Bifunctional",
}


def get_sf_name(code):
    return SUPERFAMILY_NAMES.get(code, code)


def save_figure(fig, name):
    """Save figure as PDF and PNG, print file sizes."""
    for ext in (".pdf", ".png"):
        path = os.path.join(OUT, name + ext)
        fig.savefig(path, dpi=300, bbox_inches="tight")
        size_kb = os.path.getsize(path) / 1024
        print(f"  Saved: {path}  ({size_kb:.1f} KB)")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data...")

# Domain data
cath_protein = pd.read_csv(CATH_PROTEIN, sep="\t")
cath_domains = pd.read_csv(CATH_DOMAINS, sep="\t")
domain_metrics = pd.read_csv(DOMAIN_METRICS, sep="\t")
foldseek = pd.read_csv(FOLDSEEK, sep="\t")

# Termini data
nvc_paired = pd.read_csv(NVC_PAIRED, sep="\t")
contact_order = pd.read_csv(CONTACT_ORDER, sep="\t")
region_bounds = pd.read_csv(REGION_BOUNDS, sep="\t")
seq_metrics = pd.read_csv(SEQ_METRICS, sep="\t")
struct_metrics = pd.read_csv(STRUCT_METRICS, sep="\t")

# MTS data
combined_targ = pd.read_csv(COMBINED_TARG, sep="\t")
mts_domain = pd.read_csv(MTS_DOMAIN, sep="\t")

# Homology data
rbh_annot = pd.read_csv(RBH_ANNOT, sep="\t")
orthogroups = pd.read_csv(ORTHOGROUPS, sep="\t")

# Substrate lists
groel_subs = pd.read_csv(GROEL_SUBS, sep="\t")
hsp60_subs = pd.read_csv(HSP60_SUBS, sep="\t")

# Structure index - for dataset mapping
struct_index = pd.read_csv(STRUCT_INDEX, sep="\t")

# Build a unified accession-to-datasets mapping from structure_index
acc2ds = {}
for _, row in struct_index.iterrows():
    acc = row["uniprot_accession"]
    ds_list = assign_dataset(row["source_dataset"])
    acc2ds[acc] = ds_list

print(f"  Proteins in structure index: {len(struct_index)}")
print(f"  Domain assignments: {len(cath_domains)}")
print(f"  N-vs-C paired entries: {len(nvc_paired)}")
print(f"  MTS targeting entries: {len(combined_targ)}")
print(f"  MTS-domain relationship entries: {len(mts_domain)}")
print(f"  RBH annotated pairs: {len(rbh_annot)}")
print(f"  Orthogroups: {len(orthogroups)}")
print()

# ===================================================================
# FIGURE 1: Domain Architecture Overview
# ===================================================================
print("=" * 60)
print("FIGURE 1: Domain Architecture Overview")
print("=" * 60)

fig1 = plt.figure(figsize=(16, 5.5))
gs1 = GridSpec(1, 3, figure=fig1, wspace=0.35)

# --- Fig 1A: Stacked bar of CATH class distribution ---
ax1a = fig1.add_subplot(gs1[0, 0])

# Explode domain assignments by dataset
rows_1a = []
for _, row in cath_domains.iterrows():
    acc = row["uniprot_accession"]
    cls_name = cath_class_name(row["cath_class"])
    datasets = acc2ds.get(acc, [assign_primary_dataset(row.get("source", ""))])
    for ds in datasets:
        if ds in ("GroEL", "HSP60", "Mito", "Matrix"):
            rows_1a.append({"dataset": ds, "cath_class": cls_name})

df_1a = pd.DataFrame(rows_1a)
ct_1a = pd.crosstab(df_1a["dataset"], df_1a["cath_class"], normalize="index")
# Reorder columns
class_order = ["Mainly Alpha", "Mainly Beta", "Alpha Beta", "Few SS"]
ct_1a = ct_1a.reindex(columns=[c for c in class_order if c in ct_1a.columns], fill_value=0)
# Reorder rows
ds_order = ["GroEL", "HSP60", "Mito", "Matrix"]
ct_1a = ct_1a.reindex([d for d in ds_order if d in ct_1a.index])

colors_1a = [CATH_CLASS_COLORS[c] for c in ct_1a.columns]
ct_1a.plot(kind="bar", stacked=True, ax=ax1a, color=colors_1a, edgecolor="white", linewidth=0.5)
ax1a.set_ylabel("Proportion")
ax1a.set_xlabel("")
ax1a.set_title("A. CATH Class Distribution", fontweight="bold", loc="left")
ax1a.set_ylim(0, 1.05)
ax1a.legend(title="CATH Class", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title_fontsize=9)
ax1a.set_xticklabels(ax1a.get_xticklabels(), rotation=0, ha="center")

# --- Fig 1B: Horizontal bar of top 10 superfamilies GroEL vs HSP60 ---
ax1b = fig1.add_subplot(gs1[0, 1])

# Get per-dataset superfamily counts
rows_1b = []
for _, row in cath_domains.iterrows():
    acc = row["uniprot_accession"]
    sf = str(row["cath_superfamily"])
    datasets = acc2ds.get(acc, [])
    for ds in datasets:
        if ds in ("GroEL", "HSP60"):
            rows_1b.append({"dataset": ds, "superfamily": sf})

df_1b = pd.DataFrame(rows_1b)
# Top superfamilies overall
top_sf = df_1b["superfamily"].value_counts().head(10).index.tolist()
df_1b_top = df_1b[df_1b["superfamily"].isin(top_sf)]

sf_counts = df_1b_top.groupby(["dataset", "superfamily"]).size().reset_index(name="count")
sf_counts["sf_name"] = sf_counts["superfamily"].apply(get_sf_name)
sf_pivot = sf_counts.pivot(index="sf_name", columns="dataset", values="count").fillna(0)
# Sort by total
sf_pivot["total"] = sf_pivot.sum(axis=1)
sf_pivot = sf_pivot.sort_values("total", ascending=True).drop(columns="total")
# Ensure GroEL and HSP60 columns
for c in ["GroEL", "HSP60"]:
    if c not in sf_pivot.columns:
        sf_pivot[c] = 0

y_pos = np.arange(len(sf_pivot))
bar_h = 0.35
ax1b.barh(y_pos - bar_h / 2, sf_pivot["GroEL"], bar_h, label="GroEL",
           color=DATASET_COLORS["GroEL"], edgecolor="white", linewidth=0.5)
ax1b.barh(y_pos + bar_h / 2, sf_pivot["HSP60"], bar_h, label="HSP60",
           color=DATASET_COLORS["HSP60"], edgecolor="white", linewidth=0.5)
ax1b.set_yticks(y_pos)
ax1b.set_yticklabels(sf_pivot.index, fontsize=8)
ax1b.set_xlabel("Domain Count")
ax1b.set_title("B. Top 10 CATH Superfamilies", fontweight="bold", loc="left")
ax1b.legend(fontsize=9)

# --- Fig 1C: Domain count distribution ---
ax1c = fig1.add_subplot(gs1[0, 2])

# Map each protein in cath_protein to datasets, get n_domains
rows_1c = []
for _, row in cath_protein.iterrows():
    acc = row["uniprot_accession"]
    n_dom = row["n_domains"]
    datasets = acc2ds.get(acc, [])
    for ds in datasets:
        if ds in ("GroEL", "HSP60", "Mito", "Matrix"):
            label = "Single-domain" if n_dom == 1 else ("Multi-domain" if n_dom >= 2 else "No domain")
            rows_1c.append({"dataset": ds, "domain_class": label, "n_domains": n_dom})

df_1c = pd.DataFrame(rows_1c)
# Only keep proteins with at least 1 domain
df_1c_filt = df_1c[df_1c["n_domains"] > 0]
ct_1c = pd.crosstab(df_1c_filt["dataset"], df_1c_filt["domain_class"], normalize="index")
cat_order = ["Single-domain", "Multi-domain"]
ct_1c = ct_1c.reindex(columns=[c for c in cat_order if c in ct_1c.columns], fill_value=0)
ct_1c = ct_1c.reindex([d for d in ds_order if d in ct_1c.index])

colors_1c = [CB_PALETTE[0], CB_PALETTE[1]]
ct_1c.plot(kind="bar", stacked=True, ax=ax1c, color=colors_1c, edgecolor="white", linewidth=0.5)
ax1c.set_ylabel("Proportion")
ax1c.set_xlabel("")
ax1c.set_title("C. Domain Count Distribution", fontweight="bold", loc="left")
ax1c.set_ylim(0, 1.05)
ax1c.legend(title="", fontsize=9)
ax1c.set_xticklabels(ax1c.get_xticklabels(), rotation=0, ha="center")

fig1.suptitle("Figure 1: Domain Architecture of Chaperonin Substrates", fontsize=14, fontweight="bold", y=1.02)
save_figure(fig1, "fig1_domain_architecture")

# ===================================================================
# FIGURE 2: N-domain vs C-region Comparison
# ===================================================================
print("\n" + "=" * 60)
print("FIGURE 2: N-domain vs C-region Stability Comparison")
print("=" * 60)

# Assign primary dataset to nvc_paired
nvc = nvc_paired.copy()
nvc["dataset"] = nvc["source_dataset"].apply(assign_primary_dataset)
nvc = nvc[nvc["dataset"].isin(["GroEL", "HSP60", "Mito", "Matrix"])]

fig2 = plt.figure(figsize=(18, 6))
gs2 = GridSpec(1, 3, figure=fig2, wspace=0.35)

# --- Fig 2A: Paired violin of relative contact order ---
ax2a = fig2.add_subplot(gs2[0, 0])

# Prepare long-form data for contact order
rco_long = []
for _, row in nvc.iterrows():
    ds = row["dataset"]
    n_val = row.get("relative_contact_order_n_domain", np.nan)
    c_val = row.get("relative_contact_order_c_region", np.nan)
    if pd.notna(n_val):
        rco_long.append({"Dataset": ds, "Region": "N-domain", "Relative Contact Order": n_val})
    if pd.notna(c_val):
        rco_long.append({"Dataset": ds, "Region": "C-region", "Relative Contact Order": c_val})
rco_df = pd.DataFrame(rco_long)

if len(rco_df) > 0:
    sns.violinplot(data=rco_df, x="Dataset", y="Relative Contact Order", hue="Region",
                   split=True, inner="quart", palette=NC_COLORS, ax=ax2a, order=ds_order,
                   cut=0, linewidth=1)
    # Add jittered dots (subsample if too many)
    for i, ds in enumerate(ds_order):
        for j, reg in enumerate(["N-domain", "C-region"]):
            subset = rco_df[(rco_df["Dataset"] == ds) & (rco_df["Region"] == reg)]
            if len(subset) > 0:
                jitter = np.random.normal(0, 0.03, size=len(subset))
                x_pos = i + (j - 0.5) * 0.15
                ax2a.scatter(x_pos + jitter, subset["Relative Contact Order"],
                            alpha=0.15, s=8, color=NC_COLORS[reg], zorder=3, edgecolors="none")

ax2a.set_title("A. Relative Contact Order", fontweight="bold", loc="left")
ax2a.set_xlabel("")
ax2a.legend(title="Region", fontsize=9)

# --- Fig 2B: Paired violin of mean pLDDT ---
ax2b = fig2.add_subplot(gs2[0, 1])

plddt_long = []
for _, row in nvc.iterrows():
    ds = row["dataset"]
    n_val = row.get("mean_plddt_n_domain", np.nan)
    c_val = row.get("mean_plddt_c_region", np.nan)
    if pd.notna(n_val):
        plddt_long.append({"Dataset": ds, "Region": "N-domain", "Mean pLDDT": n_val})
    if pd.notna(c_val):
        plddt_long.append({"Dataset": ds, "Region": "C-region", "Mean pLDDT": c_val})
plddt_df = pd.DataFrame(plddt_long)

if len(plddt_df) > 0:
    sns.violinplot(data=plddt_df, x="Dataset", y="Mean pLDDT", hue="Region",
                   split=True, inner="quart", palette=NC_COLORS, ax=ax2b, order=ds_order,
                   cut=0, linewidth=1)
    for i, ds in enumerate(ds_order):
        for j, reg in enumerate(["N-domain", "C-region"]):
            subset = plddt_df[(plddt_df["Dataset"] == ds) & (plddt_df["Region"] == reg)]
            if len(subset) > 0:
                jitter = np.random.normal(0, 0.03, size=len(subset))
                x_pos = i + (j - 0.5) * 0.15
                ax2b.scatter(x_pos + jitter, subset["Mean pLDDT"],
                            alpha=0.15, s=8, color=NC_COLORS[reg], zorder=3, edgecolors="none")

ax2b.set_title("B. Mean pLDDT (Confidence)", fontweight="bold", loc="left")
ax2b.set_xlabel("")
ax2b.legend(title="Region", fontsize=9)

# --- Fig 2C: Heatmap of N-C differences across metrics and datasets ---
ax2c = fig2.add_subplot(gs2[0, 2])

metrics_for_heatmap = [
    ("relative_contact_order", "Rel. Contact Order"),
    ("mean_plddt", "Mean pLDDT"),
    ("frac_helix", "Frac. Helix"),
    ("frac_strand", "Frac. Strand"),
    ("frac_hydrophobic", "Frac. Hydrophobic"),
    ("frac_plddt_gt70", "Frac. pLDDT>70"),
]

heatmap_data = []
pval_data = []
for metric_base, metric_label in metrics_for_heatmap:
    n_col = f"{metric_base}_n_domain"
    c_col = f"{metric_base}_c_region"
    row_vals = []
    row_pvals = []
    for ds in ds_order:
        sub = nvc[nvc["dataset"] == ds].dropna(subset=[n_col, c_col])
        if len(sub) > 5:
            diff = sub[n_col] - sub[c_col]
            mean_diff = diff.mean()
            # Paired t-test (or Wilcoxon)
            try:
                _, pval = stats.wilcoxon(sub[n_col], sub[c_col])
            except Exception:
                pval = 1.0
        else:
            mean_diff = np.nan
            pval = 1.0
        row_vals.append(mean_diff)
        row_pvals.append(pval)
    heatmap_data.append(row_vals)
    pval_data.append(row_pvals)

hm_df = pd.DataFrame(heatmap_data, index=[m[1] for m in metrics_for_heatmap], columns=ds_order)
pv_df = pd.DataFrame(pval_data, index=[m[1] for m in metrics_for_heatmap], columns=ds_order)

# Create annotation strings with significance stars
annot = hm_df.copy().astype(str)
for i in range(hm_df.shape[0]):
    for j in range(hm_df.shape[1]):
        val = hm_df.iloc[i, j]
        pv = pv_df.iloc[i, j]
        stars = ""
        if pv < 0.001:
            stars = "***"
        elif pv < 0.01:
            stars = "**"
        elif pv < 0.05:
            stars = "*"
        if pd.notna(val):
            annot.iloc[i, j] = f"{val:.3f}{stars}"
        else:
            annot.iloc[i, j] = "n/a"

# Diverging colormap: red = N higher, blue = C higher
vmax = max(abs(hm_df.min().min()), abs(hm_df.max().max()))
if pd.isna(vmax) or vmax == 0:
    vmax = 1
sns.heatmap(hm_df, annot=annot, fmt="", ax=ax2c, cmap="RdBu_r",
            center=0, vmin=-vmax, vmax=vmax, linewidths=0.5,
            cbar_kws={"label": "Mean Difference (N - C)", "shrink": 0.8},
            annot_kws={"fontsize": 8})
ax2c.set_title("C. N-C Difference Heatmap", fontweight="bold", loc="left")
ax2c.set_ylabel("")
ax2c.set_xlabel("")

fig2.suptitle("Figure 2: N-domain vs C-region Structural Comparison", fontsize=14, fontweight="bold", y=1.02)
save_figure(fig2, "fig2_n_vs_c_stability")

# ===================================================================
# FIGURE 3: GroEL Class Comparison
# ===================================================================
print("\n" + "=" * 60)
print("FIGURE 3: GroEL Class Comparison")
print("=" * 60)

# Filter for GroEL substrates with class info
groel_nvc = nvc_paired.copy()
groel_nvc = groel_nvc[groel_nvc["source_dataset"].str.contains("groel", case=False, na=False)]
groel_nvc = groel_nvc[groel_nvc["groel_class"].isin(["I", "II", "III"])]

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5.5))

class_order = ["I", "II", "III"]
class_colors = [CB_PALETTE[0], CB_PALETTE[1], CB_PALETTE[2]]

# --- Fig 3A: RCO N-C difference by class ---
if "relative_contact_order_diff" in groel_nvc.columns and len(groel_nvc) > 0:
    sns.boxplot(data=groel_nvc, x="groel_class", y="relative_contact_order_diff",
                order=class_order, palette=class_colors, ax=ax3a,
                fliersize=3, linewidth=1.2)
    sns.stripplot(data=groel_nvc, x="groel_class", y="relative_contact_order_diff",
                  order=class_order, color="black", alpha=0.3, size=4, jitter=True, ax=ax3a)
    ax3a.axhline(0, color="gray", linestyle="--", linewidth=0.8, zorder=0)
    ax3a.set_xlabel("GroEL Class")
    ax3a.set_ylabel("RCO Difference (N - C)")
    ax3a.set_title("A. Contact Order Asymmetry by GroEL Class", fontweight="bold", loc="left")

    # Add n counts
    for i, cls in enumerate(class_order):
        n = len(groel_nvc[groel_nvc["groel_class"] == cls])
        ax3a.text(i, ax3a.get_ylim()[0] + 0.01 * (ax3a.get_ylim()[1] - ax3a.get_ylim()[0]),
                  f"n={n}", ha="center", fontsize=9, color="gray")

# --- Fig 3B: pLDDT N-C difference by class ---
if "mean_plddt_diff" in groel_nvc.columns and len(groel_nvc) > 0:
    sns.boxplot(data=groel_nvc, x="groel_class", y="mean_plddt_diff",
                order=class_order, palette=class_colors, ax=ax3b,
                fliersize=3, linewidth=1.2)
    sns.stripplot(data=groel_nvc, x="groel_class", y="mean_plddt_diff",
                  order=class_order, color="black", alpha=0.3, size=4, jitter=True, ax=ax3b)
    ax3b.axhline(0, color="gray", linestyle="--", linewidth=0.8, zorder=0)
    ax3b.set_xlabel("GroEL Class")
    ax3b.set_ylabel("pLDDT Difference (N - C)")
    ax3b.set_title("B. pLDDT Asymmetry by GroEL Class", fontweight="bold", loc="left")

    for i, cls in enumerate(class_order):
        n = len(groel_nvc[groel_nvc["groel_class"] == cls])
        ax3b.text(i, ax3b.get_ylim()[0] + 0.01 * (ax3b.get_ylim()[1] - ax3b.get_ylim()[0]),
                  f"n={n}", ha="center", fontsize=9, color="gray")

fig3.suptitle("Figure 3: GroEL Class-Specific N/C Asymmetry", fontsize=14, fontweight="bold", y=1.02)
fig3.tight_layout()
save_figure(fig3, "fig3_groel_class_comparison")

# ===================================================================
# FIGURE 4: Mitochondrial Targeting
# ===================================================================
print("\n" + "=" * 60)
print("FIGURE 4: Mitochondrial Targeting")
print("=" * 60)

fig4 = plt.figure(figsize=(16, 5.5))
gs4 = GridSpec(1, 3, figure=fig4, wspace=0.35)

# --- Fig 4A: Stacked bar of targeting classifications ---
ax4a = fig4.add_subplot(gs4[0, 0])

# Filter for HSP60 substrates (proteins that are hsp60 substrates)
targ = combined_targ.copy()
# The combined_targeting.tsv has all mitochondrial proteins, filter for HSP60 substrates
hsp60_mask = targ["is_hsp60_substrate"] == True
if hsp60_mask.sum() == 0:
    # Try string comparison
    hsp60_mask = targ["is_hsp60_substrate"].astype(str).str.lower() == "true"

targ_hsp60 = targ[hsp60_mask] if hsp60_mask.sum() > 0 else targ

# Simplify classification labels
class_map = {
    "High-confidence matrix": "High-conf. Matrix",
    "Non-canonical matrix import": "Non-canonical Matrix",
    "Inner membrane (MIM)": "Inner Membrane",
    "Intermembrane space (IMS)": "IMS",
    "Outer membrane (MOM)": "Outer Membrane",
    "Non-mitochondrial / no targeting signal": "Non-mito.",
    "Mitochondrial (other/unspecified)": "Mito. (other)",
    "Probable matrix": "Probable Matrix",
    "Mitochondrial (with MTS)": "Mito. (with MTS)",
    "Secretory pathway (signal peptide)": "Secretory",
    "Predicted mitochondrial (transit peptide)": "Predicted Mito.",
}

if len(targ_hsp60) > 0:
    targ_counts = targ_hsp60["targeting_classification"].value_counts()
else:
    # Use all data instead
    targ_counts = targ["targeting_classification"].value_counts()

targ_counts_mapped = targ_counts.rename(index=class_map)
# Sort by count
targ_counts_mapped = targ_counts_mapped.sort_values(ascending=True)

# Color by category
targ_color_map = {
    "High-conf. Matrix": "#2ca02c",
    "Non-canonical Matrix": "#98df8a",
    "Probable Matrix": "#7fbf7f",
    "Inner Membrane": "#ff7f0e",
    "IMS": "#ffbb78",
    "Outer Membrane": "#d62728",
    "Non-mito.": "#9467bd",
    "Mito. (other)": "#c5b0d5",
    "Mito. (with MTS)": "#8c564b",
    "Secretory": "#e377c2",
    "Predicted Mito.": "#bcbd22",
}
colors_4a = [targ_color_map.get(c, "#7f7f7f") for c in targ_counts_mapped.index]

targ_counts_mapped.plot(kind="barh", ax=ax4a, color=colors_4a, edgecolor="white", linewidth=0.5)
ax4a.set_xlabel("Number of Proteins")
ax4a.set_title("A. Targeting Classification", fontweight="bold", loc="left")
# Annotate bar ends with counts
for i, (val, name) in enumerate(zip(targ_counts_mapped.values, targ_counts_mapped.index)):
    ax4a.text(val + max(targ_counts_mapped) * 0.02, i, str(val), va="center", fontsize=9)

subtitle = "HSP60 substrates" if hsp60_mask.sum() > 0 else "All mitochondrial proteins"
ax4a.set_ylabel("")

# --- Fig 4B: Histogram of gap distance ---
ax4b = fig4.add_subplot(gs4[0, 1])

mts_dom = mts_domain.copy()
gap = mts_dom["gap_length"].dropna()

if len(gap) > 0:
    ax4b.hist(gap, bins=30, color=CB_PALETTE[2], edgecolor="white", linewidth=0.5, alpha=0.85)
    ax4b.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Zero (MTS = Domain start)")
    ax4b.set_xlabel("Gap Distance (residues)")
    ax4b.set_ylabel("Count")
    ax4b.set_title("B. MTS-to-Domain Gap Distance", fontweight="bold", loc="left")
    ax4b.legend(fontsize=9)
    # Add annotations
    n_pre = (gap > 0).sum()
    n_overlap = (gap < 0).sum()
    n_zero = (gap == 0).sum()
    ax4b.text(0.95, 0.95, f"Pre-domain: {n_pre}\nOverlap: {n_overlap}\nExact: {n_zero}",
              transform=ax4b.transAxes, ha="right", va="top", fontsize=9,
              bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.7))

# --- Fig 4C: Scatter of MTS cleavage vs domain start ---
ax4c = fig4.add_subplot(gs4[0, 2])

mts_x = mts_dom["transit_peptide_end"].dropna()
mts_y = mts_dom["first_domain_start"].dropna()
# Use intersection of valid entries
valid = mts_dom.dropna(subset=["transit_peptide_end", "first_domain_start"])

if len(valid) > 0:
    ax4c.scatter(valid["transit_peptide_end"], valid["first_domain_start"],
                 alpha=0.5, s=25, color=CB_PALETTE[1], edgecolors="gray", linewidth=0.3)
    # y=x line
    lims = [min(valid["transit_peptide_end"].min(), valid["first_domain_start"].min()) - 5,
            max(valid["transit_peptide_end"].max(), valid["first_domain_start"].max()) + 5]
    ax4c.plot(lims, lims, "k--", linewidth=1, alpha=0.6, label="y = x")
    ax4c.set_xlim(lims)
    ax4c.set_ylim(lims)
    ax4c.set_xlabel("MTS Cleavage Site (residue)")
    ax4c.set_ylabel("First Domain Start (residue)")
    ax4c.set_title("C. MTS vs First Domain", fontweight="bold", loc="left")
    ax4c.legend(fontsize=9)

    # Count above/below
    n_above = (valid["first_domain_start"] > valid["transit_peptide_end"]).sum()
    n_below = (valid["first_domain_start"] < valid["transit_peptide_end"]).sum()
    ax4c.text(0.05, 0.95,
              f"Above y=x (pre-domain): {n_above}\nBelow y=x (overlap): {n_below}",
              transform=ax4c.transAxes, ha="left", va="top", fontsize=9,
              bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.7))

fig4.suptitle("Figure 4: Mitochondrial Targeting and MTS-Domain Relationship", fontsize=14, fontweight="bold", y=1.02)
save_figure(fig4, "fig4_mts_targeting")

# ===================================================================
# FIGURE 5: Orthology Overview
# ===================================================================
print("\n" + "=" * 60)
print("FIGURE 5: Orthology Overview")
print("=" * 60)

fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(14, 6))

# --- Fig 5A: Venn diagram or UpSet ---
try:
    from matplotlib_venn import venn3

    og = orthogroups.copy()
    groel_only = og[og["category"] == "groel_only"]
    hsp60_only = og[og["category"] == "hsp60_only"]
    shared = og[og["category"] == "shared"]

    # Collect all accessions per set
    groel_accs = set()
    hsp60_accs = set()
    for _, row in og.iterrows():
        g_acc = str(row.get("groel_accessions", ""))
        h_acc = str(row.get("hsp60_accessions", ""))
        if g_acc and g_acc != "nan":
            groel_accs.update(g_acc.split(","))
        if h_acc and h_acc != "nan":
            hsp60_accs.update(h_acc.split(","))

    # Build sets based on orthogroup membership
    shared_groel = set()
    shared_hsp60 = set()
    for _, row in shared.iterrows():
        g_acc = str(row.get("groel_accessions", ""))
        h_acc = str(row.get("hsp60_accessions", ""))
        if g_acc and g_acc != "nan":
            shared_groel.update(g_acc.split(","))
        if h_acc and h_acc != "nan":
            shared_hsp60.update(h_acc.split(","))

    # For Venn: use orthogroup counts
    n_groel_only = len(groel_only)
    n_hsp60_only = len(hsp60_only)
    n_shared = len(shared)

    v = venn3(subsets=(n_groel_only, n_hsp60_only, n_shared, 0, 0, 0, 0),
              set_labels=("GroEL\nOrthogroups", "HSP60\nOrthogroups", ""),
              ax=ax5a)
    # Customize: use simpler 2-set representation
    ax5a.clear()
    from matplotlib_venn import venn2
    v2 = venn2(subsets=(n_groel_only, n_hsp60_only, n_shared),
               set_labels=("GroEL\nOrthogroups", "HSP60\nOrthogroups"),
               ax=ax5a,
               set_colors=(DATASET_COLORS["GroEL"], DATASET_COLORS["HSP60"]),
               alpha=0.6)
    ax5a.set_title("A. Orthogroup Overlap", fontweight="bold", loc="left")
    ax5a.text(0.5, -0.1, f"Shared: {n_shared} orthogroups ({len(shared_groel)} GroEL, {len(shared_hsp60)} HSP60 proteins)",
              transform=ax5a.transAxes, ha="center", fontsize=9, style="italic")

except ImportError:
    # Fallback: bar chart
    og = orthogroups.copy()
    cat_counts = og["category"].value_counts()
    cat_counts.plot(kind="bar", ax=ax5a, color=[DATASET_COLORS["GroEL"], DATASET_COLORS["HSP60"], CB_PALETTE[3]],
                    edgecolor="white", linewidth=0.5)
    ax5a.set_title("A. Orthogroup Categories", fontweight="bold", loc="left")
    ax5a.set_ylabel("Number of Orthogroups")
    ax5a.set_xlabel("")
    ax5a.set_xticklabels(ax5a.get_xticklabels(), rotation=0)

# --- Fig 5B: Scatter of GroEL vs HSP60 homolog properties ---
# Use RBH pairs and merge with n_vs_c data for contact order
rbh = rbh_annot.copy()

# Merge GroEL-side metrics
groel_co = contact_order[contact_order["region"] == "n_domain"][["uniprot_accession", "relative_contact_order"]].rename(
    columns={"uniprot_accession": "groel_accession", "relative_contact_order": "groel_rco_n"})
hsp60_co = contact_order[contact_order["region"] == "n_domain"][["uniprot_accession", "relative_contact_order"]].rename(
    columns={"uniprot_accession": "hsp60_accession", "relative_contact_order": "hsp60_rco_n"})

rbh_merged = rbh.merge(groel_co, on="groel_accession", how="left").merge(hsp60_co, on="hsp60_accession", how="left")
rbh_valid = rbh_merged.dropna(subset=["groel_rco_n", "hsp60_rco_n"])

if len(rbh_valid) > 0:
    class_color_map = {"I": CB_PALETTE[0], "II": CB_PALETTE[1], "III": CB_PALETTE[2],
                       "I or II": CB_PALETTE[3]}
    colors = [class_color_map.get(str(c), "gray") for c in rbh_valid["groel_class"]]
    ax5b.scatter(rbh_valid["groel_rco_n"], rbh_valid["hsp60_rco_n"],
                 c=colors, alpha=0.7, s=50, edgecolors="gray", linewidth=0.4, zorder=3)

    # Diagonal
    lims = [0, max(rbh_valid["groel_rco_n"].max(), rbh_valid["hsp60_rco_n"].max()) * 1.1]
    ax5b.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5)
    ax5b.set_xlim(lims)
    ax5b.set_ylim(lims)

    # Correlation
    r, p = stats.pearsonr(rbh_valid["groel_rco_n"], rbh_valid["hsp60_rco_n"])
    ax5b.text(0.05, 0.95, f"r = {r:.2f}, p = {p:.2e}\nn = {len(rbh_valid)} pairs",
              transform=ax5b.transAxes, ha="left", va="top", fontsize=9,
              bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))

    # Legend
    for cls, col in class_color_map.items():
        n_cls = (rbh_valid["groel_class"] == cls).sum()
        if n_cls > 0:
            ax5b.scatter([], [], c=col, label=f"Class {cls} (n={n_cls})", s=40)
    ax5b.legend(title="GroEL Class", fontsize=8, title_fontsize=9, loc="lower right")
else:
    ax5b.text(0.5, 0.5, "Insufficient data\nfor RCO comparison", ha="center", va="center",
              transform=ax5b.transAxes, fontsize=12)

ax5b.set_xlabel("GroEL Substrate RCO (N-domain)")
ax5b.set_ylabel("HSP60 Homolog RCO (N-domain)")
ax5b.set_title("B. Conservation of N-domain Contact Order", fontweight="bold", loc="left")

fig5.suptitle("Figure 5: Cross-Species Orthology of Chaperonin Substrates", fontsize=14, fontweight="bold", y=1.02)
fig5.tight_layout()
save_figure(fig5, "fig5_orthology")

# ===================================================================
# FIGURE 6: Summary Figure
# ===================================================================
print("\n" + "=" * 60)
print("FIGURE 6: Summary Figure")
print("=" * 60)

fig6 = plt.figure(figsize=(15, 5))
gs6 = GridSpec(1, 3, figure=fig6, wspace=0.4)

# --- Fig 6A: Domain architecture key enrichment (simplified stacked bar) ---
ax6a = fig6.add_subplot(gs6[0, 0])

# Show top 5 enriched superfamilies in GroEL vs HSP60
if len(df_1b) > 0:
    groel_sf = df_1b[df_1b["dataset"] == "GroEL"]["superfamily"].value_counts()
    hsp60_sf = df_1b[df_1b["dataset"] == "HSP60"]["superfamily"].value_counts()
    # Compute enrichment ratio (GroEL / HSP60)
    all_sf = set(groel_sf.index) | set(hsp60_sf.index)
    enrichment = {}
    for sf in all_sf:
        g = groel_sf.get(sf, 0.5)
        h = hsp60_sf.get(sf, 0.5)
        enrichment[sf] = g / h
    enrich_df = pd.Series(enrichment).sort_values(ascending=False)
    # Top 5 GroEL-enriched
    top5 = enrich_df.head(5)
    top5_names = [get_sf_name(s) for s in top5.index]
    bars = ax6a.barh(range(len(top5)), top5.values, color=DATASET_COLORS["GroEL"],
                     edgecolor="white", linewidth=0.5, alpha=0.85)
    ax6a.set_yticks(range(len(top5)))
    ax6a.set_yticklabels(top5_names, fontsize=9)
    ax6a.axvline(1, color="gray", linestyle="--", linewidth=0.8)
    ax6a.set_xlabel("GroEL / HSP60 Ratio")
    ax6a.set_title("A. GroEL-Enriched Folds", fontweight="bold", loc="left")

# --- Fig 6B: N-vs-C contact order (strongest signal) ---
ax6b = fig6.add_subplot(gs6[0, 1])

# Compute mean N-C RCO difference per dataset
nvc_summary = []
for ds in ds_order:
    sub = nvc[nvc["dataset"] == ds]
    diff_col = "relative_contact_order_diff"
    if diff_col in sub.columns:
        vals = sub[diff_col].dropna()
        if len(vals) > 0:
            mean_d = vals.mean()
            sem_d = vals.sem()
            nvc_summary.append({"Dataset": ds, "Mean RCO Diff": mean_d, "SEM": sem_d})

nvc_sum_df = pd.DataFrame(nvc_summary)
if len(nvc_sum_df) > 0:
    colors_6b = [DATASET_COLORS.get(d, "gray") for d in nvc_sum_df["Dataset"]]
    ax6b.bar(nvc_sum_df["Dataset"], nvc_sum_df["Mean RCO Diff"],
             yerr=nvc_sum_df["SEM"], capsize=4, color=colors_6b,
             edgecolor="white", linewidth=0.5, alpha=0.85)
    ax6b.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax6b.set_ylabel("Mean RCO Diff (N - C)")
    ax6b.set_xlabel("")
    ax6b.set_title("B. N-C Contact Order Asymmetry", fontweight="bold", loc="left")

# --- Fig 6C: MTS is pre-domain (pie chart) ---
ax6c = fig6.add_subplot(gs6[0, 2])

if len(mts_domain) > 0:
    n_pre = mts_domain["mts_is_pre_domain"].sum()
    n_overlap = mts_domain["mts_overlaps_domain"].sum()
    n_other = len(mts_domain) - n_pre - n_overlap

    sizes = [n_pre, n_overlap]
    labels_pie = [f"Pre-domain\n(n={n_pre})", f"MTS overlaps\ndomain (n={n_overlap})"]
    colors_pie = [CB_PALETTE[2], CB_PALETTE[4]]
    if n_other > 0:
        sizes.append(n_other)
        labels_pie.append(f"Other\n(n={n_other})")
        colors_pie.append(CB_PALETTE[7])

    wedges, texts, autotexts = ax6c.pie(sizes, labels=labels_pie, colors=colors_pie,
                                         autopct="%1.1f%%", startangle=90,
                                         textprops={"fontsize": 9},
                                         pctdistance=0.6)
    for at in autotexts:
        at.set_fontsize(9)
        at.set_fontweight("bold")
    ax6c.set_title("C. MTS-Domain Relationship", fontweight="bold", loc="left")

fig6.suptitle("Figure 6: Key Findings Summary", fontsize=14, fontweight="bold", y=1.02)
save_figure(fig6, "fig6_summary")

# ===================================================================
# Done
# ===================================================================
print("\n" + "=" * 60)
print("ALL FIGURES GENERATED SUCCESSFULLY")
print("=" * 60)
print(f"\nOutput directory: {OUT}")
print(f"Total figures: 6 (12 files: 6 PDF + 6 PNG)")
