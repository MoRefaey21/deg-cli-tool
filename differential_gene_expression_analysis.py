"""
Differential Gene Expression Analysis Script
============================================

Description
-----------
Identifies differentially expressed genes (DEGs) between a control and a
treatment group from a gene expression matrix. Calculates fold change and
statistical significance and produces a volcano plot plus a structured report.

The appropriate statistical test is selected automatically based on whether
the samples are paired and whether the data meets the normality assumption.

Features
--------
- Interactive prompts for file loading, group definition, and thresholds
- Automatic statistical test selection:
    - Paired + normal     -> paired t-test
    - Paired + non-normal -> Wilcoxon signed-rank
    - Unpaired + normal   -> Welch's t-test (unequal variances)
    - Unpaired + non-normal -> Mann-Whitney U
- Shapiro-Wilk normality screening on a randomly sampled subset of genes
- Benjamini-Hochberg FDR correction across all tested genes; significance is
  called on the adjusted p-value (q-value), not the raw p-value
- Fold change and log2 fold change. If the input is already log2-scaled
  (e.g. RMA-normalized microarray data) log2 FC is the difference of group
  means; otherwise it is log2 of the ratio of group means
- Volcano plot (genes coloured only when they meet BOTH FC and FDR thresholds)
- Tab-separated DEG report documenting the test selection rationale

Requirements
------------
NumPy, SciPy, Matplotlib, tqdm

Usage notes
-----------
The script will prompt for:
- Data file path and format parameters
- Indices of control and treatment columns (0-based)
- Whether the samples are paired
- Whether the expression values are already log2-transformed
- Significance threshold (FDR) and fold change threshold
- Output report file name
"""

import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import tqdm

from log2_fc import log2_fc
from log10_p import neg_log10_pvalue
from volcano_plot import volcano_plot


print("=== Differential Gene Expression Analysis ===\n")

print("Running notes:")
print("1. Ensure your data file has gene IDs in the first column and sample")
print("   expression values in the remaining columns.")
print("2. The script will ask whether your samples are paired (e.g. the same")
print("   subjects measured before and after treatment).")
print("3. Based on the pairing structure and a normality check on a subset of")
print("   genes, the appropriate statistical test is selected automatically:")
print("   - Paired   + normal     -> paired t-test")
print("   - Paired   + non-normal -> Wilcoxon signed-rank")
print("   - Unpaired + normal     -> Welch's t-test")
print("   - Unpaired + non-normal -> Mann-Whitney U\n")


# =====================================================================
# Data loading
# =====================================================================
while True:
    print("=== Data Loading ===")
    try:
        data_file = input(
            "Enter the file name (with extension, same directory) or full path: "
        ).strip()
        delimiter = input(
            "Enter the file delimiter (e.g. ',' for CSV, '\\t' for tab-separated files): "
        ).replace("\\t", "\t").strip()
        skip_header = int(input(
            "Enter the number of lines to skip at the top of the file "
            "(1 for a single header row, 0 for none): "
        ).strip())
        comments_input = input(
            "Enter the comment symbol to skip (press Enter if none): "
        ).strip()
        comments = comments_input if comments_input else None

        GE_data = np.genfromtxt(
            data_file,
            delimiter=delimiter,
            skip_header=skip_header,
            dtype=str,
            comments=comments,
        )

        if GE_data.ndim != 2:
            raise ValueError(
                "File loaded as a 1D array. Check the delimiter and file format."
            )

        print(f"\nData loaded successfully: {GE_data.shape[0]} rows, "
              f"{GE_data.shape[1]} columns.\n")
        break

    except Exception as e:
        print(f"\nError loading data: {e}")
        print("Please re-enter the parameters correctly.\n")


# =====================================================================
# Gene IDs and expression values
# =====================================================================
gene_ids = np.char.strip(GE_data[:, 0], '"')
print(f"Extracted {len(gene_ids)} gene IDs.\n")

try:
    expression_values = GE_data[:, 1:].astype(float)
except ValueError as e:
    print(f"Error converting expression values to float: {e}")
    print("Check that all columns after the gene ID column contain numeric data only.")
    raise SystemExit(1)


# =====================================================================
# Group definition
# =====================================================================
print("=== Group Definition ===")
n_samples = expression_values.shape[1]
print(f"Detected {n_samples} sample columns (after the gene ID column).\n")

while True:
    try:
        control_group = list(map(
            int,
            input("Enter control group indices (comma-separated, 0 = first sample): ")
            .strip().split(",")
        ))
        treatment_group = list(map(
            int,
            input("Enter treatment group indices (comma-separated, 0 = first sample): ")
            .strip().split(",")
        ))

        all_idx = control_group + treatment_group
        if any(i < 0 or i >= n_samples for i in all_idx):
            raise ValueError(f"Indices must be between 0 and {n_samples - 1}.")
        if set(control_group) & set(treatment_group):
            raise ValueError("Control and treatment indices must not overlap.")

        print(f"\nControl group: {len(control_group)} samples. "
              f"Treatment group: {len(treatment_group)} samples.\n")
        break

    except ValueError as e:
        print(f"\nInvalid input: {e}\nPlease re-enter valid comma-separated integers.\n")

control_values  = expression_values[:, control_group]
treatment_values = expression_values[:, treatment_group]

# Fold change is only meaningful on the linear scale. Microarray series
# matrices (RMA output) are normally log2-scaled, so ask rather than assume.
while True:
    log_input = input(
        "Are the expression values already log2-transformed? (y/n): "
    ).strip().lower()
    if log_input in ("y", "n"):
        is_log2 = log_input == "y"
        break
    print("Please enter 'y' or 'n'.")


# =====================================================================
# Automatic statistical test selection
# =====================================================================
print("=== Statistical Test Selection ===")

# Step 1 - paired or unpaired?
while True:
    paired_input = input(
        "Are your samples paired (e.g. same subjects before and after treatment)? (y/n): "
    ).strip().lower()
    if paired_input in ("y", "n"):
        is_paired = paired_input == "y"
        break
    print("Please enter 'y' or 'n'.")

if is_paired and len(control_group) != len(treatment_group):
    print(
        f"\nYou indicated paired samples, but the group sizes differ "
        f"({len(control_group)} vs {len(treatment_group)}). "
        f"Switching to unpaired analysis.\n"
    )
    is_paired = False

# Step 2 - normality check
n_genes = expression_values.shape[0]
min_group_size = min(len(control_group), len(treatment_group))
normality_note = ""

if min_group_size < 5:
    # Shapiro-Wilk is unreliable at very small n; apply parametric by convention
    is_normal = True
    normality_note = (
        f"Sample size per group is {min_group_size}. Normality testing is unreliable "
        f"at this sample size. Parametric test applied by convention."
    )
    print(f"\n  Note: only {min_group_size} samples per group. Normality testing "
          f"requires at least 5 observations per group to be meaningful.")
    print(f"  Applying parametric test by convention.\n")
else:
    sample_size = min(200, n_genes)
    rng = np.random.default_rng(seed=42)
    sampled = rng.choice(n_genes, size=sample_size, replace=False)

    print(f"\nChecking normality on a random sample of {sample_size} genes "
          f"(Shapiro-Wilk, alpha = 0.05)...")
    non_normal, checked = 0, 0

    for i in sampled:
        if is_paired:
            diff = treatment_values[i] - control_values[i]
            if np.std(diff) > 0 and len(diff) >= 3:
                checked += 1
                _, sw_p = stat.shapiro(diff)
                if sw_p < 0.05:
                    non_normal += 1
        else:
            for group_vals in (control_values[i], treatment_values[i]):
                if np.std(group_vals) > 0 and len(group_vals) >= 3:
                    checked += 1
                    _, sw_p = stat.shapiro(group_vals)
                    if sw_p < 0.05:
                        non_normal += 1

    if checked > 0:
        pct = 100.0 * non_normal / checked
        print(f"  {pct:.1f}% of checked samples reject normality.")
        is_normal = pct <= 30
        normality_note = (
            f"{pct:.1f}% of {sample_size} sampled genes reject normality "
            f"(Shapiro-Wilk, alpha=0.05). "
            f"{'Parametric' if is_normal else 'Non-parametric'} test selected."
        )
        if is_normal:
            print("  Normality assumption acceptable. Selecting parametric test.")
        else:
            print("  Normality assumption not met. Selecting non-parametric test.")
    else:
        is_normal = True
        normality_note = (
            "Normality could not be assessed (insufficient variance in sampled genes). "
            "Parametric test applied by default."
        )
        print("  Could not assess normality (no variance in sampled genes). "
              "Using parametric test by default.")

# Step 3 - select test
if is_paired and is_normal:
    test_type = "paired_t"
elif is_paired and not is_normal:
    test_type = "wilcoxon"
elif not is_paired and is_normal:
    test_type = "welch_t"
else:
    test_type = "mannwhitney"

# Optional override, so that a pre-specified analysis plan can fix the test
# instead of letting the normality screen choose it.
override = input(
    "Test override (Enter for automatic, or paired_t / wilcoxon / welch_t / "
    "mannwhitney): "
).strip().lower()
if override:
    if override not in ("paired_t", "wilcoxon", "welch_t", "mannwhitney"):
        raise SystemExit(f"Unknown test '{override}'.")
    if override in ("paired_t", "wilcoxon") and len(control_group) != len(treatment_group):
        raise SystemExit("Paired tests need equally sized groups.")
    test_type = override
    normality_note += f" Test fixed by the user: {test_type}."

print(f"\nSelected test: {test_type}\n")


# =====================================================================
# Run the test
# =====================================================================
test_funcs = {
    "paired_t":    stat.ttest_rel,
    "welch_t":     lambda a, b: stat.ttest_ind(a, b, equal_var=False),
    "wilcoxon":    stat.wilcoxon,
    "mannwhitney": stat.mannwhitneyu,
}
test_func = test_funcs[test_type]

test_statistic = np.full(n_genes, np.nan)
p_values       = np.full(n_genes, np.nan)

for i in tqdm.tqdm(range(n_genes), desc=f"Running {test_type}"):
    try:
        result = test_func(control_values[i], treatment_values[i])
        test_statistic[i] = result[0]
        p_values[i]       = result[1]
    except Exception:
        pass

def benjamini_hochberg(p):
    """Return BH-adjusted p-values (q-values); NaN inputs stay NaN."""
    q = np.full(p.shape, np.nan)
    valid = np.where(~np.isnan(p))[0]
    order = valid[np.argsort(p[valid])]
    m = len(order)
    ranked = p[order] * m / np.arange(1, m + 1)
    # enforce monotonicity from the largest p-value downwards
    q[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return q


q_values = benjamini_hochberg(p_values)

n_significant = int(np.sum(p_values < 0.05))
n_significant_fdr = int(np.sum(q_values < 0.05))
print(f"\n{test_type} completed.")
print(f"Number of genes with raw p < 0.05: {n_significant}")
print(f"Number of genes with FDR (BH) q < 0.05: {n_significant_fdr}\n")


# =====================================================================
# Fold change calculation
# =====================================================================
print("=== Fold Change Calculation ===")

mean_control   = np.mean(control_values,  axis=1)
mean_treatment = np.mean(treatment_values, axis=1)

with np.errstate(divide="ignore", invalid="ignore"):
    if is_log2:
        # difference of log2 means = log2 of the ratio of geometric means
        log2_fold_change_array = mean_treatment - mean_control
        fold_change_array      = 2.0 ** log2_fold_change_array
    else:
        fold_change_array      = np.where(mean_control > 0,
                                          mean_treatment / mean_control,
                                          np.nan)
        log2_fold_change_array = log2_fc(fold_change_array)

n_with_fc    = int(np.sum(~np.isnan(fold_change_array)))
n_without_fc = int(np.sum( np.isnan(fold_change_array)))
print(f"  Genes with computable fold change: {n_with_fc}")
print(f"  Genes skipped (zero or negative control mean): {n_without_fc}\n")


# =====================================================================
# DEG extraction
# =====================================================================
print("=== DEGs Extraction ===")
try:
    alpha = float(input("Enter FDR threshold for q-values (standard: 0.05): "))
    fc_threshold = float(input("Enter fold change threshold (standard: 1.5): "))
    if not (0 < alpha < 1):
        raise ValueError("Alpha must be between 0 and 1.")
    if fc_threshold <= 1:
        raise ValueError("Fold change threshold must be greater than 1.")
except ValueError as e:
    print(f"Invalid input ({e}). Using defaults: alpha = 0.05, fc_threshold = 1.5.")
    alpha = 0.05
    fc_threshold = 1.5

down_threshold = 1.0 / fc_threshold

DEGs_pvalue_indices = np.where(q_values < alpha)[0]
n_raw_p = int(np.sum(p_values < alpha))

with np.errstate(invalid="ignore"):
    fc_up   = fold_change_array >= fc_threshold
    fc_down = fold_change_array <= down_threshold

DEGs_fc_indices      = np.where((fc_up | fc_down) & ~np.isnan(fold_change_array))[0]
DEGs_overlap_indices = np.intersect1d(DEGs_pvalue_indices, DEGs_fc_indices)

over_expressed  = int(np.sum(fc_up))
under_expressed = int(np.sum(fc_down))

print(f"  Genes with raw p < {alpha} (uncorrected, for reference): {n_raw_p}")
print(f"  Total DEGs by FDR q-value: {len(DEGs_pvalue_indices)}")
print(f"  Total DEGs by fold change: {len(DEGs_fc_indices)}")
print(f"  Over-expressed (FC threshold): {over_expressed}")
print(f"  Under-expressed (FC threshold): {under_expressed}")
print(f"  Overlap (both criteria): {len(DEGs_overlap_indices)}\n")


# =====================================================================
# Volcano plot
# =====================================================================
print("=== Volcano Plot Generation ===")
print("Calculating -log10 p-values...")
log10_pvalue = neg_log10_pvalue(p_values)
print("Done.\n")

finite_log2fc = log2_fold_change_array[np.isfinite(log2_fold_change_array)]
finite_log10p = log10_pvalue[np.isfinite(log10_pvalue)]

x_bound = max(np.abs(finite_log2fc).max() if finite_log2fc.size else 1.0, 1.0)
y_bound = max(finite_log10p.max()          if finite_log10p.size else 5.0, 5.0)

print("Generating the volcano plot...")
volcano_plot(
    log2_fold_change_array,
    log10_pvalue,
    log2fc_threshold=np.log2(fc_threshold),
    xlim=(-x_bound * 1.1, x_bound * 1.1),
    ylim=(0, y_bound * 1.1),
    significant=(q_values < alpha),
    alpha=alpha,
)


# =====================================================================
# Report generation
# =====================================================================
print("\n=== Report File Generation ===")
output_report = input("Enter the output report file name (without extension): ").strip() + ".txt"


def _classify(i):
    fc = fold_change_array[i]
    if np.isnan(fc):
        return "FC not computable"
    if fc >= fc_threshold:
        return "Over-expressed"
    if fc <= down_threshold:
        return "Under-expressed"
    return "Within FC thresholds"


def _row(i):
    return (
        f"{gene_ids[i]}\t"
        f"{fold_change_array[i]:.6g}\t"
        f"{log2_fold_change_array[i]:.6g}\t"
        f"{p_values[i]:.6g}\t"
        f"{q_values[i]:.6g}\t"
        f"{test_statistic[i]:.6g}\t"
        f"{_classify(i)}\n"
    )


print(f"Saving results to {output_report}...")
try:
    with open(output_report, "w") as report:
        report.write("# Differential Gene Expression Analysis Report\n\n")
        report.write(f"# Data file: {data_file}\n")
        report.write(f"# Paired analysis: {is_paired}\n")
        report.write(f"# Normality assessment: {normality_note}\n")
        report.write(f"# Statistical test selected: {test_type}\n")
        report.write(f"# Input log2-transformed: {is_log2}\n")
        report.write(f"# FDR (Benjamini-Hochberg) threshold: {alpha}\n")
        report.write(f"# Fold change threshold: {fc_threshold}\n")
        report.write(f"# Total genes analyzed: {len(gene_ids)}\n")
        report.write(f"# Genes with raw p < {alpha} (uncorrected): {n_raw_p}\n")
        report.write(f"# DEGs by FDR q-value: {len(DEGs_pvalue_indices)}\n")
        report.write(f"# DEGs by fold change: {len(DEGs_fc_indices)}\n")
        report.write(f"# Over-expressed (FC): {over_expressed}\n")
        report.write(f"# Under-expressed (FC): {under_expressed}\n")
        report.write(f"# Overlapping DEGs: {len(DEGs_overlap_indices)}\n\n")

        header = "Gene_ID\tFC\tlog2FC\tP_value\tFDR_q_value\tTest_statistic\tChange_Type\n"

        report.write("## DEGs by FDR q-value\n")
        report.write(header)
        for i in DEGs_pvalue_indices:
            report.write(_row(i))
        report.write("\n")

        report.write("## DEGs by fold change\n")
        report.write(header)
        for i in DEGs_fc_indices:
            report.write(_row(i))
        report.write("\n")

        report.write("## Overlapping DEGs\n")
        report.write(header)
        for i in DEGs_overlap_indices:
            report.write(_row(i))

    print(f"Analysis complete. Report saved to {output_report}.")

except Exception as e:
    print(f"\nError saving report: {e}")
    print("Try a different output filename or check your write permissions.")
