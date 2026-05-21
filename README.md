# Transcriptomic Analysis of B Cells in Multiple Sclerosis: Differential Gene Expression Before and After Apheresis Treatment

A Python command-line tool for differential gene expression (DEG) analysis, applied to CD19+ B cells from multiple sclerosis (MS) patients before and after apheresis treatment. Completed as the final project for CIT656: Programming to Bioinformatics, Postgraduate Bioinformatics Diploma at Nile University (2024-2025).

**Authors:** Mohammed Ashraf Ahmed, Ahmad Karam Ahmad, Alaa Hamed Mohamed, Omnia Esmail Ebrahim

## Question

Therapeutic apheresis (plasma exchange and immunoadsorption) is used to manage acute MS relapses by removing circulating autoantibodies and inflammatory mediators, but its molecular effects on immune cells remain poorly characterized. This project asks how the transcriptomic profile of CD19+ B cells changes between two timepoints in the same patients (immediately before the first apheresis session and immediately after the last), using a paired design to isolate treatment-induced changes from inter-individual variability.

## Pipeline

Input data are pre-processed gene-level expression profiles from GEO accession GSE272973 (Clariom D microarray), restricted to CD19+ B cell samples (3 paired before-and-after pairs). The user supplies the file path, delimiter, header and comment options, and the indices of control and treatment samples. The script then runs an interactive workflow consisting of: automated statistical test selection (the user answers one yes/no question about whether the samples are paired; the appropriate test is then chosen based on pairing structure and a Shapiro-Wilk normality check on a random subset of genes); gene-wise hypothesis testing with a progress bar; fold change and log2 fold change calculation with safe handling of zero-valued controls; DEG extraction by p-value, by fold change, and by their intersection; volcano plot generation (genes coloured only when they meet both the fold change and significance thresholds); and a structured tab-separated report documenting the test selection rationale.

## Results

From 135,750 genes analyzed with the paired t-test on the MS dataset, 7,193 genes met p < 0.05, 102 genes met |FC| > 1.5 (85 over-expressed and 17 under-expressed), and 16 genes met both criteria. The intersection list provides a high-confidence set of B-cell transcriptional changes following apheresis. The volcano plot visualizes the asymmetric pattern of activation over suppression in the post-treatment transcriptome, consistent with a rebound response after the removal of circulating inhibitory factors. The full annotated DEG list and figures are available in the project report PDF.

## Tools

Python 3.7+ with NumPy (data handling and vectorized math), SciPy (hypothesis testing and Shapiro-Wilk normality screening), Matplotlib (volcano plot), and tqdm (progress bars).

## Repository contents

- `differential_gene_expression_analysis.py`: main script
- `log2_fc.py`: helper for log2 fold change calculation
- `log10_p.py`: helper for -log10(p-value) calculation
- `volcano_plot.py`: helper for volcano plot generation
- `volcano_plot.png`: example output from the MS dataset
- `DEGs_report.txt`: example annotated DEG report from the MS dataset
- `Transcriptomic_Analysis_of_B_Cells_in_Multiple_Sclerosis.pdf`: full scientific report
- `Transcriptomic_Analysis_of_B_Cells_in_Multiple_Sclerosis.pptx`: project presentation

## Usage

```
python differential_gene_expression_analysis.py
```

The script will prompt for: data file path, delimiter, header and comment settings, control and treatment column indices, whether the samples are paired, alpha and fold change thresholds, and an output report file name. For the MS dataset used in the project: paired = yes, alpha = 0.05, fold change threshold = 1.5.

## Data

The dataset is publicly available from GEO under accession **GSE272973**. The B cell subset comprises 6 samples (3 paired pre-apheresis and 3 paired post-apheresis). The script accepts any tab- or comma-separated expression matrix with gene IDs in the first column and sample expression values in the remaining columns. The raw data file is not included in this repository.
