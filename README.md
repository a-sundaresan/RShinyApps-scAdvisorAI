# scAdvisorAI

![Language](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Framework](https://img.shields.io/badge/Framework-Shiny-blue?style=flat-square)
![AI](https://img.shields.io/badge/AI-Anthropic%20Claude-orange?style=flat-square)
![Input](https://img.shields.io/badge/Input-10X%20CellRanger-lightgrey?style=flat-square)
![Analysis](https://img.shields.io/badge/Analysis-scRNA--seq%20QC-green?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)

An AI-powered interactive quality control advisor for 10X Genomics single-cell RNA-seq data, built with R Shiny and the Anthropic Claude API.

Upload your CellRanger `.h5` output, explore six QC metrics interactively, and let Claude analyze your data distributions to recommend optimal filtering thresholds — all from a single browser interface.

---

## Features

- **Interactive QC plots** — six metrics visualized with real-time blue/grey highlighting as you drag the threshold sliders
- **Live filtering** — all plots update instantly as thresholds change, showing kept vs filtered cells at every step
- **AI Advisor** — powered by Claude Sonnet via the Anthropic API; analyzes your QC distributions and recommends data-driven filtering thresholds with reasoning
- **Auto-fill sliders** — AI recommendations automatically update the threshold sliders
- **Multi-sample support** — upload and merge multiple `.h5` files in a single session
- **Filtering summary** — before/after cell count table and bar chart comparison
- **Download filtered object** — export the filtered Seurat object as `.rds` for downstream analysis

---

## QC Metrics

| Tab | Metric | Description |
|-----|--------|-------------|
| Cell Counts | Cells per sample | Total cells detected per sample after CellRanger filtering |
| UMI Distribution | nUMI per cell | Number of transcripts captured; low counts suggest empty droplets |
| Genes per Cell | nGene per cell | Unique genes detected; very high counts may indicate doublets |
| UMIs vs Genes | nUMI × nGene scatter | Linear relationship expected; outliers suggest low quality |
| Mitochondrial Ratio | mitoRatio | High mito content (>20%) indicates dying or damaged cells |
| Complexity | log10(Genes/UMI) | Transcriptome complexity; low novelty score suggests poor quality |

---

## AI Advisor

Clicking **Get AI Recommendations** sends your QC summary statistics to Claude Sonnet, which analyzes the distributions and returns:

1. **Recommended thresholds** — `min_nUMI`, `min_nGene`, `max_mitoRatio`, `min_novelty_score`
2. **Reasoning** — data-driven explanation for each threshold
3. **Quality flags** — any concerns about the data quality
4. **Methods text** — a ready-to-use 2–3 sentence QC description for your Methods section

Recommendations are applied automatically to the sliders. You can accept, adjust, or override them before downloading.

---

## Installation

```r
# Install CRAN packages
install.packages(c(
  "shiny", "bslib", "dplyr", "ggplot2",
  "patchwork", "httr", "jsonlite",
  "future", "promises"
))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Seurat", "hdf5r"))
```

---

## Usage

### Local

1. Clone the repo:
```bash
git clone https://github.com/a-sundaresan/RShinyApps-scAdvisorAI.git
cd RShinyApps-scAdvisorAI
```

2. Create `config.json` with your Anthropic API key:
```json
{
  "ANTHROPIC_API_KEY": "sk-ant-your-key-here"
}
```

3. Add `config.json` to `.gitignore`:
```bash
echo "config.json" >> .gitignore
```

4. Run the app:
```r
shiny::runApp("app.R")
```

### Deployed (shinyapps.io)

Set your API key as an environment variable in the shinyapps.io dashboard:

**App → Settings → Environment Variables → `ANTHROPIC_API_KEY`**

The app automatically detects whether to read from `config.json` (local) or the environment variable (deployed).

---

## Input

CellRanger filtered feature barcode matrix in `.h5` format. Multiple samples can be uploaded simultaneously — the app merges them automatically and labels each by filename.

---

## Default Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_nUMI` | 500 | Minimum UMI count per cell |
| `min_nGene` | 300 | Minimum genes detected per cell |
| `max_mitoRatio` | 0.20 | Maximum mitochondrial ratio |
| `min_novelty` | 0.80 | Minimum log10(genes/UMI) score |

Defaults follow published single-cell QC best practices and can be adjusted via sliders or overridden by AI recommendations.

---

## Output

- **Filtered Seurat object** (`.rds`) — ready for downstream clustering and analysis
- **AI-generated Methods text** — copy directly into your manuscript

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `shiny` + `bslib` | UI framework |
| `Seurat` | Single-cell data structure and QC metrics |
| `hdf5r` | Reading 10X `.h5` files |
| `ggplot2` + `patchwork` | Visualization |
| `dplyr` | Data manipulation |
| `httr` + `jsonlite` | Anthropic API calls |
| `future` + `promises` | Async API calls — keeps UI responsive during AI request |

---

## Project Structure

```
RShinyApps-scAdvisorAI/
├── app.R           # Full Shiny application
├── config.json     # API key — local only, NOT committed (add to .gitignore)
├── .gitignore      # Excludes config.json
└── README.md
​```
---

## Related Projects

- [RShinyApps-BulkRNASeqDEAnalysis](https://github.com/a-sundaresan/RShinyApps-BulkRNASeqDEAnalysis) — Interactive bulk RNA-seq differential expression analysis
- [AutoAnnotSC](https://github.com/a-sundaresan/AutoAnnotSC) — Agentic scRNA-seq cell type annotation pipeline

---

## Author

**Aishwarya Sundaresan**
[![Portfolio](https://img.shields.io/badge/Portfolio-a--sundaresan.github.io-black?style=flat-square)](https://a-sundaresan.github.io)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-aishwarya--sundaresan-blue?style=flat-square&logo=linkedin)](https://www.linkedin.com/in/aishwarya-sundaresan/)
