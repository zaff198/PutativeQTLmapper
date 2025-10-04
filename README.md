# PutativeQTLmapper 🧬

An interactive **Shiny app** for visualizing *putative QTLs* identified from **GWAS** results.  
This tool integrates TASSEL-like outputs (p-values, R², models) with published marker position data to create **genetic (cM)** and **physical (Mb)** QTL distribution maps.

---

## 🧠 Features

- Upload TASSEL GWAS outputs or manually enter data  
- Automatically generate **QTL maps** based on genetic (cM) and physical (Mb) positions  
- Highlight significant loci with `-log10(p)` marker sizing  
- Supports multiple traits and models (GLM / MLM)  
- Download combined high-resolution QTL maps  
- Includes dummy data for demonstration  

---

## 📦 Requirements

You’ll need **R ≥ 4.1** and the following packages:

```r
install.packages(c(
  "shiny", "DT", "ggplot2", "ggrepel",
  "dplyr", "scales", "gridExtra",
  "readxl", "writexl"
))
