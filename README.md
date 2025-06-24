# Annotation_and_enrichment_for_miRNA-

Here is a clean and structured **English documentation** for the Bash script `Enrich_miRNA.sh` you provided:

---

# 📄 Script Documentation: `Enrich_miRNA.sh`

## 🧾 Purpose

This script automates the **miRNA target gene prediction** and **functional enrichment analysis** based on a list of mature miRNAs. It supports **automatic or manual selection of interaction databases**, allows filtering by confidence scores, and interfaces with an R script (`Annot_genes.R`) to perform enrichment analysis using tools such as **EnrichR**, **ClusterProfiler**, and **ReactomePA**.

---

## 🏁 Input

* A file containing **mature miRNA identifiers**, or manual entry of miRNA names.
* User choice to manually provide interaction **database files** or to automatically download from public sources.
* Optional confidence score thresholds for prediction databases.
* Enrichment analysis method(s) selection.

---

## ⚙️ Workflow Summary

1. **Read miRNA list**:

   * From input file (argument `$1`), or interactively via prompt.

2. **Database setup**:

   * **Manual**: User provides paths to custom database files.
   * **Automatic**: User selects from a list of supported databases:

     * `ENCORI`, `miRDB`, `TargetScan`, `TarBase`, `DianaTarbase`
   * Optional customization of ENCORI parameters (`CLIP-seq`, `Degradome`, `Pan-Cancer`, etc.)
   * Confidence thresholds are requested for `miRDB`, `TargetScan`, and `TarBase`.

3. **Download and extract databases**:

   * Files are stored in a local `DB/` directory.
   * Compressed files are automatically extracted.
   * Clean-up and renaming are handled for compatibility.

4. **Common configuration**:

   * User selects interaction mode: `INTER` (intersection) or `UNION`.
   * Minimum number of interactions to retain a target gene.
   * Output directory name (default: `Results_enrichment`).
   * Enrichment tools: choose from `EnrichR`, `ClusterProfiler`, `ReactomePA`, or all.

5. **Run R script** (`Annot_genes.R`):

   * All configured inputs are passed to the R script for enrichment analysis.

---

## 📥 Example Usage

```bash
./Enrich_miRNA.sh miRNA_list.txt
```

The script will prompt the user interactively for:

* Database source (manual or automatic),
* Database selection,
* Confidence thresholds (if applicable),
* Field mode and thresholds,
* Enrichment method(s).

---

## 🧰 Requirements

* Unix/Linux system with:

  * `bash`
  * `wget`, `gunzip`, `unzip`, `tar`
  * `Rscript` and a properly set up `Annot_genes.R` script
* Internet access if downloading databases

---

## 📦 Output

* A result folder (e.g., `Results_enrichment/`) containing:

  * The list of regulated target genes
  * Results from selected enrichment analysis tools

---

## ✅ Supported Databases

| Name         | Type                     | Source                                                |
| ------------ | ------------------------ | ----------------------------------------------------- |
| ENCORI       | Experimental + Predicted | [ENCORI API](https://rnasysu.com/encori)              |
| miRDB        | Predicted                | [miRDB](https://mirdb.org/)                           |
| TargetScan   | Predicted                | [TargetScan](https://www.targetscan.org/)             |
| TarBase      | Experimental             | [DIANA TarBase](https://dianalab.e-ce.uth.gr)         |
| DianaTarbase | Legacy                   | [DIANA Tools](http://diana.imis.athena-innovation.gr) |

---

## 🔧 Notes

* The script handles basic validation (e.g., numeric input, valid database names).
* Files larger than 100MB should be managed outside GitHub or using Git LFS if versioned.
* Database file names are normalized internally for compatibility with downstream analysis.

---

Let me know if you’d like a **README.md** version of this documentation for your GitHub repository.
