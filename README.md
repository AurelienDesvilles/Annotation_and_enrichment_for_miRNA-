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

Here is a detailed **English documentation** for your R script `Annot_genes.R`:

---

# 📄 Script Documentation: `Annot_genes.R`

## 🎯 Purpose

This script performs **target gene annotation** and **functional enrichment analysis** for a given list of **mature miRNAs**, using one or more miRNA-target databases. It supports filtering by interaction count and field mode (`INTER` or `UNION`) and executes downstream enrichment analysis using **EnrichR**, **ClusterProfiler**, and **ReactomePA**.

---

## 🧾 Input Arguments

The script must be run via command line and accepts **six positional arguments**:

```bash
Rscript Annot_genes.R <miRNA_list> <databases> <output_directory> <field_mode> <nb_interact> <enrichment_databases>
```

### 🔢 Arguments Explained

| Position | Argument               | Description                                                                                            |
| -------- | ---------------------- | ------------------------------------------------------------------------------------------------------ |
| 1        | `miRNA_list`           | Comma-separated list of mature miRNA IDs                                                               |
| 2        | `databases`            | Comma-separated paths to database files                                                                |
| 3        | `output_directory`     | Output folder for results                                                                              |
| 4        | `field_mode`           | Mode for intersecting results: `"INTER"` or `"UNION"`                                                  |
| 5        | `nb_interact`          | Minimum number of distinct miRNAs regulating a gene                                                    |
| 6        | `enrichment_databases` | Comma-separated list of enrichment tools: `EnrichR`, `ClusterProfiler`, `ReactomePA`, or a combination |

---

## 📚 Database Handling

Each selected database is parsed and filtered based on format and specific thresholds. Supported databases include:

| Database     | Filtering Method                                                                  |
| ------------ | --------------------------------------------------------------------------------- |
| ENCORI       | Direct import (optional filtering commented)                                      |
| miRDB        | Score threshold (read from `_score.txt` file); gene names converted via `biomaRt` |
| TargetScan   | Filter by PCT score and species; miRNA family parsing                             |
| TarBase      | Filter by `microt_score`                                                          |
| DianaTarbase | Filter by `species == Homo sapiens`                                               |
| miRTarBase   | Filter by support type `"Functional MTI"`                                         |

Each database is transformed into a common format: pairs of `miRNA` and `geneName`.

---

## 🔗 Interaction Aggregation

* All interactions across databases are merged.
* If `field_mode` is `INTER`, only miRNAs present in all databases are retained.
* Genes targeted by at least `nb_interact` unique miRNAs are selected.

---

## 📄 Output Files

The following result files are generated in the `output_directory`:

| File                              | Description                         |
| --------------------------------- | ----------------------------------- |
| `Selected_Genes.txt`              | Final list of selected gene symbols |
| `Interactions_miRNA_filtered.txt` | miRNA-to-gene filtered interactions |
| Enrichment result files           | See below                           |

---

## 🧬 Functional Enrichment Analysis

For the selected genes, the script performs enrichment analysis using the tools specified:

### 🔹 EnrichR

* Databases used: `GO_Biological_Process_2021`, `KEGG_2021_Human`, `Reactome_2022`
* Output:

  * `GO_Results.txt`
  * `KEGG_Results.txt`
  * `Reactome_Results.txt`

### 🔹 ClusterProfiler

* Gene Ontology Biological Process analysis
* Visual output: `ClusterProfiler.pdf` (barplot of top 10 terms)
* Data: `ClusterProfiler_results.csv`

### 🔹 ReactomePA

* Reactome pathway enrichment
* Visual output: `ReactomePA.pdf`
* Data: `ReactomePA_results.csv`

> All visual outputs are saved as PDF using the **Cairo** graphics library.

---

## 💻 Requirements

### R Packages (loaded silently)

* `dplyr`, `readr`, `tidyr`, `purrr`
* `biomaRt`, `org.Hs.eg.db`, `DOSE`
* `enrichR`, `clusterProfiler`, `ReactomePA`
* `Cairo` (for PDF plots)

Make sure they are installed via:

```r
install.packages("BiocManager")
BiocManager::install(c("biomaRt", "org.Hs.eg.db", "DOSE", "clusterProfiler", "ReactomePA", "enrichR"))
install.packages("Cairo")
```

---

## ✅ Example

```bash
Rscript Annot_genes.R "hsa-miR-21,hsa-miR-34a" "DB/ENCORI.txt,DB/miRDB.txt" Results_enrichment UNION 2 "EnrichR,ClusterProfiler"
```

---

## 📝 Notes

* Gene symbols are internally converted to Entrez IDs when required (e.g., for ReactomePA).
* Intermediate and final filtering ensures only biologically meaningful genes are retained.
* If any enrichment tool fails to return results, the script notifies the user and skips plot generation.

---



