# 🧬 FastSTR

**FastSTR** — Ultra-fast and accurate identification of Short Tandem Repeats (STRs) from long-read DNA sequences. Developed for genome-wide STR detection, consensus construction, and comparative STR analysis.

---

## 📘 Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Command Line Options](#command-line-options)
5. [Input & Output](#input--output)
6. [Usage](#usage)
7. [Performance](#performance)
8. [Citation](#citation)
9. [License](#license)
10. [Changelog](#changelog)


---

## 🌍 Overview

**FastSTR** is a novel and efficient tool for de novo detection of short tandem repeats (STRs) in genomic sequences. It combines fast motif recognition with accurate sequence alignment to achieve both high precision and completeness in STR identification. FastSTR is optimized for large-scale genomic datasets and enables rapid detection of repetitive elements without relying on predefined motif libraries or fixed repeat-length thresholds.

Compared to classical tools like **TRF**, **T-reks**, and **TRASH**, FastSTR achieves:

- ⚡ **High-speed parallel processing** — Processes genomic fragments in parallel, achieving **up to 10× faster runtime**.  
- 🧠 **Context-aware motif recognition** — Uses an **N-gram + Markov** model to identify representative motifs without predefined motif libraries.  
- 🧩 **Segmented global alignment** — Efficiently handles **ultra-long or complex STRs** while maintaining base-level precision.  
- 🔍 **Smart interval merging** — Applies an **interval-gain decision** strategy to accurately resolve overlapping STRs.  
- 🧬 **Enhanced detection in complex regions** — Identifies **confounding or nested repeat regions** (e.g., centromeric satellites) with a novel **density-based concentration test**.  
- 💾 **Lightweight & scalable** — Requires few dependencies, easy to install and run, and supports **multiple operating systems**.  

---

## ⚙️ Installation

### Option 1: Install via `pip`
 
```bash
pip install faststr
```

### Option 2: Install via `conda`
*(coming soon)*  
```bash
conda install -c bioconda faststr
```

### Option 3: Local installation (development)
```bash
git clone https://github.com/yourname/faststr.git
cd faststr
pip install -e .
```

---

## 🚀 Quick Start

### Basic Command

```bash
faststr [--strict | --normal | --loose] [--default] genome.fa
```

### Example

```bash
faststr --strict --default genome.fa
```

This runs FastSTR in **strict mode** using the **default model** to identify STRs in the `genome.fa` file.

---

## ⚡ Command Line Options

| Argument | Type | Default | Description |
|-----------|------|----------|-------------|
| `match` | int | 2 | Match score |
| `mismatch` | int | 5 | Mismatch score |
| `gap_open` | int | 7 | Gap opening penalty |
| `gap_extend` | int | 3 | Gap extension penalty |
| `p_indel` | int | 15 | Indel percentage threshold |
| `p_match` | int | 80 | Match percentage threshold |
| `score` | int | 50 | Alignment score threshold |
| `quality_control` | bool | False | Enable read-level quality control |
| `DNA_file` | str | — | Path to DNA FASTA input |
| `-f` | str | — | Output directory |
| `-s` | int | 1 | Start index |
| `-e` | int | 0 | End index |
| `-l` | int | 15000 | Sub-read length |
| `-o` | int | 1000 | Overlap length |
| `-p` | int | 1 | Number of CPU cores |
| `-b` | float | 0.045 | Motif coverage threshold |

---

## 🧠 Alignment Modes

| Mode | Description |
|------|--------------|
| `--strict` | High precision, recommended for curated assemblies |
| `--normal` | Balanced mode, suitable for most datasets |
| `--loose` | High sensitivity, tolerant of mismatches |

---

## 🧬 Model Presets

| Preset | Description |
|---------|-------------|
| `--default` | Standard scoring model |
| *(future)* `--sensitive` | Optimized for noisy long reads |
| *(future)* `--speed` | Optimized for large-scale detection |

---

## 📥 Input & Output

### Input
- DNA sequences in **FASTA** format

### Output
| File Pattern   | Description |
|----------------|-------------|
| `*detail.dat`  | Contains all STR positions and motifs, quality statistics for each STR, and STR counts per chromosome. |
| `*align.dat`   | Detailed alignment of all STRs against reference STRs, including mismatches and indels. |
| `*.csv`        | Merged STR intervals with representative motifs and summary statistics for each interval. |
| `*.log`        | Processing logs. |

---

## 🧪 Usage

### 1️⃣ Identify STRs in a genome
```bash
faststr --normal --default human_genome.fa
```

### 2️⃣ Use multiple cores
```bash
faststr --strict --default genome.fa -p 8
```

---

## 📈 Performance

| Dataset             | Genome Size | Tool   | Runtime   | Recall | Precision |
|--------------------|------------|--------|-----------|--------|-----------|
| Human (T2T)        | 2.94 G     | TRF    | 18 h 31 min | -      | -         |
|                    |            | FastSTR| 1 h 13 min  | 0.950  | 0.994     |
| Mouse (GRCm39)     | 2.57 G     | TRF    | 1 h 41 min  | -      | -         |
|                    |            | FastSTR| 38 min      | 0.966  | 0.997     |
| Zebrafish (GRCz11) | 1.58 G     | TRF    | 2 h 51 min  | -      | -         |
|                    |            | FastSTR| 25 min      | 0.945  | 0.998     |

*Note: TRF is used as the ground-truth. FastSTR runs based on 72 CPUs.*


---

## 📚 Citation

If you use **FastSTR** in your research, please cite:

> Xingyu Liao *et al.*,  
> **Efficient Identification of Short Tandem Repeats via Context-Aware Motif Discovery and Ultra-Fast Sequence Alignment**,  
> *Nat. Methods*, 2025.

---


## 📄 License

This project is licensed under the **MIT License**.  
See [LICENSE](LICENSE) for more details.

---

## 🧾 Changelog

### v1.0.0 (2025)
- Initial release of FastSTR
- Supports three alignment modes and one default model
- Implemented parallel computation
- Added `.csv`, `.dat`, `.log` outputs
---





