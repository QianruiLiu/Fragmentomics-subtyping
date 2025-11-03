# Fragmentomics-subtyping
The github repository from research project of BINP39

**Workflow overview**

> 1)Data Availability & Download

> 2)Data preprocessing

> 3)Tumor Fraction analysis

> 4)Feature substraction

> 5)Model training

> 6)Analysis with R(TF analysis & survival analysis)

For an easier management of the installed programs and their dependencies, all analyses were performed in a Conda environment. Mamba is also needed.

Mamba can be installed with

```bash
conda install -n base -c conda-forge mamba -y
```

## 1) Data Availability & Download

### LuCaP PDX cfDNA

* **Source:** NCBI (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA900550)

  The accession list can be get from the website as SraAccList.csv.

* **Download:**

```bash
mkdir -p ~/sra_raw && cd ~/sra_raw #Put SraAccList.csv file in this directory at the same time

prefetch --option-file SraAccList.csv --progress 

mkdir -p ~/fastq && cd ~/fastq 

cat ~/sra_raw/SraAccList.csv | xargs -n1 -P4 -I{} \
  fasterq-dump {} -O ~/fastq -t ~/tmp -e 8 -p --split-files  #Get fastq files

find ~/fastq -name "*.fastq" -print0 | xargs -0 -n1 -P8 pigz -p 8   #Compressing fastq files

```

### 64 lpWGS samples

* **Source:** EGA dataset (https://ega-archive.org/datasets/EGAD00001008462) (granted access required).
* **Download (pyEGA3):**

```bash
mkdir -p ~/ega_lpWGS

pyega3 -cf credential_file.json fetch EGAD00001008462 --output-dir ~/ega_lpWGS
#The example of creating a credential file can be seen at https://github.com/EGA-archive/ega-download-client/blob/master/pyega3/config/default_credential_file.json
```

### Clinical Cohort 

* Place raw FASTQs under `data/clinical_fastq/` (or BAMs under `results/bam/`). Provide a `samples/clinical.tsv`.

## 2) Data preprocessing

### LuCaP PDX cfDNA


### 3.4 Tumor Fraction Estimation (ichorCNA)

* Run ichorCNA on BAMs; ensure `*.params.txt` produced. Then:

```bash
python scripts/ichor_tf_collect.py --root results/ichor --out tables/ichor_tf.tsv
```

---

## 4) Software Stack

* **Snakemake:** >=7.x (pinned in `envs/*`)
* **Conda/Mamba:** for reproducible envs
* **BWA‑MEM2 / Samtools / Picard:** alignment + duplicates + indexing
* **FinaleToolkit:** WPS, cleavage, fragment length, motif features (external repo)
* **Griffin:** nucleosome/coverage profiles over AD‑ATAC‑TF and NE‑ATAC‑TF loci
* **R (>=4.2) + packages:** tidymodels, glmnet, survival, survminer, ggplot2, etc.

Environment creation:

```bash
mamba env create -f envs/finaletoolkit.yml
mamba env create -f envs/griffin.yml
mamba env create -f envs/rstats.yml
```

---

## 5) Pipeline A — FinaleToolkit (External Snakemake)

> This is maintained in a separate repository. We call it **Pipeline A**.

1. **Clone external repo** (ideally as a sibling directory):

   ```bash
   git clone https://github.com/epifluidlab/finaletoolkit_workflow ../finaletoolkit_workflow
   ```
2. **Configure** its `config.yaml` per that repo’s instructions (references, loci, output paths).
3. **Run:**

   ```bash
   cd ../finaletoolkit_workflow
   snakemake -c 16 --use-conda
   ```
4. **Outputs (used downstream in this repo):**

   * Per‑sample **WPS** and **cleavage** tidy profiles (AD/NE)
   * **Fragment‑length** summaries/histograms
   * **End‑motif** fractions and **motif entropy**
5. **Optional post‑processing:** Convert tidy meta‑profiles into metrics

   ```bash
   # Example: WPS (AD loci) central_mean, shoulder_mean, contrast, amplitude, phasing_strength
   python ../<this-repo>/scripts/profile_metrics.py \
     results/<sample>/wps/AD/tidy.tsv \
     ../<this-repo>/tables/<sample>.wps.AD --binsize 5
   ```

---

## 6) Pipeline B — This Repo’s Snakemake (Griffin + Integrations)

> **Pipeline B** orchestrates: alignment/QC → Griffin meta‑profiles → site‑wise / family‑wise feature tables → joins with FinaleToolkit outputs.

### 6.1 Alignment & Duplicates

```bash
bash scripts/align_markdup.sh \
  --fastq1 data/fastq/<S>_R1.fastq.gz \
  --fastq2 data/fastq/<S>_R2.fastq.gz \
  --out results/bam/<S>.bam \
  --ref /path/to/hg38.fa \
  --remove-dup true
```

* Produces indexed BAMs with RG tags; optional duplicate removal.

### 6.2 Griffin Nucleosome/Coverage Profiles

* Run Griffin at **AD‑ATAC‑TF** and **NE‑ATAC‑TF** loci; export per‑sample **central_mean / shoulder_mean / contrast / amplitude** and optionally **GC‑corrected** counterparts.
* **Site‑wise expansion:**

```bash
python scripts/griffin_sitewise.py \
  --input results/griffin/<S>/profiles.tsv \
  --output tables/griffin_sitewise/<S>.wide.tsv \
  --gc-correct true
```

### 6.3 Feature Family Merges

```bash
python scripts/merge_features.py \
  --ad tables/ad_features.tsv \
  --ne tables/ne_features.tsv \
  --out tables/family_merged.tsv
```

### 6.4 Join FinaleToolkit + Griffin

```bash
python scripts/join_final_griffin.py \
  tables/finaletoolkit_features.tsv \
  tables/griffin_features.tsv \
  -o tables/merged_features.tsv
```

**Make it one command:**

```bash
snakemake -s workflow/Snakefile -c 16 --use-conda
```

---

## 7) Analyses (R) — Each Script Generates Visualizations

All three scripts accept `--in`/`--out` style CLI args (see headers in each script). Typical defaults assume:

* `tables/merged_features.tsv` (features + labels + ctF)
* `figs/` for plots, `tables/` for exported TSVs

### 7.1 `R/01_ml_classifier.R`

* **Goal:** Train elastic‑net classifier (ARPC vs NEPC), CV with 1‑SE model selection, export probabilities.
* **ctF handling:** optional soft‑correction and/or ctF‑aware calibration (see 7.3 for survival usage).
* **Run:**

  ```bash
  Rscript R/01_ml_classifier.R \
    --features tables/merged_features.tsv \
    --out_dir results/ml \
    --fig_dir figs/ml
  ```
* **Outputs:** `tables/predictions.tsv` (columns: sample_id, .pred_ARPC, .pred_NEPC, ctF, etc.), ROC/PR curves, coefficient paths, feature family VIPs.

### 7.2 `R/02_tf_analysis.R`

* **Goal:** Rank TFs by signed AUC (sAUC_dir), generate top‑k barplots and sample×TF heatmaps; optionally export TF‑centric profile snapshots.
* **Run:**

  ```bash
  Rscript R/02_tf_analysis.R \
    --griffin tables/griffin_sitewise/ \
    --labels tables/labels.tsv \
    --out_dir results/tf \
    --fig_dir figs/tf
  ```
* **Outputs:** `figs/tf_rank_top50.pdf`, `figs/tf_heatmap_top12.pdf`, `tables/tf_rank.tsv`.

### 7.3 `R/03_survival_projection.R`

* **Goal:** Calibrate NEPC probability vs ctF (regress logit(p) on ctF → inverse‑logit), fit Cox PH on development cohort, and **project absolute S(t)** at 12/24/36 months for the clinical cohort.
* **Run:**

  ```bash
  Rscript R/03_survival_projection.R \
    --pred tables/predictions.tsv \
    --clin tables/clinical.tsv \
    --out_dir results/surv \
    --fig_dir figs/surv
  ```
* **Outputs:** `figs/km_curves.pdf`, `tables/surv_projection.tsv` (per‑sample S_12m, S_24m, S_36m).

---

## 8) Configuration Contracts

`config/config.yaml` keys (suggested):

```yaml
references:
  fasta: /path/to/hg38.fa
  blacklist: /path/blacklist.bed
  ad_atac_bed: loci/AD-ATAC-TF.bed
  ne_atac_bed: loci/NE-ATAC-TF.bed

paths:
  raw_fastq: data/fastq
  bam_out:   results/bam
  griffin_results: results/griffin
  finaletoolkit_results: ../finaletoolkit_workflow/results
  tables:    tables
  figs:      figs

cohorts:
  lucap:     config/samples/lucap.tsv
  lpwgs64:   config/samples/lpwgs64.tsv
  clinical:  config/samples/clinical.tsv

ega:
  dataset: EGAD00001008462
  credential_file: ~/.ega/credential_file.json

threads: 16
```

**Sample sheet contract** (`config/samples/lucap.tsv`):

```
sample_id	group	fastq1	fastq2
LuCaP_35CR	ARPC	data/fastq/LuCaP_35CR_R1.fq.gz	data/fastq/LuCaP_35CR_R2.fq.gz
...
```

---

## 9) Reproducibility & Provenance

* Pin **conda envs** in `envs/` and record exact versions (`snakemake --version`, `mamba list` → `logs/versions.txt`).
* Record **commit SHAs** for external repos (FinaleToolkit workflow, Griffin).
* Export **frozen inputs** (sample sheets, loci BEDs, reference FASTAs checksums) into `docs/provenance/`.
* Keep **random seeds** fixed in R scripts; log CV splits.

---

## 10) Troubleshooting & FAQs

* **pyEGA3 on Windows?** Use **WSL** (Ubuntu), ensure Java installed, store `credential_file.json` under WSL home (`/home/<user>/.ega/`).
* **FinaleToolkit asks for BED when doing genome‑wide?** Provide **locus panels** explicitly (AD/NE ATAC‑TF) to ensure biological contrast; genome‑wide WPS can be large and memory intensive—use tiling or predefined panels.
* **ctF bias in probabilities?** Use the calibration step (logit regression on ctF) before survival modeling.
* **Some Griffin profiles show cross‑signals on AD vs NE panels?** A few borderline samples may display cross‑signals; rely on rank‑based sAUC_dir and calibrated probabilities rather than a hard 0.5 threshold.

---

## 11) Citations

* FinaleToolkit & workflow repo
* Griffin nucleosome profiling
* ichorCNA for ctDNA fraction
* Any public cohorts (SRA/EGA) used
* Reference: https://doi.org/10.1101/2022.06.21.496879

---

## 12) Changelog

* `YYYY‑MM‑DD` — Initial public release; added Snakemake B, linked FinaleToolkit A; uploaded R scripts.
* `YYYY‑MM‑DD` — Added survival projection; updated config schema; improved TF visualizations.

---

### How to Customize This Template

* Replace all angle‑bracket placeholders (`<...>`) with your actual paths, URLs, and sample sheets.
* If your R scripts accept different CLI flags, update Section 7 accordingly.
* If you keep FinaleToolkit as a submodule, add instructions (`git submodule add ...`).

> This README is designed for **printability** (for supervisors) and **hands‑on reproducibility** (for labmates). Keep the **file contracts** in Sections 6–8 stable to avoid breaking downstream steps.
