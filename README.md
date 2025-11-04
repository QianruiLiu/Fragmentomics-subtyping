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

* **Environment**

```bash
mamba create -y -n PDX_downloading sra-tools pigz

conda activate PDX_downloading
```

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

* **Environment**

```bash
mamba  create -y -n lpWGS_downloading pyega3

conda activate lpWGS_downloading
```

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
* **Tools and Environment**
  
  bwa (0.7.17)

  samtools (1.14)

  picard (3.4.0)

  gatk (4.6.2.0)

The whole environment can be accessed in env/lucap_preprocessing.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/lucap_preprocessing.yaml

  conda activate pdx_prepoc
  ```
* **Running the pipeline**
  
  The fastq files downloaded in ~/fastq should firstly be aligned with the concatenated hg38+mm10 reference.

  The script concatenate_reference.py is offered by Ha Lab in https://github.com/GavinHaLab/PDX_mouseSubtraction/blob/main/scripts/concatenate_reference.py, can be run as:

  ```bash
  python concatenate_reference.py --humanRef hg38.fa --mouseRef Mus_musculus_NCBI_GRCm38.fa --concatRef ~/fastq/hg38_mm10.fa --tag _mm10
  ```
  Where hg38.fa and Mus_musculus_NCBI_GRCm38.fa can be get by wget:

  ```bash
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/NCBI/GRCm38/Mus_musculus_NCBI_GRCm38.tar.gz
  ```
  Then, do the alignment with bwa mem.

  ```bash
  cd ~/fastq
  mkdir -p ../lucap_original
  for R1 in *_R1.fastq.gz *_1.fastq.gz; do
    [ -e "$R1" ] || continue  
    R2="$(echo "$R1" | sed 's/_R1\.fastq\.gz/_R2.fastq\.gz/; s/_1\.fastq\.gz/_2.fastq\.gz/')" #seek for R2 with same name with R1
    # The output filename
    OUT="../lucap_original/$(basename "$R1" | sed -E 's/(_R1|_1)\.fastq\.gz$//').concat.bam"
    # 1) align with the concatenated reference
    bwa mem -t 16 hg38_mm10.fa "$R1" "$R2" \
      | samtools view -b - \
      | samtools sort -@8 -o "$OUT"
    # 2) build the index
    samtools index "$OUT"
  done
  ```
  Run mouse_substraction snakemake pipeline in https://github.com/GavinHaLab/PDX_mouseSubtraction offered by Ha Lab.
  The pipeline contains mouse_substraction, realignment, Picard MarkDuplicates and base-quality, recalibration (GATK).
  
  To run this pipeline:

    a.clone the mouseSubtraction repository on your computer with
    ```bash
    git clone https://github.com/GavinHaLab/PDX_mouseSubtraction.git

    cd PDX_mouseSubtraction
    ```
    b.Create a samples.yaml with your list of bam files and place it in config/samples.yaml

    c.Add the path of reference, tags, path of VCF files and path of tools in config/config.yaml

    d.Run snakemake with
    ```bash
    snakemake -s subtract_mouse_and_realign.snakefile --cores 48 #Change the cores you use according to your computer
    ```

  The vcf file used in GATK part can be accessed with:
  ```bash
   wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
   wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  ```

### lpWGS
* **Tools and Environment**
  
  bwa-mem2 (2.3)

  samtools(1.22.1)

  picard (3.4.0)

  The whole environment can be accessed in env/lucap_preprocessing.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/lpWGS_preprocessing.yaml

  conda activate 5hmc
  ```
* **Usage of script**
  lpWGS preprocessing can be done with env/align_markdup.sh in this repository.
  ```bash
   bash scripts/align_markdup.sh \
  --ref /refs/hg38.fa \
  --in  ~/ega_lpWGS \
  --out ~/ega_lpWGS/bam_out \
  --threads 28 \
  --remove-dup true #Here you can choose not to remove duplicates
  ```
### Clinical Cohort
* **Tools and Environment**

  fastqc （0.12.1）
  ```bash
  mamba create -n fastqc -c bioconda -c conda-forge fastqc=0.12.1 -y
  mamba activate fastqc
  ```
  ```bash
  mkdir -p ~/clinical_samples/fastqc_out
  for f in ~/clinical_samples/bam_files/*.bam; do
    fastqc -o ~/clinical_samples/fastqc_out -t 8 --noextract "$f"
  done
  ```
    
## 3) Tumor Fraction Estimation (ichorCNA)

* **Environment**

  R (4.2.3)

  BiocManager (1.30.22)

  GenomeInfoDb (1.34.9)

  GenomicRanges (1.48.0)

  ichorCNA (0.5.1)

  HMMcopy(Bioconductor) (1.42.0)

  hmmcopy (0.1.1)

  The whole environment can be accessed in env/lucap_preprocessing.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/ichorCNA.yaml

  conda activate ichor
  ```

Tumor Fraction Estimation is done using ichorCNA on three data cohorts. The snakemake pipeline provided by Ha Lab can be found in https://github.com/GavinHaLab/ichorCNA/tree/v0.4.0/scripts/snakemake

* **To run this pipeline:**
  
  a.clone the mouseSubtraction repository on your computer with
    
  ```bash
    git clone https://github.com/GavinHaLab/ichorCNA.git

    cd ichoCNA/scripts/snakemake/
  ```
  b.Create a samples.yaml with your list of bam files and place it in config/samples.yaml

  c.Add the path of readcounter, path of wig files of hg38, modify genome build and genome style in config/config.yaml.

  Notice: The chromosome naming format should be consistent with the bam files. e.g. chr1, chr2 for UCSC format used in this workflow.

  d.Run snakemake with
  ```bash
    snakemake -s ichorCNA.snakefile --cores 48 #Change the cores you use according to your computer
  ```


## 4) Feature substraction
### Griffin
* **Environment**
  argparse (1.1)

  pysam (0.15.4)
  
  pyBigWig (0.3.17)
  
  pandas (1.3.2)

  numpy (1.21.2)
  
  scipy (1.7.1)
  
  pyyaml (5.3.1)

  matplotlib (3.4.1)

  snakemake (5.5.4)
  
  python (3.7.4)

  The whole environment can be accessed in env/griffin.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/griffin.yaml

  conda activate griffin_demo
  ```

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
