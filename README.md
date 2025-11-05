# Fragmentomics-subtyping
The github repository from research project of BINP39

**Workflow overview**

> 1)Data Availability & Download

> 2)Data preprocessing

> 3)Tumor Fraction analysis

> 4)Feature substraction

> 5)Model training

> 6)Other analysis with R(TF analysis & survival analysis)

For an easier management of the installed programs and their dependencies, all analyses were performed in a Conda environment. Mamba is also needed.

Mamba can be installed with

```bash
conda install -n base -c conda-forge mamba -y
```

**Note: The complete pipelines for (3) Tumor Fraction analysis and (4) Feature subtraction must be run end-to-end on each of the three datasets: LuCaP, the 64-sample lpWGS cohort, and the clinical cohort to derive tumor-fraction estimates and feature matrices for each dataset.**

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
### Griffin (coverage features)
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
* **To run Griffin:**

  The total Griffin workflow consists of 3 parts described in https://github.com/adoebley/Griffin. All of them can be done with snakemake. Get the pipeline with:
  ```bash
  git clone https://github.com/adoebley/Griffin.git
  ```
  For hg38(which is also used in our workflow), there is no need to run part 1: griffin_genome_GC_frequncy.

  For Part 2: griffin_GC_and_mappability_correction

  ```bash
  cd Griffin/snakemakes/griffin_GC_and_mappability_correction/
  ```
  a.Create a samples.yaml with your list of bam files and place it in config/samples.yaml

  b.Edit config.yaml to provide the path to the reference genome (hg38)

  c.Run snakemake with
  ```bash
    snakemake -s griffin_GC_and_mappability_correction.snakefile --cores 48 #Change the cores you use according to your computer
  ```

  For Part 3: griffin_nucleosome_profiling

  ```bash
  cd Griffin/snakemakes/griffin_nucleosome_profiling/
  ```
  a.Copy the samples.GC.yaml from the griffin_GC_correction step(located in Griffin\snakemakes\griffin_GC_and_mappability_correction\results) into the config directory

  b.Make a sites.yaml containing paths to sites you want to focus on to caculate the coverage.

  In this workflow, two kinds of sites are used:

    &emsp;&emsp;(1)AD-ATAC-TF.bed and NE.ATAC.TF.bed which can be generated by
    ```bash
    curl -L -o AD-ATAC-TF.bed https://raw.githubusercontent.com/GavinHaLab/CRPCSubtypingPaper/main/ctdPheno/AD-ATAC-TF.bed
    curl -L -o NE-ATAC-TF.bed https://raw.githubusercontent.com/GavinHaLab/CRPCSubtypingPaper/main/ctdPheno/NE-ATAC-TF.bed
    ```
    &emsp;&emsp;(2)381 transcription factor binding sites which can be generated by
    ```
    curl -sL "https://api.github.com/repos/adoebley/Griffin_analyses/contents/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2?ref=main"
    ```

  c.Edit config.yaml to provide the path to the reference genome (hg38) and other config settings as needed

  d.Run snakemake with
  ```bash
    snakemake -s griffin_nucleosome_profiling.snakefile --cores 48 #Change the cores you use according to your computer
  ```
  
### FinaleToolkit features (WPS, Cleavage profile, Fragment-length, End-motif composition)

* **Environment**
  
  finaletoolkit （0.11.0）

  snakemake （8.26.0）

  bedtools （2.31.1）

  htslib （1.21）

  samtools （1.21）

  The whole environment can be accessed in env/finaletoolkit.yaml in this repository. The environment can be directly built by

    ```bash
    conda env create -f env/finaletoolkit.yaml

    conda activate finaletoolkit_workflow

    pip install finaletoolkit #intall finaltoolkit separately with pip

    ```
* **To run Finaletoolkit:**
  
  a.clone the finaletoolkit workflow repository on your computer with
    
  ```bash
    git clone https://github.com/epifluidlab/finaletoolkit_workflow

    cd finaletoolkit_workflow
  ```
  b.Modify the param.yaml to enter the input and output directory path， set NE-ATAC-TF.sorted.bed or AD-ATAC-TF.sorted.bed as interval_file. Set the features you want to generate. (Here WPS, Cleavage profile, Fragment-length, End-motif composition)

  NE-ATAC-TF.sorted.bed and AD-ATAC-TF.sorted.bed are the two sites file sorted by chromosome name from the to bed files generated in Griffin (Part 3: griffin_nucleosome_profiling). Finaletoolkit workflow can only recognize sorted bed files.
  ```bash
  sort -k1,1V -k2,2n AD-ATAC-TF.bed > AD-ATAC-TF.sorted.bed

  sort -k1,1V -k2,2n NE-ATAC-TF.bed > NE-ATAC-TF.sorted.bed
  ```
  The example param.yaml is in example_config_file/ in this repository.

  c.Run snakemake with
  ```bash
    snakemake --configfile params.yaml --cores 32 --jobs 4 #set different cores and jobs according to your computer
  ```
### Feature extraction and feature matrix generating

* **Environment**

  deeptools（3.5.4）

  bedtools（2.31.1）

  pysam（0.15.3）

  pybigwig（0.3.18）

  py2bit（0.3.0）

  numpy（1.21.6）

  pandas（1.3.5）

  matplotlib（3.5.3）

  scipy（1.7.3）

  python（3.7.12）

  The whole environment can be accessed in env/getFeatures.yaml in this repository. The environment can be directly built by

    ```bash
    conda env create -f env/getFeatures.yaml

    conda activate deeptools

    ```
* **To run getFeatures pipeline:**
 
  a. Clone this repository to your computer with
  ```bash
  git clone https://github.com/QianruiLiu/Fragmentomics-subtyping.git

  cd scripts/getFeatures_snakemake
  ```
  b. Modify the config.yaml, specify the output directory, create a sample list with the output of Finaletoolkit, set the site list (AD-ATAC-TF.sorted.bed AND NE-ATAC-TF.sorted.bed) and the parameters.

  The example config.yaml is in scripts/getFeatures_snakemake/config.yaml in this repository.

  c.Run snakemake with
  ```bash
  snakemake -j 48 #set different cores and jobs according to your computer
  ```
  The ouput of this step is a merged featurematrix "_finaletoolkit_features.tsv".

* **To generate final featurematrix**

  As the whole worflow in this repository should be run on all three types of cohorts, here only take LuCaP data as an example:
  
  a. Merge all the griffin features using scripts/merge_griffin.py

  ```bash
  python scripts/merge_griffin.py -b /Griffin/snakemakes/griffin_nucleosome_profiling/results/   -o /merged/LuCaP_griffin_features.tsv   --only_gccorr
  ```
  b. Merge griffin features with finaltoolkit features (generated by getFeatures pipeline) into one feature matrix with scripts/merge_griff_with_final.py

  ```bash
  python scripts/merge_griff_with_final.py --finaletoolkit merged/LuCaP_finaletoolkit_features.tsv --griffin  merged/LuCaP_griffin_features.tsv --out merged/LuCaP_final_feature_matrix.tsv
  ```
  The output LuCaP_final_feature_matrix.tsv will be used as a training dataset in the machine learning part as an input featurematrix.

* **To generate final Tumor fraction form**
  
  Still take LuCaP data as an example, the script in scripts/collect_ichor_tf.py can be run:
  ```bash
  python scripts/collect_ichor_tf.py --root ~/ichoCNA/scripts/snakemake/results/ichorCNA/ --out /merged/LuCaP_TF.txt
  ```

## 5) Model training
* **Environment and scripts**
  
  The Model training part is done with R. So make sure you have Rscript in your computer. If not, install it with
  ```bash
  sudo apt install -y r-base r-base-dev build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
  ```
  The whole environment can be accessed in env/model_training.yaml in this repository. The environment can be directly built by

    ```bash
    conda env create -f env/model_training.yaml

    conda activate model_training

    cd ~/Fragmentomics_subtyping/model

    ```
    Before training, the following files should be directly in the ~/Fragmentomics_subtyping/model directory:
  
    &emsp;&emsp;LuCaP_final_feature_matrix.tsv, lpWGS_final_feature_matrix.tsv, clinic_final_feature_matrix, healthy_final_features.tsv generated in 4)Feature substraction part;

    &emsp;&emsp;LuCaP_TF.txt, lpWGS_TF.txt, clinic_TF generated by scripts/collect_ichor_tf.py.

    &emsp;&emsp;LuCaP_subtypes.tsv generated from https://github.com/GavinHaLab/CRPCSubtypingPaper/blob/main/Data/LuCaP_subtypes.txt

    Then run the model for prediction. Both lpWGS_final_feature_matrix.tsv and clinic_final_feature_matrix are used as test set.

    ```bash
    Rscript Final_model.R   --train LuCaP_final_feature_matrix.tsv   --test lpWGS_final_feature_matrix.tsv   --labels LuCaP_subtypes.tsv   --tf_train LuCaP_TF.txt   --tf_test lpWGS_TF.txt   --healthy healthy_final_features.tsv   --folds 3 --repeats 3 --grid_size 25 --seed 2025   --gamma 0.5 --arpc_low_tf_q 0.30   --out_prefix results/lpWGS_predicted

    Rscript Final_model.R   --train LuCaP_final_feature_matrix.tsv   --test clinic_final_feature_matrix.tsv   --labels LuCaP_subtypes.tsv   --tf_train LuCaP_TF.txt   --tf_test clinic_TF.txt   --healthy healthy_final_features.tsv   --folds 3 --repeats 3 --grid_size 25 --seed 2025   --gamma 0.5 --arpc_low_tf_q 0.30   --out_prefix results/clinic_predicted
    ```
* **Outputs**
  
  The output will be predicted_combined_LLR_softTF_predictions (All prediction scores), predicted_top10_nepc_prob_samples (Top10 NEPC scores), predicted_cv_auc_per_resample_corrected,
  predicted_cv_pr_auc_per_resample_corrected and two figures (Family-level importance and Coefficient paths of top features)


## 6) Other analysis with R(TF analysis & survival analysis)

### TF analysis

* **Environment and scripts**
  
  The whole environment can be accessed in env/TF_analysis.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/TF_analysis.yaml

  conda activate tfanalysis

  cd ~/Fragmentomics_subtyping/TF_analysis

  ```
  Before running the R script, The following files should be prepared:

  &emsp;&emsp;LuCaP_TF.txt, LuCaP_subtypes.tsv directly in the TF_analysis directory;

  &emsp;&emsp;The griffin results on 381 transcription factor binding sites should be put in a directory. Here take /mnt/f/Project_LuCap/Griffin_results/griffin_nucleosome_profiling/ as an example.

  The script should be run as:
  ```bash
  Rscript TF_analysis.R --subtypes LuCaP_subtypes.tsv --tf_file LuCaP_TF.txt --sites_dir /mnt/f/Project_LuCap/Griffin_results/griffin_nucleosome_profiling/ --out_prefix results/TF_onLucap
  ```
* **Outputs**

  TF_ranking_contrast.tsv, TF_topN_contrast.tsv, TF_ranking_contrast_AR_axis.tsv, TF_ranking_contrast_NE_axis.tsv, TF_topK_discriminability.tsv and four figures
  (TF_top50_contrast.png, TF_top50_bar_by_direction.png, TF_top50_dot_by_direction.png and Top 12 TF × sample z-score heatmap)


### Survival analysis

* **Environment and scripts**
  
   The whole environment can be accessed in env/survival_analysis.yaml in this repository. The environment can be directly built by

  ```bash
  conda env create -f env/survival_analysis.yaml

  conda activate survival

  cd ~/Fragmentomics_subtyping/survival_analysis

  ```
  Before running the R script, The following files should be prepared directly in survival_analysis directory:

  &emsp;&emsp;clinic_predicted_combined_LLR_softTF_predictions.tsv, lpWGS_predicted_combined_LLR_softTF_predictions.tsv which contains prediction score of lpWGS and clinical data generated from 5) Model training.

  &emsp;&emsp;lpWGS_clinical_ubc.txt. A clinical file contains survival period and tumor fraction of 64 lpWGS samples

  &emsp;&emsp;clinic_TF.txt

  The scripts should be run as:
  ```bash
  Rscript Survival_analysis.R --dev_clin lpWGS_clinical_ubc.txt --dev_pred lpWGS_predicted_combined_LLR_softTF_predictions.tsv --clinic_pred clinic_predicted_combined_LLR_softTF_predictions.tsv --clinic_tf clinic_TF.txt --times 12,24,36 --out_prefix results/survival
  ```

* **Outputs**

  survival_clinic_survival_predictions.tsv and two figures (survival_OS_by_pNEPCcal_Top10_lpWGS and survival_OS_by_pNEPCcal_Top10_lpWGS_risktable)

  


  
  
    


























