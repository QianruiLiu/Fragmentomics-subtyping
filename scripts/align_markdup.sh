#!/usr/bin/env bash
set -euo pipefail

# =========================
# align_hg38_simple.sh — single input dir
# =========================
# Requirements: bash, awk, sed, coreutils, bwa-mem2 (or bwa), samtools, picard
#
# Supported filename patterns (any of):
#   1) <SM>__<ID>_R{1,2}.fastq.gz
#   2) <SM>_<ANY>__<ID>_R{1,2}.fastq.gz   # still compatible if a "type" token exists
#   3) <SM>_R{1,2}.fastq.gz               # when no ID, a hash will be used as ID
#
# RG inference:
#   - SM = the part left of "__" in the filename stem (or the whole stem if no "__")
#   - ID = the part right of "__" (or a hash of the stem if no "__")
#   - LB = "lib"
#   - PU = <ID>.1
#
# Example:
#   bash align_hg38_simple.sh --ref /refs/hg38.fa --in /data/fastq --out /data/bam_out --threads 28 --remove-dup false --aligner bwa-mem2

# -------- default parameters --------
REF=""
IN_DIR=""
OUT_ROOT=""
THREADS=16
REMOVE_DUP="false"         # true=hard dedup; false=mark duplicates only
ALIGNER="bwa-mem2"         # or "bwa"

print_usage() {
  cat <<'EOF'
Usage: align_hg38_simple.sh --ref REF.fa --in FASTQ_DIR --out OUT_DIR [options]

Required:
  --ref PATH        hg38 reference FASTA (must be indexed for bwa; samtools faidx + .dict recommended)
  --in  PATH        Input FASTQ directory (contains paired *_R1.fastq.gz / *_R2.fastq.gz)
  --out PATH        Output root directory

Options:
  --threads INT     Number of threads (default: 16)
  --remove-dup [true|false]
                    Picard REMOVE_DUPLICATES (default: false)
  --aligner [bwa-mem2|bwa]
                    Aligner (default: bwa-mem2)

Outputs:
  <OUT_DIR>/{logs,tmp}/
  <OUT_DIR>/<SM>.<ID>.sorted.bam (+ .bai)
  <OUT_DIR>/<SM>.<ID>.markdup.bam or .dedup.bam (+ .bai)
  <OUT_DIR>/<SM>.<ID>.markdup.metrics.txt
EOF
}

# -------- parse args --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref)        REF="${2:-}"; shift 2;;
    --in)         IN_DIR="${2:-}"; shift 2;;
    --out)        OUT_ROOT="${2:-}"; shift 2;;
    --threads)    THREADS="${2:-}"; shift 2;;
    --remove-dup) REMOVE_DUP="${2:-}"; shift 2;;
    --aligner)    ALIGNER="${2:-}"; shift 2;;
    -h|--help)    print_usage; exit 0;;
    *) echo "Unknown option: $1"; print_usage; exit 1;;
  esac
done

# -------- validations --------
[[ -z "$REF" || -z "$IN_DIR" || -z "$OUT_ROOT" ]] && { echo "[ERR] --ref/--in/--out are required."; print_usage; exit 1; }
[[ -d "$IN_DIR" ]] || { echo "[ERR] Input directory not found: $IN_DIR"; exit 1; }

if [[ "$ALIGNER" == "bwa-mem2" ]]; then
  command -v bwa-mem2 >/dev/null || { echo "[ERR] bwa-mem2 not found in PATH"; exit 1; }
elif [[ "$ALIGNER" == "bwa" ]]; then
  command -v bwa >/dev/null || { echo "[ERR] bwa not found in PATH"; exit 1; }
else
  echo "[ERR] --aligner must be bwa-mem2 or bwa"; exit 1
fi
command -v samtools >/dev/null || { echo "[ERR] samtools not found in PATH"; exit 1; }
command -v picard   >/dev/null || { echo "[ERR] picard not found in PATH"; exit 1; }
[[ -f "$REF" ]] || { echo "[ERR] Reference FASTA not found: $REF"; exit 1; }

# Index hints (not enforced)
[[ -f "${REF}.fai" ]] || echo "[WARN] Missing ${REF}.fai; consider: samtools faidx $REF"
if [[ "$ALIGNER" == "bwa-mem2" ]]; then
  ls "${REF}".* 2>/dev/null | grep -Eq '\.(0123|amb|ann|bwt\.[^/]+|pac)$' || \
    echo "[WARN] bwa-mem2 index not detected; consider: bwa-mem2 index $REF"
else
  ls "${REF}".* 2>/dev/null | grep -Eq '\.(amb|ann|bwt|pac|sa)$' || \
    echo "[WARN] bwa index not detected; consider: bwa index $REF"
fi

# -------- directories --------
mkdir -p "$OUT_ROOT"/{logs,tmp}
LOG_DIR="${OUT_ROOT}/logs"

log() { printf "[%(%F %T)T] %s\n" -1 "$*" >&2; }

shopt -s nullglob
# Support both *_R1.fastq.gz and *_1.fastq.gz
for R1 in "${IN_DIR}"/*_R1.fastq.gz "${IN_DIR}"/*_1.fastq.gz; do
  [[ -e "$R1" ]] || continue
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
  [[ "$R1" == "$R2" ]] && R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  if [[ ! -f "$R2" ]]; then
    log "[WARN] Missing pair for $R1 — skipped"
    continue
  fi

  fname="$(basename "$R1")"
  stem="${fname%_R1.fastq.gz}"
  stem="${stem%_1.fastq.gz}"

  # Support <SM>__<ID> or <SM>_<ANY>__<ID>; otherwise fallback to hash
  SM="$stem"; ID=""
  if [[ "$stem" == *"__"* ]]; then
    left="${stem%%__*}"   # left of "__"
    right="${stem##*__}"  # right of "__"
    SM="$left"
    ID="$right"
  else
    ID="$(echo -n "$stem" | md5sum | cut -c1-8)"
    log "[INFO] Filename has no '__'; fallback RG: SM=${SM}, ID=${ID}"
  fi

  RG="@RG\tID:${ID}\tSM:${SM}\tLB:lib\tPL:ILLUMINA\tPU:${ID}.1"
  OUT_PREFIX="${OUT_ROOT}/${SM}.${ID}"
  SORT_BAM="${OUT_PREFIX}.sorted.bam"
  METRICS="${OUT_PREFIX}.markdup.metrics.txt"
  if [[ "$REMOVE_DUP" == "true" ]]; then
    MD_BAM="${OUT_PREFIX}.dedup.bam"
  else
    MD_BAM="${OUT_PREFIX}.markdup.bam"
  fi
  LOGFILE="${LOG_DIR}/${SM}.${ID}.log"

  log "[RUN] ${SM}.${ID}"
  {
    set -x
    if [[ "$ALIGNER" == "bwa-mem2" ]]; then
      bwa-mem2 mem -t "${THREADS}" -R "${RG}" "${REF}" "${R1}" "${R2}" \
        | samtools view -@ "${THREADS}" -b - \
        | samtools sort -@ "${THREADS}" -o "${SORT_BAM}" -
    else
      bwa mem -t "${THREADS}" -R "${RG}" "${REF}" "${R1}" "${R2}" \
        | samtools view -@ "${THREADS}" -b - \
        | samtools sort -@ "${THREADS}" -o "${SORT_BAM}" -
    fi
    samtools index -@ "${THREADS}" "${SORT_BAM}"

    if [[ "$REMOVE_DUP" == "true" ]]; then
      picard MarkDuplicates I="${SORT_BAM}" O="${MD_BAM}" M="${METRICS}" \
        REMOVE_DUPLICATES=true CREATE_INDEX=true
    else
      picard MarkDuplicates I="${SORT_BAM}" O="${MD_BAM}" M="${METRICS}" \
        REMOVE_DUPLICATES=false CREATE_INDEX=true
    fi
    set +x
  } >"${LOGFILE}" 2>&1 || { log "[ERR] Failed: ${SM}.${ID} (see ${LOGFILE})"; continue; }

  log "[OK ] ${SM}.${ID}"
done
shopt -u nullglob

log "All done."
