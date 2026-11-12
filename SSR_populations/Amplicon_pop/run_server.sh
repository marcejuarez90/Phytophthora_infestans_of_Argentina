#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════════
#  run_server.sh — full amplicon phylogeny, to run on the server
#
#  Both analyses (all samples; Argentine + European reference set), without the
#  samples in samples_exclude.txt, with
#    * NJ tree, 1000 bootstrap replicates          -> <analysis>_nj.nwk
#    * per-sample QC and genotype-label check      -> <analysis>_sample_qc.tsv
#
#  Usage (from this folder):
#      bash run_server.sh [VCF] [OUTDIR]
#  Defaults: VCF=/u03/mjuarez/marce/Pi_pop/total/totalpop_final_filter.vcf
#            OUTDIR=results_excl_08_7054A
#  Environment variables: NBOOT (default 1000).
#
#  Everything printed is also saved in <OUTDIR>/run.log (with the versions used).
# ══════════════════════════════════════════════════════════════════════════════
cd "$(dirname "$0")" || exit 1

VCF="${1:-/u03/mjuarez/marce/Pi_pop/total/totalpop_final_filter.vcf}"
OUT="${2:-results_excl_08_7054A}"
NBOOT="${NBOOT:-1000}"

PY=$(command -v python3 || command -v python)
[ -z "$PY" ] && { echo "ERROR: python not found"; exit 1; }
"$PY" -c "import numpy" 2>/dev/null || { echo "ERROR: numpy is not available for $PY"; exit 1; }
[ -f "$VCF" ] || { echo "ERROR: VCF not found: $VCF"; exit 1; }

mkdir -p "$OUT"
{
  echo "== run_server.sh  $(date)  host: $(hostname)"
  echo "python: $("$PY" --version 2>&1)  numpy: $("$PY" -c 'import numpy; print(numpy.__version__)')"
  echo "vcf: $VCF   md5: $(md5sum "$VCF" | cut -d' ' -f1)"
  echo "excluded: $(grep -v '^#' samples_exclude.txt | tr '\n' ' ')"
  echo
} | tee "$OUT/run.log"

"$PY" amplicon_nj_phylogeny.py \
    --vcf "$VCF" \
    --outdir "$OUT" \
    --analyses all,ar_eu \
    --exclude samples_exclude.txt \
    --genotypes sample_genotypes.tsv \
    --nboot "$NBOOT" \
    --max-site-missing 0.2 \
    --max-sample-missing 0.5 2>&1 | tee -a "$OUT/run.log"

echo
echo "Done. Bring back (or just re-download this folder):"
echo "  $OUT/ar_eu/  and  $OUT/all/   (*_nj.nwk, *_sample_qc.tsv, *_report.txt, run.log)"
