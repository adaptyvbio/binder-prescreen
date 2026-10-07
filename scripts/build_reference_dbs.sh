#!/usr/bin/env bash
# Rebuild every prior-art reference database the submission prescreen searches.
#
# Usage:  bash build_reference_dbs.sh OUT_DIR [arm ...]
#   with no arm names, builds all eight public arms, cheapest first, ending with
#   `uspto` - a ~3.4 GB download that unpacks to ~9.3 GB, then a ~25 min build.
#   Budget a couple of hours and ~15 GB of free space for a full run, or skip
#   the patent arm with
#     SKIP_USPTO=1 bash build_reference_dbs.sh OUT_DIR
#   which leaves a ~1.7 GB set that covers everything but patents.
#
# An arm that is already built is skipped, so an interrupted run resumes where it
# stopped; FORCE=1 rebuilds everything.
#
# Needs only `mmseqs`, `curl`, `python3` and network access - no credentials for
# any arm. Point the prescreen at the result with
#   export PRESCREEN_DB_ROOT=OUT_DIR
#
# Each arm lands in OUT_DIR/<arm>/<arm> as an MMseqs2 target database, which is
# the layout prescreen.refdb.resolve() expects. No `createindex` is run: the
# index triples the on-disk size and is a local speed optimisation, not part of
# the shipped reference set. Add it per deployment with
#   mmseqs createindex OUT_DIR/<arm>/<arm> tmp --split-memory-limit 2G
set -euo pipefail

OUT=${1:?usage: build_reference_dbs.sh OUT_DIR [arm ...]}; shift || true
ARMS=("$@")
if [ ${#ARMS[@]} -eq 0 ]; then
  # cheapest first: a full run is dominated by the last two, and the small arms
  # being usable early is worth more than finishing in size order
  ARMS=(thpdb therasabdab plabdab_nano proteinbase_public plabdab pdb swissprot uspto)
  [ -n "${SKIP_USPTO:-}" ] && ARMS=("${ARMS[@]/uspto}")
fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MMSEQS=${MMSEQS:-mmseqs}
PYTHON=${PYTHON:-python3}
export KMP_AFFINITY=disabled   # else MMseqs2 aborts with OMP error #179

mkdir -p "$OUT" "$OUT/.src"
log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
fetch(){ curl -fSL --retry 3 -o "$2" "$1"; }
# For the 3.4 GB patent download: resume a partial transfer rather than starting over,
# falling back to a clean fetch if the server will not serve a range.
fetch_big(){ curl -fSL --retry 3 -C - -o "$2" "$1" || curl -fSL --retry 3 -o "$2" "$1"; }

for arm in "${ARMS[@]}"; do
  [ -z "$arm" ] && continue
  d="$OUT/$arm"
  # `createdb` leaves both of these; their presence is what resolve() looks for
  if [ -z "${FORCE:-}" ] && [ -f "$d/$arm.dbtype" ] && [ -f "$d/$arm.index" ]; then
    log "$arm: already built ($(wc -l < "$d/$arm.index") entries) - skipping, FORCE=1 to rebuild"
    continue
  fi
  mkdir -p "$d"
  if [ "$arm" = uspto ]; then
    log "uspto: ~3.4 GB download (~9.3 GB unpacked) + ~25 min build + ~2.6 GB of database. SKIP_USPTO=1 omits it."
  fi
  log "building $arm"
  case "$arm" in
  pdb)
    # `mmseqs databases` records the upstream release in pdb.version. It keeps
    # every seqres chain, nucleic-acid ones included (~5.7%); they can never
    # match a protein query, so they are left in rather than filtered.
    $MMSEQS databases PDB "$d/pdb" "$OUT/.src/tmp_pdb" ;;
  swissprot)
    $MMSEQS databases UniProtKB/Swiss-Prot "$d/swissprot" "$OUT/.src/tmp_sprot" ;;
  plabdab)
    f="$OUT/.src/paired_sequences.csv"
    fetch https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz "$f.gz"
    gunzip -f "$f.gz"
    # one FASTA record per chain, named plabdab_<ID>, heavy and light of the
    # same entry sharing the ID, so a hit maps back to the entry either way
    $PYTHON - "$f" "$d/plabdab.fasta" <<'PY'
import csv, sys
AA=set("ACDEFGHIKLMNPQRSTVWY"); src,dst=sys.argv[1],sys.argv[2]; n=0
with open(dst,"w") as out:
    for r in csv.DictReader(open(src, encoding="utf-8-sig")):
        t=(r.get("reference_title") or "").strip()
        for col,chain in (("heavy_sequence","heavy"),("light_sequence","light")):
            s=(r.get(col) or "").strip().upper()
            if s and set(s)<=AA:
                out.write(f">plabdab_{r.get('ID','').strip()} {chain}|{t}\n{s}\n"); n+=1
assert n>100000, f"only {n} chains"
print(f"chains={n}")
PY
    $MMSEQS createdb "$d/plabdab.fasta" "$d/plabdab" ;;
  plabdab_nano)
    f="$OUT/.src/plabdab_nano_all.csv"
    fetch https://opig.stats.ox.ac.uk/webapps/plabdab-nano/static/downloads/all_sequences.csv.gz "$f.gz"
    gunzip -f "$f.gz"
    $PYTHON "$HERE/build/plabdab_nano_to_fasta.py" "$f" "$d/plabdab_nano.fasta"
    $MMSEQS createdb "$d/plabdab_nano.fasta" "$d/plabdab_nano" ;;
  therasabdab)
    f="$OUT/.src/therasabdab.csv"
    fetch https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv "$f"
    $PYTHON "$HERE/build/therasabdab_to_fasta.py" "$f" "$d/therasabdab.fasta"
    $MMSEQS createdb "$d/therasabdab.fasta" "$d/therasabdab" ;;
  thpdb)
    fetch https://webs.iiitd.edu.in/raghava/thpdb/sequences/allseq "$d/thpdb.fasta"
    $MMSEQS createdb "$d/thpdb.fasta" "$d/thpdb" ;;
  proteinbase_public)
    # The published designs, from the public API - no credentials. This is the
    # designed-binder prior art: a submission matching it closely is an existing
    # design rather than a new one.
    $PYTHON "$HERE/build/fetch_proteinbase.py" --out "$d/proteins.fasta"
    $MMSEQS createdb "$d/proteins.fasta" "$d/proteinbase_public" ;;
  uspto)
    # EBI's dump is an EMBL-style flat file. createdb has no parser for it and
    # ingests the whole corpus as ONE entry, which is how the deployed patent
    # database came to hold a single record against a 9 GB data file.
    # USPTO_DAT lets you point at an existing unpacked dump instead of re-downloading.
    f=${USPTO_DAT:-"$OUT/.src/uspto_prt.dat"}
    if [ ! -f "$f" ]; then
      fetch_big https://ftp.ebi.ac.uk/pub/databases/patentdata/uspto_prt.dat.gz "$f.gz"
      gunzip -f "$f.gz"
    else
      log "using existing dump $f"
    fi
    $PYTHON "$HERE/build/uspto_to_fasta.py" "$f" "$d/uspto.fasta"
    # --shuffle 0 --compressed 1 keep createdb within ~2 GB RAM; the default
    # shuffle buffers the whole 10.2 M-sequence corpus and OOMs a 16 GB box.
    $MMSEQS createdb "$d/uspto.fasta" "$d/uspto" --shuffle 0 --compressed 1
    rm -f "$d/uspto.fasta" ;;
  *) log "unknown arm: $arm"; exit 2 ;;
  esac
  n=$(wc -l < "$d/$arm.index")
  log "$arm: $n entries"
  # The failure this check exists for is silent: `createdb` reports success on a
  # file it could not parse, and only a real search ever touches the data.
  printf '>probe\nVDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK\n' > "$OUT/.src/probe.fasta"
  $MMSEQS easy-search "$OUT/.src/probe.fasta" "$d/$arm" "$OUT/.src/probe.m8" \
      "$OUT/.src/tmp_probe" --max-seqs 10 -e 1e-3 --split-memory-limit 2G >/dev/null
  log "$arm: searchable"
  rm -rf "$OUT/.src/tmp_probe"
done
log "done; set PRESCREEN_DB_ROOT=$OUT"
