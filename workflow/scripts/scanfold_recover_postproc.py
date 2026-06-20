#!/usr/bin/env python3
"""
Recovery script for interrupted ScanFold2.0 post-processing.

Use when the scan CSV and FinalPartners.txt are complete but the output files
(CT, DBN, WIG, BP) were never written because the run was interrupted.

Skips the expensive RNAfold scanning step by reading the existing CSV,
then re-runs the fold step and writes all output files.

Usage:
    conda activate scanfold2
    python workflow/scripts/scanfold_recover_postproc.py \
        --fasta results/.../scanfold_retry_batches/batch_NNN.fa \
        --csv   results/.../scanfold2_gpu/batch_NNN/SEQID.win_120.stp_1.csv \
        --seq   SEQID \
        --outdir results/.../scanfold2_gpu/batch_NNN
"""

import argparse
import os
import re
import statistics
import sys
import time

import numpy as np
import pandas as pd
from Bio import SeqIO

start_time = time.time()

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True, help="Batch FASTA file")
parser.add_argument(
    "--csv", required=True, help="Existing scan CSV for the target sequence"
)
parser.add_argument(
    "--seq", required=True, help="Sequence ID to recover (e.g. ENST00000412264.2)"
)
parser.add_argument(
    "--outdir", required=True, help="Output directory (same as original batch dir)"
)
parser.add_argument("--step", type=int, default=1)
parser.add_argument("--window", type=int, default=120)
parser.add_argument("--temp", type=int, default=37)
args = parser.parse_args()

fasta_file = os.path.abspath(args.fasta)
csv_file = os.path.abspath(args.csv)
seq_id = args.seq
output_dir = os.path.abspath(args.outdir)
step_size = args.step
window_size = args.window
temperature = args.temp

# ---------------------------------------------------------------------------
# Bootstrap ScanFoldFunctions from the software directory
# ---------------------------------------------------------------------------

script_dir = os.path.dirname(os.path.abspath(__file__))
scanfold_dir = os.path.join(script_dir, "..", "..", "software", "ScanFold2.0")
scanfold_dir = os.path.normpath(scanfold_dir)
sys.path.insert(0, scanfold_dir)

from ScanFoldFunctions import (
    NucPair,
    NucZscore,
    best_basepair,
    competing_pairs,
    makedbn,
    write_bp,
    write_ct,
    write_wig,
    write_wig_dict,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def elapsed():
    return round(time.time() - start_time, 1)


def log(msg):
    print(f"[{elapsed()}s] {msg}", flush=True)


def get_path(filename):
    return os.path.join(output_dir, filename)


# ---------------------------------------------------------------------------
# Determine output file base names (mirrors ScanFold2.0.py naming logic)
# ---------------------------------------------------------------------------

basename = os.path.basename(fasta_file).split(".")[0]  # e.g. batch_00002_retry_00000

# Detect multi-record input (retry batches always are)
with open(fasta_file) as f:
    n_records = sum(1 for line in f if line.startswith(">"))
multi_record = n_records > 1

structure_basename = f"{basename}.{seq_id}" if multi_record else basename
no_filter_base = f"{structure_basename}.no_filter"
minus_1_base = f"{structure_basename}.minus_1"
minus_2_base = f"{structure_basename}.minus_2"
outname = f"{seq_id}.win_{window_size}.stp_{step_size}.csv"

log(f"Recovering: {seq_id}")
log(f"Output dir: {output_dir}")
log(f"structure_basename: {structure_basename}")

# ---------------------------------------------------------------------------
# Load sequence from FASTA
# ---------------------------------------------------------------------------

seq = None
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id == seq_id or record.name == seq_id:
        seq = record.seq.transcribe()
        break

if seq is None:
    sys.exit(f"ERROR: sequence '{seq_id}' not found in {fasta_file}")

cur_record_length = len(seq)
log(f"Sequence length: {cur_record_length} nt")

# Build nuc_dict (coordinate → NucZscore)
nuc_dict = {}
for i, nuc in enumerate(seq, start=1):
    x = NucZscore(nuc, i)
    nuc_dict[x.coordinate] = x

start_coordinate = str(list(nuc_dict.keys())[0])
end_coordinate = str(list(nuc_dict.keys())[-1])

# ---------------------------------------------------------------------------
# Load existing scan CSV
# ---------------------------------------------------------------------------

log(f"Reading CSV: {csv_file}")
df = pd.read_csv(csv_file, sep="\t", index_col=0)

# The CSV has a typo column "Sequeunce" (empty) and the real "Sequence"
sequences = df["Sequence"].values
structures = df["Structure"].values
starts = df["Start"].values.astype(int)
zscores = df["Z-score"].values
mfes = df["NativeMFE"].values
eds = df["ED"].values

mfe_list = mfes.tolist()
zscore_list = zscores.tolist()
ed_list = eds.tolist()
minz = df["Z-score"].min()
log(f"Windows loaded: {len(df)}, min Z-score: {minz}")

# ---------------------------------------------------------------------------
# ScanFold-Fold: build bp_dict from existing per-window data
# ---------------------------------------------------------------------------

log("Building bp_dict from scan data...")
bp_dict = {}

for idx in range(len(df)):
    sequence = sequences[idx]
    structure = structures[idx]
    start = starts[idx]
    zscore = zscores[idx]
    mfe = mfes[idx]
    ed = eds[idx]

    unpaired_positions = np.array([i for i, s in enumerate(structure) if s == "."])
    for pos in unpaired_positions:
        nucleotide = sequence[pos]
        coordinate = pos + start
        x = NucPair(nucleotide, coordinate, nucleotide, coordinate, zscore, mfe, ed)
        bp_dict.setdefault(coordinate, []).append(x)

    base_pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i + 1)
        elif char == ")":
            if stack:
                left_pos = stack.pop()
                base_pairs.extend([left_pos, i + 1])

    for l in range(0, len(base_pairs), 2):
        lbp = base_pairs[l]
        rbp = base_pairs[l + 1]
        lb = sequence[lbp - 1]
        rb = sequence[rbp - 1]
        lbp_coord = lbp + start - 1
        rbp_coord = rbp + start - 1
        x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, mfe, ed)
        z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, mfe, ed)
        bp_dict.setdefault(lbp_coord, []).append(x)
        bp_dict.setdefault(rbp_coord, []).append(z)

log(f"bp_dict built ({len(bp_dict)} entries). Computing best base pairs...")

# ---------------------------------------------------------------------------
# Select best base pairs (mirrors ScanFold2.0.py)
# ---------------------------------------------------------------------------

best_bps = {}
best_sum_bps = {}
best_sum_bps_means = {}
best_total_window_mean_bps = {}

log_total = open(get_path(outname + ".ScanFold.log"), "w")
log_win = open(get_path(outname + ".ScanFold.FinalPartners.txt"), "w")
log_win.write(
    "i\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
    "\t*Indicates most favorable bp has competition; bp(j) has more favorable partner or is "
    "more likely to be unpaired\n"
)

for k, v in sorted(bp_dict.items()):
    zscore_dict = {}
    pair_dict = {}
    mfe_dict = {}
    ed_dict = {}

    for pair in v:
        partner_key = str(pair.jnucleotide) + "-" + str(pair.jcoordinate)
        x = NucPair(
            pair.inucleotide,
            pair.icoordinate,
            pair.jnucleotide,
            pair.jcoordinate,
            pair.zscore,
            pair.mfe,
            pair.ed,
        )
        try:
            zscore_dict[partner_key].append(pair.zscore)
            mfe_dict[partner_key].append(pair.mfe)
            ed_dict[partner_key].append(pair.ed)
            pair_dict[partner_key].append(x)
        except KeyError:
            zscore_dict[partner_key] = [pair.zscore]
            mfe_dict[partner_key] = [pair.mfe]
            ed_dict[partner_key] = [pair.ed]
            pair_dict[partner_key] = [x]

    sum_z = {k1: sum(v1) for k1, v1 in zscore_dict.items()}
    mean_z = {k1: statistics.mean(v1) for k1, v1 in zscore_dict.items()}
    mean_mfe = {k1: statistics.mean(v1) for k1, v1 in mfe_dict.items()}
    mean_ed = {k1: statistics.mean(v1) for k1, v1 in ed_dict.items()}

    total_windows = sum(len(v1) for v1 in zscore_dict.values())
    num_bp = sum(1 for k1 in zscore_dict if int(k) != int(k1[2:]))

    k_nuc = str(nuc_dict[k].nucleotide)
    log_total.write(
        f"\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED\tSumZ\tSumZ/#TotalWindows\tBPs= {num_bp}\n"
    )
    log_total.write(f"nt-{k}\t-\t{k_nuc}\t{total_windows}\t-\t-\t-\t-\t-\n")

    total_window_mean_z = {}
    for k1, v1 in zscore_dict.items():
        key_i = int(k1[2:])
        total_window_mean_z[k1] = sum(v1) / total_windows
        z_sum = str(round(sum(v1), 2))
        z_avg = str(round(statistics.mean(v1), 2))
        test = str(round(total_window_mean_z[k1], 2))
        k1_mean_mfe = str(round(mean_mfe[k1], 2))
        k1_mean_ed = str(round(mean_ed[k1], 2))
        if int(k) == key_i:
            log_total.write(
                f"{k}\tNoBP\t{k1[:1]}\t{len(v1)}\t{k1_mean_mfe}\t{z_avg}\t{k1_mean_ed}\t{z_sum}\t{test}\n"
            )
        else:
            log_total.write(
                f"{k}\t{key_i}\t{k1[:1]}\t{len(v1)}\t{k1_mean_mfe}\t{z_avg}\t{k1_mean_ed}\t{z_sum}\t{test}\n"
            )
        best_bp_key = min(total_window_mean_z, key=total_window_mean_z.get)

    best_bp_mean_z = mean_z[best_bp_key]
    best_bp_sum_z = sum_z[best_bp_key]
    best_bp_mean_mfe = mean_mfe[best_bp_key]
    best_bp_mean_ed = mean_ed[best_bp_key]
    best_total_window_mean_z = total_window_mean_z[best_bp_key]

    best_bp_data = re.split("-", best_bp_key)
    best_nucleotide = best_bp_data[0]
    best_coordinate = best_bp_data[1]

    best_total_window_mean_bps[k] = NucPair(
        nuc_dict[k].nucleotide,
        nuc_dict[k].coordinate,
        best_nucleotide,
        best_coordinate,
        best_total_window_mean_z,
        best_bp_mean_mfe,
        best_bp_mean_ed,
    )
    best_bps[k] = NucPair(
        nuc_dict[k].nucleotide,
        nuc_dict[k].coordinate,
        best_nucleotide,
        best_coordinate,
        best_bp_mean_z,
        best_bp_mean_mfe,
        best_bp_mean_ed,
    )

log(f"best_bps computed. Detecting competing pairs...")

# ---------------------------------------------------------------------------
# Detect competing partners (mirrors ScanFold2.0.py)
# ---------------------------------------------------------------------------

final_partners = {}
length = len(seq)

for k, v in sorted(best_bps.items()):
    test_k = int(k)
    if sum(test_k == int(vv.jcoordinate) for vv in best_bps.values()) >= 0:
        if len(best_bps) < length * 4:
            subdict = best_total_window_mean_bps
        elif (v.icoordinate - length * 2) >= int(start_coordinate) and (
            v.icoordinate + length * 2
        ) <= int(end_coordinate):
            keys = range(
                int(v.icoordinate - length * 2), int(v.icoordinate + length * 2)
            )
            subdict = {
                k2: best_total_window_mean_bps[k2]
                for k2 in keys
                if k2 in best_total_window_mean_bps
            }
        elif int(v.icoordinate + length * 2) <= int(end_coordinate):
            keys = range(int(start_coordinate), int(v.icoordinate + length * 2) + 1)
            subdict = {
                k2: best_total_window_mean_bps[k2]
                for k2 in keys
                if k2 in best_total_window_mean_bps
            }
        elif v.icoordinate + length * 2 >= int(end_coordinate):
            keys = range(max(1, int(v.icoordinate - length * 2)), int(end_coordinate))
            subdict = {
                k2: best_total_window_mean_bps[k2]
                for k2 in keys
                if k2 in best_total_window_mean_bps
            }
        elif len(best_bps) < length:
            subdict = best_total_window_mean_bps
        else:
            print("Sub-dictionary error")
            raise ValueError("Sub-dictionary error")

        if len(subdict) >= 0:
            comp_pairs_i = competing_pairs(subdict, v.icoordinate)
            comp_pairs_j = competing_pairs(subdict, v.jcoordinate)
            total_pairs = []
            for key, pair in comp_pairs_i.items():
                total_pairs.append(competing_pairs(subdict, pair.jcoordinate))
                total_pairs.append(competing_pairs(subdict, pair.icoordinate))
            for key, pair in comp_pairs_j.items():
                total_pairs.append(competing_pairs(subdict, pair.jcoordinate))
                total_pairs.append(competing_pairs(subdict, pair.icoordinate))

            merged_dict = {}
            i = 0
            for d in total_pairs:
                for k1, v1 in d.items():
                    merged_dict[i] = v1
                    i += 1

            if len(merged_dict) > 0:
                bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
            else:
                bp = NucPair(
                    v.inucleotide,
                    v.icoordinate,
                    v.inucleotide,
                    v.icoordinate,
                    v.zscore,
                    v.mfe,
                    v.ed,
                )

            if int(k) != bp.icoordinate and int(k) != int(bp.jcoordinate):
                log_win.write(
                    f"nt-{k}*:\t{v.icoordinate}\t{v.jcoordinate}\t"
                    f"{round(v.mfe, 2)}\t{round(v.zscore, 2)}\t{round(v.ed, 2)}\n"
                )
                final_partners[k] = NucPair(
                    v.inucleotide,
                    v.icoordinate,
                    v.inucleotide,
                    v.icoordinate,
                    best_bps[bp.icoordinate].zscore,
                    bp.mfe,
                    bp.ed,
                )
            else:
                log_win.write(
                    f"nt-{k}:\t{bp.icoordinate}\t{bp.jcoordinate}\t"
                    f"{round(best_bps[k].mfe, 2)}\t{round(best_bps[k].zscore, 2)}\t{round(best_bps[k].ed, 2)}\n"
                )
                final_partners[k] = NucPair(
                    bp.inucleotide,
                    bp.icoordinate,
                    bp.jnucleotide,
                    bp.jcoordinate,
                    best_bps[bp.icoordinate].zscore,
                    best_bps[bp.icoordinate].mfe,
                    best_bps[bp.icoordinate].ed,
                )
        else:
            continue
    else:
        final_partners[k] = NucPair(
            v.inucleotide,
            v.icoordinate,
            v.jnucleotide,
            v.jcoordinate,
            best_bps[k].zscore,
            best_bps[k].mfe,
            best_bps[k].ed,
        )

log_total.close()
log_win.close()
log(f"Competing pairs resolved. Writing output files...")

# ---------------------------------------------------------------------------
# Write output files (mirrors ScanFold2.0.py)
# ---------------------------------------------------------------------------

strand = 1
name = seq_id

write_bp(
    best_bps, get_path(structure_basename + ".ALL.bp"), start_coordinate, name, minz
)
write_ct(
    final_partners,
    get_path(no_filter_base + ".ct"),
    float(10),
    strand,
    name,
    start_coordinate,
)
write_ct(
    final_partners,
    get_path(minus_1_base + ".ct"),
    float(-1),
    strand,
    name,
    start_coordinate,
)
write_ct(
    final_partners,
    get_path(minus_2_base + ".ct"),
    float(-2),
    strand,
    name,
    start_coordinate,
)
makedbn(get_path(no_filter_base), "NoFilter")
makedbn(get_path(minus_1_base), "Zavg_-1")
makedbn(get_path(minus_2_base), "Zavg_-2")
write_bp(final_partners, get_path(outname + ".bp"), start_coordinate, name, minz)
write_wig_dict(
    final_partners, get_path(outname + ".Zavg.wig"), name, step_size, "zscore"
)
write_wig(mfe_list, step_size, seq_id, get_path(outname + ".scan-MFE.wig"))
write_wig(zscore_list, step_size, seq_id, get_path(outname + ".scan-zscores.wig"))
write_wig(ed_list, step_size, seq_id, get_path(outname + ".scan-ED.wig"))

log(f"All files written. Done in {elapsed()}s")

# Print summary of what was written
for f in [
    f"{no_filter_base}.ct",
    f"{no_filter_base}.dbn",
    f"{minus_1_base}.ct",
    f"{minus_1_base}.dbn",
    f"{minus_2_base}.ct",
    f"{minus_2_base}.dbn",
    f"{structure_basename}.ALL.bp",
    f"{outname}.bp",
    f"{outname}.Zavg.wig",
    f"{outname}.scan-MFE.wig",
    f"{outname}.scan-zscores.wig",
    f"{outname}.scan-ED.wig",
]:
    full = get_path(f)
    exists = "OK" if os.path.exists(full) else "MISSING"
    print(f"  [{exists}] {f}")
