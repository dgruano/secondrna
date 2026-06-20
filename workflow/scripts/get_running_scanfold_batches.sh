#!/usr/bin/env bash
# Print running scanfold2_run_batch jobs from SLURM.
#
# Usage:
#   ./get_running_scanfold_batches.sh [-y [FILE]] [SAMPLE]
#
#   -y [FILE]   Write YAML config instead of TSV. FILE defaults to running_batches.yaml.
#   SAMPLE      Filter to a single sample (outputs bare batch IDs, or YAML for that sample).
#
# Examples:
#   # All samples, TSV to stdout
#   ./get_running_scanfold_batches.sh
#
#   # All samples, write running_batches.yaml
#   ./get_running_scanfold_batches.sh -y
#
#   # All samples, write custom path
#   ./get_running_scanfold_batches.sh -y /tmp/running.yaml
#
#   # Single sample, bare batch IDs to stdout
#   ./get_running_scanfold_batches.sh gencode.v47.repeat.simple
#
#   # Single sample, write YAML
#   ./get_running_scanfold_batches.sh -y gencode.v47.repeat.simple
#
#   snakemake scanfold2_retry_all_mikael --profile profiles/default --configfile running_batches.yaml

PREFIX="rule_scanfold2_run_batch_wildcards_"
YAML=0
YAML_FILE="running_batches.yaml"

# Parse -y [FILE]
if [[ "${1:-}" == "-y" ]]; then
    YAML=1
    shift
    # Next arg is the output file only if it doesn't look like a sample name
    if [[ -n "${1:-}" && "${1}" == *.yaml || "${1}" == *.yml || "${1}" == /* || "${1}" == ./* ]]; then
        YAML_FILE="$1"
        shift
    fi
fi

SAMPLE="${1:-}"

tsv_to_yaml() {
    awk -F'\t' '
        BEGIN { print "scanfold_running_batches:" }
        !seen[$1]++ { print "  " $1 ":" }
        { print "    - \"" $2 "\"" }
    '
}

if [[ -n "$SAMPLE" ]]; then
    RAW=$(squeue -u "$USER" -h -o "%100k" \
        | grep "${PREFIX}${SAMPLE}_" \
        | sed "s/.*${PREFIX}${SAMPLE}_//")

    if [[ $YAML -eq 1 ]]; then
        echo "$RAW" | awk -v sample="$SAMPLE" '
            BEGIN { print "scanfold_running_batches:"; print "  " sample ":" }
            NF    { print "    - \"" $0 "\"" }
        ' > "$YAML_FILE"
        echo "Written: $YAML_FILE"
    else
        echo "$RAW"
    fi
else
    RAW=$(squeue -u "$USER" -h -o "%100k" \
        | grep "${PREFIX}" \
        | sed "s/.*${PREFIX}//" \
        | sed 's/_\([0-9]\{5\}\)/\t\1/')

    if [[ $YAML -eq 1 ]]; then
        echo "$RAW" | tsv_to_yaml > "$YAML_FILE"
        echo "Written: $YAML_FILE"
    else
        echo "$RAW"
    fi
fi
