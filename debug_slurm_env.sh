#!/bin/bash
# Script to diagnose SLURM environment inheritance on different clusters

# Usage: ./debug_slurm_env.sh <partition> [time_limit]
# Example: ./debug_slurm_env.sh xeon24el8_test 5

if [ $# -lt 1 ]; then
    echo "ERROR: Partition name required"
    echo ""
    echo "Usage: $0 <partition> [time_limit_minutes]"
    echo ""
    echo "Example:"
    echo "  $0 xeon24el8_test"
    echo "  $0 xeon24el8_test 5"
    echo ""
    echo "Available partitions on this system:"
    sinfo -o "%P" | tail -n +2 || echo "  (Could not list partitions)"
    exit 1
fi

PARTITION="$1"
TIME_LIMIT="${2:-5}"  # Default to 5 minutes if not specified

# Output directory - use current directory (shared filesystem)
OUTPUT_DIR="$(pwd)/slurm_diagnostic_output"
mkdir -p "$OUTPUT_DIR"

echo "=== SLURM Environment Inheritance Diagnostic ==="
echo "Partition: $PARTITION"
echo "Time limit: $TIME_LIMIT minutes"
echo "Output directory: $OUTPUT_DIR"
echo ""

echo "1. Check SLURM configuration for environment defaults:"
echo "---------------------------------------------------"
if command -v scontrol &>/dev/null; then
    scontrol show config | grep -i "propagate\|export" || echo "No relevant config found"
else
    echo "scontrol not available"
fi
echo ""

echo "2. Submit a test job to see what environment it receives:"
echo "---------------------------------------------------"

# Create inline diagnostic script (avoids /tmp issues between login and compute nodes)
DIAGNOSTIC_SCRIPT='#!/bin/bash
echo "=== Job Environment Report ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

echo "--- Key Environment Variables ---"
echo "CONDA_PREFIX: ${CONDA_PREFIX:-<not set>}"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "CONDA_SHLVL: ${CONDA_SHLVL:-<not set>}"
echo "PYTHONPATH: ${PYTHONPATH:-<not set>}"
echo "VIRTUAL_ENV: ${VIRTUAL_ENV:-<not set>}"
echo "LOADEDMODULES: ${LOADEDMODULES:-<not set>}"
echo ""

echo "--- All CONDA_* variables ---"
env | grep CONDA || echo "None"
echo ""

echo "--- Module command available? ---"
type module &>/dev/null && echo "YES - module is available" || echo "NO - module NOT available"
echo ""

echo "--- PATH (first 5 entries) ---"
echo "$PATH" | tr ":" "\n" | head -5
'

# Common sbatch options
SBATCH_OPTS="--parsable --partition=$PARTITION --time=$TIME_LIMIT --nodes=1 --ntasks=1 --account=hamilmater"

# Submit with DEFAULT export behavior
SUBMIT_OUTPUT=$(sbatch $SBATCH_OPTS --output="${OUTPUT_DIR}/default_%j.out" --wrap="$DIAGNOSTIC_SCRIPT" 2>&1)
JOB1=$(echo "$SUBMIT_OUTPUT" | grep -oE '^[0-9]+$' | tail -1)
if [ -n "$JOB1" ]; then
    echo "Job submitted with DEFAULT export: $JOB1"
    echo "  Output will be in: ${OUTPUT_DIR}/default_${JOB1}.out"
else
    echo "Job submission FAILED with DEFAULT export"
    echo "  Error: $SUBMIT_OUTPUT"
    JOB1=""
fi

# Submit with --export=ALL
SUBMIT_OUTPUT=$(sbatch $SBATCH_OPTS --export=ALL --output="${OUTPUT_DIR}/exportall_%j.out" --wrap="$DIAGNOSTIC_SCRIPT" 2>&1)
JOB2=$(echo "$SUBMIT_OUTPUT" | grep -oE '^[0-9]+$' | tail -1)
if [ -n "$JOB2" ]; then
    echo "Job submitted with --export=ALL: $JOB2"
    echo "  Output will be in: ${OUTPUT_DIR}/exportall_${JOB2}.out"
else
    echo "Job submission FAILED with --export=ALL"
    echo "  Error: $SUBMIT_OUTPUT"
    JOB2=""
fi

# Submit with --export=NONE
SUBMIT_OUTPUT=$(sbatch $SBATCH_OPTS --export=NONE --output="${OUTPUT_DIR}/exportnone_%j.out" --wrap="$DIAGNOSTIC_SCRIPT" 2>&1)
JOB3=$(echo "$SUBMIT_OUTPUT" | grep -oE '^[0-9]+$' | tail -1)
if [ -n "$JOB3" ]; then
    echo "Job submitted with --export=NONE: $JOB3"
    echo "  Output will be in: ${OUTPUT_DIR}/exportnone_${JOB3}.out"
else
    echo "Job submission FAILED with --export=NONE"
    echo "  Error: $SUBMIT_OUTPUT"
    JOB3=""
fi

echo ""
echo "3. Wait for jobs to complete..."

# Wait for jobs with a timeout
WAIT_TIME=30
echo "Waiting up to $WAIT_TIME seconds for jobs to complete..."

for i in $(seq 1 $WAIT_TIME); do
    sleep 1

    # Check if all submitted jobs are done
    ALL_DONE=true
    for JOB in "$JOB1" "$JOB2" "$JOB3"; do
        if [ -n "$JOB" ]; then
            if squeue -j "$JOB" &>/dev/null; then
                ALL_DONE=false
                break
            fi
        fi
    done

    if $ALL_DONE; then
        echo "All jobs completed after $i seconds"
        break
    fi

    # Show progress every 5 seconds
    if [ $((i % 5)) -eq 0 ]; then
        echo "  Still waiting... ($i/$WAIT_TIME seconds)"
    fi
done

echo ""
echo "4. Results:"
echo "---------------------------------------------------"

if [ -n "$JOB1" ] && [ -f "${OUTPUT_DIR}/default_${JOB1}.out" ]; then
    echo "DEFAULT export behavior:"
    cat "${OUTPUT_DIR}/default_${JOB1}.out"
    echo ""
else
    echo "DEFAULT export: No output file found"
    [ -z "$JOB1" ] && echo "  (Job was not submitted)"
    echo ""
fi

if [ -n "$JOB2" ] && [ -f "${OUTPUT_DIR}/exportall_${JOB2}.out" ]; then
    echo "--export=ALL behavior:"
    cat "${OUTPUT_DIR}/exportall_${JOB2}.out"
    echo ""
else
    echo "--export=ALL: No output file found"
    [ -z "$JOB2" ] && echo "  (Job was not submitted)"
    echo ""
fi

if [ -n "$JOB3" ] && [ -f "${OUTPUT_DIR}/exportnone_${JOB3}.out" ]; then
    echo "--export=NONE behavior:"
    cat "${OUTPUT_DIR}/exportnone_${JOB3}.out"
    echo ""
else
    echo "--export=NONE: No output file found"
    [ -z "$JOB3" ] && echo "  (Job was not submitted)"
    echo ""
fi

echo "=== Diagnostic Complete ==="
echo ""
echo "INTERPRETATION:"
echo "- If DEFAULT and --export=ALL show CONDA_* variables: Cluster exports environment by default"
echo "- If DEFAULT and --export=NONE are similar: Cluster uses clean environment by default"
echo "- If module command is NOT available in any case: Need to initialize modules in preamble"