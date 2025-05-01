#!/bin/bash

# Script: run_get_clusters.sh
# Version: 2.3
# Purpose: Execute get_clusters.sh to full completion before get_histograms.sh

# ---------------------------
# 1. Enhanced Configuration
# ---------------------------
LOG_FILE="cgd_fingerprint.log"
TIMESTAMP=$(date +"%Y-%m-%d %T")

# Load modules
module load python

# Initialize logging
{
echo "=== Starting execution at $TIMESTAMP ==="
echo "Current directory: $(pwd)"
echo "User: $(whoami)"
echo "Host: $(hostname)"
echo "---------------------------------------"
} | tee -a "$LOG_FILE"

# ---------------------------
# 2. Source Directory Validation
# ---------------------------
validate_src_directory() {
    if [ ! -f "setup.in" ]; then
        echo "ERROR: setup.in not found in $(pwd)" | tee -a "$LOG_FILE"
        return 1
    fi

    src_directory=$(grep -oP '^CGD_SRCDIR\s*=\s*\K.+' setup.in | sed 's/"//g' | sed "s/'//g")
    if [ -z "$src_directory" ]; then
        echo "ERROR: CGD_SRCDIR not properly defined in setup.in" | tee -a "$LOG_FILE"
        return 1
    fi

    src_directory="${src_directory%/}/"
    echo "Source directory: $src_directory" | tee -a "$LOG_FILE"

    if [ ! -d "$src_directory" ]; then
        echo "ERROR: Source directory does not exist: $src_directory" | tee -a "$LOG_FILE"
        ls -la "$(dirname "$src_directory")" | tee -a "$LOG_FILE"
        return 1
    fi

    return 0
}

if ! validate_src_directory; then
    exit 1
fi

# ---------------------------
# 3. Script Validation
# ---------------------------
validate_script() {
    local script_path="${src_directory}$1"
    echo "Validating script: $script_path" | tee -a "$LOG_FILE"

    if [ ! -f "$script_path" ]; then
        echo "ERROR: $1 not found" | tee -a "$LOG_FILE"
        ls -la "$src_directory" | tee -a "$LOG_FILE"
        return 1
    fi

    return 0
}

# Validate both scripts exist before proceeding
if ! validate_script "get_clusters.sh"; then
    exit 1
fi

if ! validate_script "get_histograms.sh"; then
    exit 1
fi

# ---------------------------
# 4. Execution with Completion Waiting
# ---------------------------
execute_with_completion() {
    local script_name=$1
    local script_path="${src_directory}${script_name}"
    local pid_file="/tmp/${script_name}.pid"

    echo "---------------------------------------" | tee -a "$LOG_FILE"
    echo "Starting ${script_name} with completion monitoring..." | tee -a "$LOG_FILE"
    echo "Full command: sh \"$script_path\"" | tee -a "$LOG_FILE"

    # Start the script and track its process tree
    sh "$script_path" >> "$LOG_FILE" 2>&1 &
    local main_pid=$!
    echo "${script_name} main PID: $main_pid" | tee -a "$LOG_FILE"

    # Get all child PIDs
    pstree -p $main_pid | grep -oP '\(\K\d+' > "$pid_file"
    echo "Tracking PIDs: $(tr '\n' ' ' < "$pid_file")" | tee -a "$LOG_FILE"

    # Wait for main process and all children
    while kill -0 $main_pid 2>/dev/null || [ -s "$pid_file" ]; do
        # Check if any child processes are still running
        local children_running=0
        while read pid; do
            if kill -0 $pid 2>/dev/null; then
                children_running=1
                break
            fi
        done < "$pid_file"

        [ $children_running -eq 0 ] && break
        sleep 1
    done

    # Clean up
    rm -f "$pid_file"
    wait $main_pid
    local exit_code=$?

    echo "---------------------------------------" | tee -a "$LOG_FILE"
    echo "${script_name} fully completed with exit code: $exit_code" | tee -a "$LOG_FILE"
    return $exit_code
}

# ---------------------------
# 5. Main Execution Flow
# ---------------------------
if execute_with_completion "get_clusters.sh"; then
    echo "Starting get_histograms.sh after confirmed completion of get_clusters.sh" | tee -a "$LOG_FILE"
    execute_with_completion "get_histograms.sh"
    echo "Submitted job to HPC!"
    exit $?
else
    echo "Skipping get_histograms.sh due to get_clusters.sh failure" | tee -a "$LOG_FILE"
    exit 1
fi
