#!/bin/bash

# Script: run_get_clusters.sh
# Version: 2.1
# Purpose: Robust execution of get_clusters.sh with advanced debugging

# ---------------------------
# 1. Enhanced Configuration
# ---------------------------
LOG_FILE="cgd_fingerprint.log"
ERROR_FILE="cgd_fingerprint.error"
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
    # Extract src_directory from setup.in
    if [ ! -f "setup.in" ]; then
        echo "ERROR: setup.in not found in $(pwd)" | tee -a "$LOG_FILE"
        return 1
    fi

    src_directory=$(grep -oP '^CGD_SRCDIR\s*=\s*\K.+' setup.in | sed 's/"//g' | sed "s/'//g")
    if [ -z "$src_directory" ]; then
        echo "ERROR: CGD_SRCDIR not properly defined in setup.in" | tee -a "$LOG_FILE"
        return 1
    fi

    # Normalize path (remove trailing slash then add one)
    src_directory="${src_directory%/}/"
    echo "Source directory: $src_directory" | tee -a "$LOG_FILE"

    # Verify directory exists
    if [ ! -d "$src_directory" ]; then
        echo "ERROR: Source directory does not exist: $src_directory" | tee -a "$LOG_FILE"
        echo "DEBUG: Directory contents at $(dirname "$src_directory"):" | tee -a "$LOG_FILE"
        ls -la "$(dirname "$src_directory")" | tee -a "$LOG_FILE"
        return 1
    fi

    return 0
}

if ! validate_src_directory; then
    exit 1
fi

# ---------------------------
# 3. Script Location Debugging
# ---------------------------
locate_script() {
    script_path="${src_directory}get_clusters.sh"
    echo "Looking for script at: $script_path" | tee -a "$LOG_FILE"

    if [ ! -f "$script_path" ]; then
        echo "ERROR: get_clusters.sh not found at expected location" | tee -a "$LOG_FILE"
        echo "DEBUG: Contents of source directory:" | tee -a "$LOG_FILE"
        ls -la "$src_directory" | tee -a "$LOG_FILE"
        
        # Alternative search in common locations
        echo "DEBUG: Searching for get_clusters.sh in common locations..." | tee -a "$LOG_FILE"
        find "$src_directory" -name "get_clusters.sh" -print | tee -a "$LOG_FILE"
        
        return 1
    fi

    echo "Found script at: $script_path" | tee -a "$LOG_FILE"
    return 0
}

if ! locate_script; then
    exit 1
fi

# ---------------------------
# 4. Execution with Debugging
# ---------------------------
echo "---------------------------------------" | tee -a "$LOG_FILE"
echo "Starting get_clusters.sh execution..." | tee -a "$LOG_FILE"
full_cmd="sh \"${src_directory}get_clusters.sh\""
echo "Full command: $full_cmd" | tee -a "$LOG_FILE"

start_time=$(date +%s)
set -x  # Enable command tracing
eval "$full_cmd" >> "$LOG_FILE" 2>> "$ERROR_FILE"
exit_code=$?
set +x  # Disable command tracing
end_time=$(date +%s)
duration=$((end_time - start_time))

# ---------------------------
# 5. Enhanced Post-Execution
# ---------------------------
{
echo "---------------------------------------"
echo "=== Execution Completed ==="
echo "Exit code: $exit_code"
[ $exit_code -eq 0 ] && status="SUCCESS" || status="FAILED"
echo "Status: $status"
echo "Duration: $duration seconds"
echo "End time: $(date +"%Y-%m-%d %T")"
echo "Log file: $LOG_FILE"
echo "Error file: $ERROR_FILE"
echo "---------------------------------------"
} | tee -a "$LOG_FILE"

# Exit with proper code
exit $exit_code