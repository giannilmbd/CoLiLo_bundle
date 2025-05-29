#!/bin/bash

# Array of global variables to search for
declare -a globals=(
    "v1" "v2" "v1_test" "v1_new" "v2_new" "pf_sols" "phv"
    "Ev1" "Ev2" "Ev1m1" "Ev2m1" "C1v" "Y1v" "L1v" "Es1" "OmegaSplusv"
    "matrix_csv_v1" "matrix_csv_v2" "matrix_csv_s1" "matrix_csv_pf"
    "matrix_csv_aut_v1" "matrix_csv_aut_v2" "matrix_csv_aut_pf"
    "matrix_coefs_csv_v1" "matrix_coefs_csv_v2" "matrix_coefs_csv_s1"
    "resh_dummy" "coeff_v1" "coeff_v2" "coeff_s1"
    "tx_v1" "ty_v1" "tz_v1" "tq_v1" "bcoef4_v1"
    "tx_v2" "ty_v2" "tz_v2" "tq_v2" "bcoef4_v2"
    "C1k" "C2k" "Y1k" "Y2k" "pfk" "phk" "Qk"
    "C1v_aut" "C2v_aut" "L1v_aut" "L2v_aut" "Y1v_aut" "Y2v_aut" "Qv_aut"
    "pfv_aut" "phv_aut" "v1_aut" "v1_aut_new" "pf_sols_aut"
    "pkh" "pkf" "pkh_new" "pkf_new" "pkh_aut" "pkf_aut" "pkh_aut_new" "pkf_aut_new"
    "v2_aut" "v2_aut_new" "s1" "s1_new" "s_kappa" "A1" "A2" "OmegaS"
    "coeff_v1_aut" "coeff_v2_aut"
)

# Create a temporary file to store results
temp_file=$(mktemp)

# Loop through each global variable
for global in "${globals[@]}"; do
    echo "Checking usage of $global..."
    found=0
    
    # Search in all .f90 and .F90 files, excluding globals files
    for file in $(find . -type f \( -name "*.f90" -o -name "*.F90" \) -not -name "globals*"); do
        if grep -q "\b$global\b" "$file"; then
            found=1
            echo "  Found in: $file"
        fi
    done
    
    # If not found in any file, add to results
    if [ $found -eq 0 ]; then
        echo "$global" >> "$temp_file"
    fi
done

# Print results
echo -e "\nPotentially unused global variables:"
if [ -s "$temp_file" ]; then
    cat "$temp_file"
else
    echo "None found - all variables appear to be used somewhere in the codebase."
fi

# Clean up
rm "$temp_file" 