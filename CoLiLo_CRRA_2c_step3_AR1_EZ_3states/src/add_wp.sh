#!/bin/bash

# File to process
FILE="my_toolbox.f90"

# Backup the original file
cp "$FILE" "${FILE}.bak"

# Use sed to add _wp to numeric literals
sed -i -E '
    # Match numeric literals and add _wp if not already present
    # Ensure it does not modify variable names or strings
    s/([^a-zA-Z0-9_])([0-9]+(\.[0-9]*)?([eEdD][+-]?[0-9]+)?)([^a-zA-Z0-9_]|$)/\1\2_wp\5/g
' "$FILE"

echo "Processed $FILE. A backup was saved as ${FILE}.bak."