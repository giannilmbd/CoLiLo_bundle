#!/bin/bash

# File to process
FILE="my_toolbox.f90"

# Backup the original file
cp "$FILE" "${FILE}.bak"

# Process the file
awk '
    # Match function or subroutine declarations
    /^\s*(function|subroutine)\s+/ {
        # Save the current line
        line = $0
        # Read the next line to check for "use my_kinds_mod"
        getline next_line
        if (next_line !~ /^\s*use\s+my_kinds_mod/) {
            # If "use my_kinds_mod" is not present, add it
            print line
            print "    use my_kinds_mod"
            print next_line
        } else {
            # Otherwise, keep the original lines
            print line
            print next_line
        }
        next
    }
    # Print all other lines as-is
    { print }
' "$FILE" > "${FILE}.tmp"

# Replace the original file with the modified file
mv "${FILE}.tmp" "$FILE"

echo "Processed $FILE. A backup was saved as ${FILE}.bak."