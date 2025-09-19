#!/bin/bash

# ./get_celldimension.sh $1
# $1 = path to position.lammpstrj


# File paths
input_file=$1/position.lammpstrj
output_file=$1/celldimensions.txt

# Read the file, find the lines after "pp pp pp", and extract the numbers
lines=$(grep "pp pp pp" -A3 -m 1 "$input_file")

# Initialize result variable
result=""
skip_first_line=true

# Loop through each line and perform subtraction
while read -r line; do
    if $skip_first_line; then
        skip_first_line=false
        continue
    fi
    # Extract two numbers from the line (assuming the numbers are separated by spaces)
    num1=$(echo "$line" | awk '{print $1}')
    num2=$(echo "$line" | awk '{print $2}')
    num1_dec=$(echo "$num1" | awk '{ printf "%.16f", $1 }')
    num2_dec=$(echo "$num2" | awk '{ printf "%.16f", $1 }')

    subtraction=$(echo "$num2_dec - $num1_dec" | bc -l)

    # Append result to the variable
    result+=" $subtraction"
done <<< "$lines"

# Print the result to the output file
echo "$result" > "$output_file"