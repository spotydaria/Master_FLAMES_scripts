#!/bin/bash
#SBATCH --job-name=pops_split       # Job name
#SBATCH --time=05:00:00              # Time limit (3 hours)
#SBATCH --mem=32G

input_file="/home/dkramarenk/projects/tools/POPs/POPs_f_Joel/raw/PoPS.features.txt" 
delimiter="\t"      
output_file_prefix="/home/dkramarenk/projects/tools/POPs/POPs_f_Joel/to_munge0/PoPS.features_"

# Count the number of columns in the input file
num_columns=$(head -n 1 "$input_file" | awk -F "$delimiter" '{print NF}')

# Calculate the number of columns for each output file
columns_per_file=$(( (num_columns - 1) / 3 ))

# Split the input file into three output files
for i in $(seq 0 2); do
    start_column=$(( 2 + i * columns_per_file ))
    end_column=$(( 1 + (i + 1) * columns_per_file ))

    # Handle the last output file if the columns don't divide evenly by 3
    if [ $i -eq 2 ] && [ $(( (num_columns - 1) % 3 )) -ne 0 ]; then
        end_column=$num_columns
    fi

    awk -F "$delimiter" -v start_col="$start_column" -v end_col="$end_column" -v OFS="$delimiter" '{
        # Print the first column and the range of columns for the current output file
        printf $1
        for (i = start_col; i <= end_col; ++i) {
            printf OFS $i
        }
        printf "\n"
    }' "$input_file" > "${output_file_prefix}$((i + 1)).txt"
done
