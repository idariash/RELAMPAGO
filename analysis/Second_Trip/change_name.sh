#!/bin/bash

# Specify the directory where your files are located
directory="./CSU/CHIVO/res/data/test_data_quality"

# Change to the directory
cd "~" || exit
cd "$directory" || exit

# Specify the part of the filename you want to replace
old_text="Link to "

# Specify the new part of the filename
new_text=""

# Loop through all files in the directory
for file in *
do
  # Check if it's a file (not a directory)
  if [ -f "$file" ]; then
    # Perform the filename replacement
    new_name=$(echo "$file" | sed "s/$old_text/$new_text/")
    
    # Rename the file
    mv "$file" "$new_name"
    
    echo "Renamed $file to $new_name"
  fi
done
