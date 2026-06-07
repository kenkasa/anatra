#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Split flist and funp into directory-based files."
)
parser.add_argument(
    "--flist",
    required=True,
    help="Input flist file"
)
parser.add_argument(
    "--funp",
    required=True,
    help="Input funp file"
)

args = parser.parse_args()

# Input file names
flist_file = args.flist
funp_file = args.funp

# Read flist
with open(flist_file, "r") as f:
    flist_lines = [line.rstrip("\n") for line in f]

# Read funp
with open(funp_file, "r") as f:
    funp_lines = [line.rstrip("\n") for line in f]

# Check the number of lines
if len(flist_lines) != len(funp_lines):
    raise ValueError("The number of lines in flist and funp does not match")

# Counter for each directory
counter = defaultdict(int)

# Process each line
for flist_line, funp_line in zip(flist_lines, funp_lines):

    cols = funp_line.split()

    if len(cols) < 1:
        continue

    # Read the value in the first column
    dir_id = int(float(cols[0]))

    # Zero padding
    dir_name = f"{dir_id:04d}"

    # Create directory
    os.makedirs(dir_name, exist_ok=True)

    # Increment file index within the directory
    counter[dir_name] += 1
    idx = counter[dir_name]

    # Output file names
    flist_out = os.path.join(dir_name, f"{idx:04d}.flist")
    funp_out = os.path.join(dir_name, f"{idx:04d}.funp")

    # Write the corresponding flist line
    with open(flist_out, "w") as f:
        f.write(flist_line + "\n")

    # Write the corresponding funp line
    with open(funp_out, "w") as f:
        f.write(funp_line + "\n")

print("Processing completed successfully")
