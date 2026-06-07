#!/usr/bin/env python3

######################################################################
# @package cvlist_generator
# @brief   CV file list generator for IEPDYN 
# @author  Kento Kasahara (KK)
#
# (c) Copyright 2026 Kento Kasahara. All rights reserved.
######################################################################
# usage: ./cvlist_generator.py
#        (No arguments are required)
#
# You will be asked
# Enter the root directory to search
# Enter the target extension (e.g., txt)
# Enter the maximum search directory depth
# Enter directory level counted from the deepest directory for grouping (default: 1)
# Did these files come from perturbed dynamics? (y/n)
# [if yes] Enter state numbers used for Rij calculation (space-separated, empty input allowed)
#       
# Enter fraction of files to be listed from each state (default: 1.0 (corresponding to all the files))
# Enter filename to save the file list
# [if yes] Enter filename to save state_id and Rij_flag



from pathlib import Path
from collections import defaultdict
import math


def collect_files(root_dir, extension, max_depth):
    """
    Search files with the specified extension
    up to the given directory depth.
    """

    root = Path(root_dir).resolve()
    matched_files = []

    for path in root.rglob(f"*{extension}"):

        # Calculate depth from the root directory
        rel_parts = path.relative_to(root).parts
        depth = len(rel_parts) - 1  # exclude the file itself

        if depth > max_depth:
            continue

        matched_files.append(path)

    return sorted(matched_files)


def group_files(files, group_level):
    """
    Group files by directory level counted
    from the deepest directory.
    """

    grouped_files = defaultdict(list)

    for path in files:

        parent_parts = path.parent.parts

        # Determine grouping directory
        if len(parent_parts) >= group_level:
            group_dir = parent_parts[-group_level]
        else:
            group_dir = path.parent.name

        grouped_files[group_dir].append(path)

    return grouped_files


def select_fraction(file_list, fraction):
    """
    Select the last fraction of files.

    fraction = 1.0 -> all files
    fraction = 0.9 -> last 90%
    fraction = 0.1 -> last 10%
    """

    n_files = len(file_list)

    if fraction >= 1.0:
        return file_list

    n_select = max(1, math.ceil(n_files * fraction))

    return file_list[-n_select:]


def main():

    print("=== File Search Utility ===\n")

    # (1) Root directory
    root_dir = input("(1) Enter the root directory to search: ").strip()

    # (2) Extension
    extension = input(
        "(2) Enter the target extension (e.g., txt): "
    ).strip()

    # (3) Maximum depth
    max_depth = int(input(
        "(3) Enter the maximum search directory depth: "
    ).strip())

    print("\nSearching files...\n")

    matched_files = collect_files(
        root_dir,
        extension,
        max_depth
    )

    if not matched_files:
        print("No matching files were found.")
        return

    print("=== Matched Files ===\n")

    for path in matched_files:
        print(path)

    print("\n----------------------------------------\n")

    # (4) Grouping level
    group_level = int(input(
        "(4) Enter directory level counted from the deepest directory for grouping "
        "(default: 1): "
    ).strip() or "1")

    grouped_files = group_files(
        matched_files,
        group_level
    )

    print("\n=== Grouped Files ===")

    sorted_groups = sorted(grouped_files.keys())

    # Assign state numbers
    state_map = {}

    for idx, group in enumerate(sorted_groups, start=1):

        state_map[group] = idx

        print(f"\n--- State {idx}: {group} ---")

        for filepath in sorted(grouped_files[group]):
            print(filepath)

    print("\n----------------------------------------\n")

    # Ask whether files are perturbed
    perturbed = input(
        "Did these files come from perturbed dynamics? (y/n): "
    ).strip().lower()

    rij_state_set = set()

    if perturbed == "y":

        print("\nAvailable states:")

        for group in sorted_groups:
            print(f"  {state_map[group]} : {group}")

        user_input = input(
            "\nEnter state numbers used for Rij calculation "
            "(space-separated, empty input allowed): "
        ).strip()

        if user_input:
            rij_state_set = {
                int(x) for x in user_input.split()
            }

    print("\n----------------------------------------\n")

    # Fraction selection
    fraction = float(input(
        "Enter fraction of files to be listed from each state "
        "(default: 1.0 (corresponding to all the files)): "
    ).strip() or "1.0")

    # Prepare output lines
    file_lines = []
    state_lines = []

    for group in sorted_groups:

        state_id = state_map[group]

        rij_flag = 1 if state_id in rij_state_set else 0

        selected_files = select_fraction(
            sorted(grouped_files[group]),
            fraction
        )

        for filepath in selected_files:

            file_lines.append(str(filepath))

            state_lines.append(
                f"{state_id} {rij_flag}"
            )

    print("\n=== Output Preview ===\n")

    for filepath, stateinfo in zip(file_lines, state_lines):
        print(filepath)
        print(stateinfo)

    print("\n----------------------------------------\n")

    # Ask output filename for file list
    filelist_name = input(
        "Enter filename to save the file list: "
    ).strip()

    with open(filelist_name, "w") as f:
        for line in file_lines:
            f.write(line + "\n")

    print(f"File list saved to: {filelist_name}")

    # Ask output filename for state/Rij data only if perturbed
    if perturbed == "y":

        statefile_name = input(
            "Enter filename to save state_id and Rij_flag: "
        ).strip()

        with open(statefile_name, "w") as f:
            for line in state_lines:
                f.write(line + "\n")

        print(f"State/Rij data saved to: {statefile_name}")


if __name__ == "__main__":
    main()
