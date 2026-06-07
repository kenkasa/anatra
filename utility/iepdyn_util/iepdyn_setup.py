#!/usr/bin/env python3

######################################################################
# @package iepdyn_setup
# @brief   IEPDYN input generator 
# @author  Kento Kasahara (KK)
#
# (c) Copyright 2026 Kento Kasahara. All rights reserved.
######################################################################
# usage: ./iepdyn_setup.py
#        (No arguments are required)

def ask_bool(prompt, desc, default=False):
    print(f"{prompt} : {desc}")
    default_str = "y" if default else "n"

    while True:
        val = input(f"(y/n) [{default_str}]: ").strip().lower()

        if val == "":
            return default
        if val in ["y", "yes"]:
            return True
        if val in ["n", "no"]:
            return False

        print("Please enter 'y' or 'n'.")


def ask_float(prompt, desc, default):
    print(f"{prompt} : {desc}")
    val = input(f"[{default}]: ").strip()
    return float(val) if val else default


def ask_int_required(prompt, desc):
    print(f"{prompt} : {desc}")
    while True:
        val = input("(required): ").strip()
        try:
            return int(val)
        except ValueError:
            print("Invalid input. Please enter an integer.")


def ask_int(prompt, desc, default):
    print(f"{prompt} : {desc}")
    val = input(f"[{default}]: ").strip()
    return int(val) if val else default


def ask_str_list(prompt, desc):
    print(f"{prompt} : {desc}")
    val = input("(space separated): ").strip()
    return val.split()


def ask_int_list(prompt, desc):
    print(f"{prompt} : {desc}")
    while True:
        val = input("(space separated integers): ").strip()
        try:
            return [int(x) for x in val.split()]
        except ValueError:
            print("Invalid input. Please enter integers only.")


def format_str_list(lst):
    return " ".join([f'"{x}"' for x in lst])


def format_int_list(lst):
    return " ".join(map(str, lst))


def bool_to_fortran(val):
    return ".true." if val else ".false."


def print_aligned_namelist(name, items):
    """
    items: list of (key, value_string)
    """
    print(f"&{name}")

    if not items:
        print("/")
        return

    maxlen = max(len(k) for k, _ in items)

    for k, v in items:
        print(f" {k.ljust(maxlen)} = {v}")

    print("/")


def main():
    print("=== IEPDYN input generator ===")

    # =========================
    # input_param
    # =========================
    print("\n[input_param]")
    print("Choose input method:")
    print("1: fcv      ! Time-series CV data")
    print("2: flist_cv ! CV-file list")

    while True:
        choice = input("Select 1 or 2: ").strip()
        if choice in ["1", "2"]:
            break
        print("Please enter 1 or 2.")

    if choice == "1":
        while True:
            fcv = ask_str_list("fcv", "Time-series CV data")
            if fcv:
                break
            print("fcv cannot be empty.")
        flist_cv = None
    else:
        while True:
            print("flist_cv : CV-file list")
            flist_cv = input("File name: ").strip()
            if flist_cv:
                break
            print("flist_cv cannot be empty.")
        fcv = None

    # =========================
    # output_param
    # =========================
    print("\n[output_param]")
    print("fhead : header of output file")
    fhead = input("[filehead]: ").strip() or "filehead"

    # =========================
    # option_param
    # =========================
    print("\n[option_param]")

    use_perturbed_traj   = ask_bool        ("use_perturbed_traj",   "whether input files come from perturbed dynamics", False) 
    use_dissociate_state = ask_bool        ("use_dissociate_state", "define dissociate state or not", False)
    use_reflection_state = ask_bool        ("use_reflection_state", "define reflection state or not", False)
    use_product_state    = ask_bool        ("use_product_state",    "define product (absorbing) state or not", False)

    calc_steady          = ask_bool        ("calc_steady",          "calculate steady-state populations or not", False)
    calc_Pint            = ask_bool        ("calc_Pint",            "calculate time integral of Pj analytically or not", False)
    extrapolate          = ask_bool        ("extrapolate",          "calculate the time development of Pj based on the integral equations")
                                                                    
    nstate               = ask_int_required("nstate",               "# of states")
    ndim                 = ask_int         ("ndim",                 "# of dimensions", 1)
    nmol                 = ask_int         ("nmol",                 "# of target molecules", 1)
                                                                    
    dt                   = ask_float       ("dt",                   "Time grid for input CV files", 1.0)
    t_sparse             = ask_float       ("t_sparse",             "Sparse time-grid for computing K-, M-, R-, and P0-functions", 1.0)
    t_range              = ask_float       ("t_range",              "Timescale for K-, M-, R-, and P0-functions", 10.0)

    if extrapolate:
        t_extend         = ask_float       ("t_extend",             "Extended timescale for P- and Q-functions", 100.0)
        dt_tcfout        = ask_float       ("dt_tcfout",            "Time grid for outputting P- and Q-functions", 2.0)

    if use_perturbed_traj:
        f_unperturbed_id = input("[list file that contains unperturbed state ID and Rij-flag for each input file]: ").strip() or ""

    initial_state_ids = ask_int_list("initial_state_ids", "Initial state IDs")

    if use_reflection_state:
        reflection_state_ids = ask_int_list("reflection_state_ids", "Reflection state IDs")
    else:
        reflection_state_ids = None

    if use_dissociate_state:
        dissociate_state_ids = ask_int_list("dissociate_state_ids", "Dissociation state IDs")
    else:
        dissociate_state_ids = None

    if use_product_state:
        product_state_ids = ask_int_list("product_state_ids", "Product (absorbing) state IDs")
    else:
        product_state_ids = None

    # =========================
    # state mode selection
    # =========================
    use_interactive_state = ask_bool(
        "state_input_mode",
        "Do you want to define states interactively now?",
        True
    )

    states = []

    if use_interactive_state:
        print("\n[state]")

        expected_cols = 2 * ndim + 1
        print(f"Each state requires {expected_cols} values.")
        print("Format: xmin xmax (for each CV) + weight")

        while True:
            states = []
            for i in range(nstate):
                while True:
                    line = input(f"State {i+1}: ").strip()
                    tokens = line.split()

                    if len(tokens) != expected_cols:
                        print(f"Error: Expected {expected_cols} values.")
                        continue

                    try:
                        values = [float(x) for x in tokens]
                    except ValueError:
                        print("Error: numeric values required.")
                        continue

                    valid = True
                    for d in range(ndim):
                        if values[2*d] >= values[2*d+1]:
                            print(f"Error: xmin >= xmax in dimension {d+1}.")
                            valid = False
                            break

                    if not valid:
                        continue

                    states.append(line)
                    break

            if len(states) == nstate:
                break
            else:
                print("Mismatch in number of states. Re-enter.")

    # =========================
    # output
    # =========================
    print("\n=== Generated input ===\n")

    # input_param
    items = []
    if fcv is not None:
        items.append(("fcv", format_str_list(fcv)))
    if flist_cv is not None:
        items.append(("flist_cv", f'"{flist_cv}"'))
    print_aligned_namelist("input_param", items)

    # output_param
    items = [
        ("fhead", f'"{fhead}"')
    ]
    print_aligned_namelist("output_param", items)

    # option_param
    items = [
        ("use_perturbed_traj",   bool_to_fortran(use_perturbed_traj)),
        ("use_dissociate_state", bool_to_fortran(use_dissociate_state)),
        ("use_reflection_state", bool_to_fortran(use_reflection_state)),
        ("use_product_state",    bool_to_fortran(use_product_state)),
        ("calc_steady",          bool_to_fortran(calc_steady)),
        ("calc_Pint",            bool_to_fortran(calc_Pint)),
        ("extrapolate",          bool_to_fortran(extrapolate)),
        ("nstate",               str(nstate)),
        ("ndim",                 str(ndim)),
        ("nmol",                 str(nmol)),
        ("dt",                   f"{dt}"),
        ("t_sparse",             f"{t_sparse}"),
        ("t_range",              str(t_range)),
        ("initial_state_ids",    format_int_list(initial_state_ids)),
    ]

    if use_perturbed_traj:
        items.append(("f_unperturbed_id",     f'"{f_unperturbed_id}"'))

    if extrapolate:
        items.append(("t_extend",             str(t_extend)))
        items.append(("dt_tcfout",            f"{dt_tcfout}"))

    if use_reflection_state:
        items.append(("reflection_state_ids", format_int_list(reflection_state_ids)))

    if use_dissociate_state:
        items.append(("dissociate_state_ids", format_int_list(dissociate_state_ids)))

    if use_product_state:
        items.append(("product_state_ids", format_int_list(product_state_ids)))

    print_aligned_namelist("option_param", items)

    # state
    print("&state")
    if use_interactive_state:
        for s in states:
            print(s)
    else:
        print(f" ! Please define {nstate} states manually")
        print(f" ! Each line must contain {2*ndim+1} values")
        print(f" ! Format: xmin xmax (per CV) + weight")
    print("/")


if __name__ == "__main__":
    main()

