#!/usr/bin/env python
"""
select_residues.py

Define a PyMOL command  select_residues
that builds a selection called  my_sel  from an arbitrary list of
residue IDs.  Call it like:

    select_residues 946 1048 1063 ... 1280
    select_residues [946 1048 1063 ... 1280]
    select_residues (946 1048 1063 ... 1280)

If called with no residue list, prints this header.
"""

from pymol import cmd
import sys

def select_residues(*args, **kw):
    """
    Implementation of the PyMOL command.

    *args : (_self, "946", "1048", ...)
    **kw  : may contain another '_self' on some old PyMOL builds
    """
    # --------------------------------------------
    # 1.  Pull out only the string tokens
    # --------------------------------------------
    tokens = [a for a in args if isinstance(a, str)]
    if not tokens:
        print(__doc__)
        return

    # --------------------------------------------
    # 2.  If user typed a single *bracketed* blob,
    #     re-split it on whitespace.
    #     Otherwise keep the individual tokens.
    # --------------------------------------------
    if len(tokens) == 1:
        blob = tokens[0].strip()
        if (blob.startswith('[') and blob.endswith(']')) or \
           (blob.startswith('(') and blob.endswith(')')):
            blob = blob[1:-1].strip()
        tokens = blob.split()

    # --------------------------------------------
    # 3.  Strip leftover bracket chars on first / last token
    # --------------------------------------------
    tokens[0]  = tokens[0].lstrip('[(')
    tokens[-1] = tokens[-1].rstrip('])')

    # --------------------------------------------
    # 4.  Convert to integers → build “+” string
    # --------------------------------------------
    try:
        resi_list = [int(t) for t in tokens]
    except ValueError:
        print("\nError: every item must be an integer residue ID.\n")
        print(__doc__)
        return

    resi_str = "+".join(str(r) for r in resi_list)

    # --------------------------------------------
    # 5.  Make the selection
    # --------------------------------------------
    sel_name = "my_sel"
    cmd.select(sel_name, f"resi {resi_str}")
    print(f"Running: select {sel_name}, resi {resi_str}")

# register with PyMOL
cmd.extend("select_residues", select_residues)

# ------------------------------------------------------------
# Optional: run from a normal terminal for a “dry run”
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(__doc__)
    else:
        # mimic PyMOL: pass everything after the script name as strings
        dummy_self = object()
        select_residues(dummy_self, *sys.argv[1:])
