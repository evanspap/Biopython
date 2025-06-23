#!/usr/bin/env python
"""
labelobject.py - PyMOL plugin to create a pseudoatom at the center of mass of an object and label it by name.

Usage:
    Within PyMOL: run labelobject.py
    Then at the PyMOL prompt: labelobject <object_name>

Example:
    run labelobject.py
    labelobject myProt
"""
from pymol import cmd, stored

# Gather CA atom colors for each residue in a selection
def get_residue_colors(sele):
    stored.colors = []
    cmd.iterate(sele, "stored.colors.append((chain, resi, name, color))")
    res_colors = {}
    for chain, resi, name, color in stored.colors:
        if name == 'CA':  # c-alpha atom
            res_colors[(chain, resi)] = cmd.get_color_tuple(color)
    return res_colors

# Compute center of mass manually from atomic coordinates
def compute_center_of_mass(obj):
    model = cmd.get_model(obj)
    if not model.atom:
        raise RuntimeError(f"No atoms found in '{obj}'")
    x_sum = y_sum = z_sum = 0.0
    for atom in model.atom:
        x_sum += atom.coord[0]
        y_sum += atom.coord[1]
        z_sum += atom.coord[2]
    n = len(model.atom)
    return (x_sum / n, y_sum / n, z_sum / n)

# Main labelobject command
def labelobject(obj=""):
    """
    labelobject <object_name>
    Creates a pseudoatom at the center of mass of 'object_name' and labels it with the object name.
    The label color will match the object's carbon atom color.
    """
    if not obj:
        print(__doc__)
        return

    print(f"Running: Creating pseudoatom and label for object '{obj}'")
    # compute center of mass manually
    try:
        coords = compute_center_of_mass(obj)
    except Exception as e:
        print(f"Error: Could not compute center of mass for '{obj}': {e}")
        return

    label_name = f"{obj}_label"
    # remove existing label atom inside object
    try:
        cmd.remove(f"{obj} and name {label_name}")
    except Exception:
        pass

    # create pseudoatom inside the existing object
    cmd.pseudoatom(name=label_name, object=obj, pos=list(coords), label=obj)
    # hide everything for that atom and show only its label
    cmd.hide("everything", f"{obj} and name {label_name}")
    cmd.show("label", f"{obj} and name {label_name}")

    # match label color to object's carbon atom colors (CA atoms)
    try:
        sele = f"{obj} and elem C and name CA"
        res_colors = get_residue_colors(sele)
        if res_colors:
            rgb = next(iter(res_colors.values()))
            color_name = f"{obj}_label_color"
            cmd.set_color(color_name, list(rgb))
            cmd.set("label_color", color_name, f"{obj} and name {label_name}")
        else:
            cmd.set("label_color", "carbon", f"{obj} and name {label_name}")
    except Exception:
        try:
            cmd.set("label_color", "carbon", f"{obj} and name {label_name}")
        except Exception as e:
            print(f"Warning: Could not set fallback label color: {e}")

# register the command with PyMOL
cmd.extend("labelobject", labelobject)
