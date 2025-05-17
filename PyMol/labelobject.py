#!/usr/bin/env python
"""
labelobject.py - PyMOL plugin to create a pseudoatom at the center of mass of an object and label it by name within the original object.

Usage:
    Within PyMOL: run labelobject.py
    Then at the PyMOL prompt: labelobject <object_name>
"""
import pymol
from pymol import cmd, stored


def get_residue_colors(sele):
    """
    get_residue_colors(sele) -> dict
    Iterates over selection, collects C-alpha atom colors, and returns a mapping
    from (chain, resi) to RGB color tuple.
    """
    pymol.stored.colors = []
    cmd.iterate(sele, "stored.colors.append((chain, resi, name, color))")
    res_colors = {}
    for chain, resi, name, color in pymol.stored.colors:
        if name == 'CA':  # C-alpha atom
            res_colors[(chain, resi)] = cmd.get_color_tuple(color)
    return res_colors


def labelobject(obj=""):
    """
    labelobject <object_name>
    Creates a pseudoatom at the center of mass of 'object_name', labels it with the object name,
    and colors the label to match the object's carbon atom color, all within the same object.
    """
    if not obj:
        print(__doc__)
        return
    print(f"Running: Creating pseudoatom and label for object '{obj}'")

    # compute center of mass
    try:
        coords = cmd.centerofmass(obj)
    except Exception as e:
        print(f"Error: Could not compute center of mass for '{obj}': {e}")
        return

    # define label atom name
    label_name = f"{obj}_label"

    # create pseudoatom inside original object using named arguments
    cmd.pseudoatom(name=label_name, object=obj, pos=list(coords), label=obj)

    # hide sphere representation of the pseudoatom, show only its label
    sel = f"{obj} and name {label_name}"
    cmd.hide("everything", sel)
    cmd.show("label", sel)

    # color label to match carbon atom color
    try:
        res_colors = get_residue_colors(f"{obj} and elem C")
        # rgb_tuple = next(iter(res_colors.values()))
        rgb_tuple = next(iter(res_colors.values()), None)
        if rgb_tuple:
            custom_color = f"{label_name}_color"
            rgb_list = list(rgb_tuple)
            cmd.set_color(custom_color, rgb_list)
            # set the label_color property for the object
            cmd.set("label_color", custom_color, obj)
        else:
            # fallback if no color found
            cmd.set("label_color", "carbon", obj)
    except Exception:
        # fallback if error
        cmd.set("label_color", "carbon", obj)

# register the command
cmd.extend("labelobject", labelobject)
