#!/usr/bin/env python3
"""
Script Name: calc_binding_energy_no_min.py

Description:
  This script calculates a crude, single-point potential energy for two protein chains
  (A and B) in a PDB file. It uses the OpenMM Python API to:
    1. Parse a PDB file with two chains.
    2. Extract each chain into a separate PDB file (chain A and chain B).
    3. Create an OpenMM system for the complex, chain A, and chain B (no minimization, no dynamics).
    4. Compute the static potential energies for each and estimate binding energy.

Usage:
  ./calc_binding_energy_no_min.py input_complex.pdb

Example:
  ./calc_binding_energy_no_min.py my_complex.pdb
"""

import sys
import os
from openmm import app
import openmm
from openmm import unit

def print_header():
    print(__doc__)

def single_point_energy(pdb_file):
    """
    Create an OpenMM system from a PDB file and retrieve a single-point
    (static) potential energy without any minimization or MD steps.
    """
    # Load coordinates/topology
    pdb = app.PDBFile(pdb_file)
    
    # Force field choice
    # (Using amber14 as an example — may be changed as needed)
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Create system with NoCutoff just for demonstration
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    
    # Create an integrator (required, but we won't run dynamics)
    integrator = openmm.LangevinIntegrator(
        300*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds
    )
    
    # Set up simulation
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    
    # Get and return single-point potential energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return energy

def extract_chain(pdb_file, chain_id, out_file):
    """Extract a single chain from a PDB and write it to out_file."""
    with open(pdb_file, 'r') as f_in, open(out_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # PDB chain ID is in column 22 (index 21 in 0-based)
                if line[21].strip() == chain_id:
                    f_out.write(line)
            elif line.startswith('TER') and line[21].strip() == chain_id:
                f_out.write(line)
        f_out.write('END\n')

def main():
    # Check arguments
    if len(sys.argv) < 2:
        print_header()
        sys.exit(1)
    
    complex_pdb = sys.argv[1]
    if not os.path.isfile(complex_pdb):
        print(f"Error: File '{complex_pdb}' not found.")
        sys.exit(1)

    # For demonstration, assume chain IDs 'A' and 'B'.
    chainA_pdb = "chainA.pdb"
    chainB_pdb = "chainB.pdb"

    # Extract chain A and B
    print(f"Running: Extracting chain A into {chainA_pdb}")
    extract_chain(complex_pdb, 'A', chainA_pdb)
    
    print(f"Running: Extracting chain B into {chainB_pdb}")
    extract_chain(complex_pdb, 'B', chainB_pdb)
    
    # Calculate single-point energies
    print(f"Running: Calculating single-point energy for the complex {complex_pdb} ...")
    E_complex = single_point_energy(complex_pdb)
    
    print(f"Running: Calculating single-point energy for chain A {chainA_pdb} ...")
    E_A = single_point_energy(chainA_pdb)
    
    print(f"Running: Calculating single-point energy for chain B {chainB_pdb} ...")
    E_B = single_point_energy(chainB_pdb)
    
    # Compute approximate binding energy
    dE_bind = E_complex - (E_A + E_B)
    
    print("---- Results ----")
    print(f"E_complex (kJ/mol): {E_complex:.2f}")
    print(f"E_chainA  (kJ/mol): {E_A:.2f}")
    print(f"E_chainB  (kJ/mol): {E_B:.2f}")
    print(f"Approx. Binding Energy (kJ/mol): {dE_bind:.2f}")

if __name__ == "__main__":
    main()
