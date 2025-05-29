import numpy as np
import os
import math

#Import own conversion functions
from glycine_monomer_topology import convert_topology

#Import own crystal functions
from glycine_crystal_topology import count_molecules

#Import position restrain functions
from position_restraints import read_gro_file, get_posre_indices, atoms_within_range

def write_gro_file(filename, atom_coordinates, atom_names, distance, box_size):
    """
    Writes the coordinates of the atoms to a gro file. Works slightly different from the 
    same named function defined in 'glycine_crystal_topology.py', because of the input datastructure.

    PARAMETERS:
    -----------
        filename: str
            The name of the output gro file.
        atom_coordinates: list
            List of np.arrays containing the xyz coordinates of the atoms.
        atom_names: list
            List containing the atom_names .
        distance: float
            The radius of the sphere in nanometers.
        box_size: float
            Length of the simulation box in nanometers.
    """

    #initialize index count
    index = 1

    #Open the file
    with open(filename, 'w') as file:
        #Calculate number of atoms
        num_elements = len(atom_coordinates)

        #Write title of the title of the file
        file.write(f"Spherical Glycine cluster of radius {distance} nm\n")

        #Write the number of atoms
        file.write('   ' + str(num_elements) + '\n')

        #Loop through the cell_grid and write the coordinates and atom names to file
        for coord, atom in zip(atom_coordinates, atom_names):
                    
            #Use f string to get correct spacing between entries
            line = f"    1GLY{atom:>7}{index:>5}{float(coord[0]):>8.3f}{float(coord[1]):>8.3f}{float(coord[2]):>8.3f}\n"
            file.write(line)
            #Increment index
            index+=1

        #Write box size to file
        file.write(f"   {box_size[0]:<9.5f}{box_size[1]:<9.5f}{box_size[2]:<9.5f}")


def remove_molecules(posre_indices, atom_coordinates, atom_names):
    """Removes atom coordinates based on a list of indices.

    PARAMETERS:
    -----------
        posre_indices: list
            List of atom indices to keep.
        atom_coordinates: list
            List of np.arrays containing the xyz coordinates of the atoms.
        atom_names: list
            List containing the atom_names .
    
    RETURNS:
    --------
        atom_coordinates: list
            Shortened list of np.arrays containing the xyz coordinates of the atoms.
        atom_names: list
            Shortened list containing the atom_names. 
    """

    # Create a mask to keep track of which indices to keep
    mask = [False] * len(atom_coordinates)

    # Set the mask for indices in posre_indices to True
    for idx in posre_indices:
        mask[idx] = True

    # Filter out coordinates and atom names based on the mask
    atom_coordinates = [coord for i, coord in enumerate(atom_coordinates) if mask[i]]
    atom_names = [name for i, name in enumerate(atom_names) if mask[i]]

    return atom_coordinates, atom_names


def calc_concentration(box_size, radius, num_insert):
    """Estimates the concentration of dispersed glycine around a spherical crystal.

    PARAMETERS:
    -----------
        box_size: float
            Length of the simulation box in nanometers.
        radius: float
            The radius of the sphere in nanometers.
        num_inserts: int
            Number of monomer glycine that is inserted around the sphere.
    
    RETURNS:
    --------
        concentration: float
            Concentration of glycine monomers in solution in mol/L. 
    """ 
    #Calculate the volume of the box in nm^3
    box_volume = box_size ** 3

    #Calculate the volume of the glycine sphere in nm^3
    sphere_volume = (4/3) * math.pi * radius ** 3

    #Calculate the volume of water and dispersed glycine in liters ((1 nm^3 = 1e-24 L))
    water_volume = (box_volume - sphere_volume) * 1e-24

    #Calculate the number of moles of dispersed glycine (1 mole = 6.022e23 molecules)
    num_insert_mol = num_insert / 6.022e23

    #Calculate the concentration of glycine in mol/L
    concentration = num_insert_mol / water_volume

    return round(concentration,2)


def build_spherical_cluster(system, distance, point = None):
    """Builds a molecular system of a spherical glycine cluster of certain diameter,
    based on a crystal system and saves it to a .gro file. 

    PARAMETERS:
    -----------
        system: str
            Path to the crystal .gro system to use as basis.
        distance: float
            Radius (in nanometer) within molecules are kept.
        point: np.array
            Coordinate from which to draw the radius.
    """
    #Path to crystal folder
    input_path = "Data/Output/System/Crystal/"

    #Path to sphere folder
    output_path = "Data/Output/System/Sphere/"

    #Read the square crystal system
    atom_coordinates, atom_names, box_size, num_atoms = read_gro_file(input_path+system)
    # print(atom_coordinates)
    
    print(atom_names)

    #By default the point from which to include molecules is the center of the box
    if point == None:
        point = box_size/2

    #Get the atoms that are within the range of the point
    atom_range_bools = atoms_within_range(atom_coordinates, point, distance)

    #Get the indices of the molecules within range
    posre_indices = get_posre_indices(atom_range_bools, atom_names)

    #Remove the molecules that are not within range
    atom_coordinates, atom_names = remove_molecules(posre_indices, atom_coordinates, atom_names)

    #Write the coordinates to file
    write_gro_file(f"{output_path}{system}_{distance}nm_sphere.gro", atom_coordinates, atom_names, distance, box_size)



def build_spherical_system(system, distance, num_insert=np.inf, point=None, add_ions=False,box_size = 5.0, ions=0):
    """Builds a molecular system of a spherical glycine cluster from a perfect crystal
    that is surrounded by free floating glycine molecules and water, with an option to add ions.
    
    PARAMETERS:
    -----------
        system: str
            Path to the crystal .gro system to use as basis.
        distance: float
            Radius (in nanometers) within which molecules are kept.
        num_insert: int, optional
            Number of glycine molecules to insert (default is the difference between original and spherical cluster).
        point: np.array, optional
            Coordinate from which to draw the radius.
        add_ions: bool, optional
            If True, ions (Na+ and Cl-) will be added to neutralize the system and/or adjust concentration.
        box_size: float
            Length of the simulation box in nanometers.
        ions: int, optional
            Ion concentration to add (default is 0.15 mol/L).
    """
    input_path = "Data/Output/System/Crystal/"
    output_path = "Data/Output/System/Sphere/"

    # Build a spherical cluster by cutting a sphere out of the crystal system
    build_spherical_cluster(system, distance, point)

    # Get number of glycine molecules in original and spherical cluster systems
    num_gly_or, _ = count_molecules(f"{input_path}{system}.gro")
    num_gly_sp, _ = count_molecules(f"{output_path}{system}_{distance}nm_sphere.gro")

    # Determine number of glycine molecules to insert
    if num_insert == np.inf:
        num_insert = num_gly_or - num_gly_sp

    # Insert free-floating glycine if necessary
    if num_insert > 0:
        os.system(
            f"gmx insert-molecules "
            f"-f {output_path}{system}_{distance}nm_sphere.gro "
            f"-ci Data/Input/System/glycine_match.pdb "
            f"-nmol {num_insert} "
            f"-o {output_path}{system}_{distance}nm_sphere_insert_{num_insert}.gro"
        )
        system_path = f"{output_path}{system}_{distance}nm_sphere_insert_{num_insert}.gro"
    else:
        system_path = f"{output_path}{system}_{distance}nm_sphere.gro"

    # Solvate the system
    os.system(
        f"gmx solvate "
        f"-cp {system_path} "
        f"-cs spc216.gro "
        f"-o {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv.gro"
    )

    # Solvated system path
    solvated_system = f"{output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv.gro"

    # Create or update topology file
    num_gly, num_sol = count_molecules(solvated_system)
    convert_topology(f"{output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv.top", num_gly, num_sol)

    # Add ions if specified
    if add_ions:
        # Prepare for ion insertion by generating a .tpr file
        os.system(
            f"gmx grompp "
            f"-f Data/Input/Simulation/minim.mdp "  # Use a simple mdp file for energy minimization or ionization
            f"-c {solvated_system} "
            f"-p {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv.top "
            f"-o {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv_ions.tpr"
        )

        # Add Na+ and Cl- ions to neutralize or adjust concentration
        os.system(
            f"echo 13 | gmx genion "
            f"-s {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv_ions.tpr "
            f"-p {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv.top "
            f"-pname NA -nname CL -nn {ions} -np {ions} "
            f"-o {output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv_ions.gro "
            f"-rmin {ions}"
        )

        # Update the solvated system path to include ions
        solvated_system = f"{output_path}{system}_{distance}nm_sphere_insert_{num_insert}_solv_ions.gro"

    # Calculate and report on final system
    gly_conc = calc_concentration(box_size, distance, num_insert)
    print(f"Created a glycine crystal with a {distance}nm radius consisting of {num_gly_sp} glycine molecules.\n"
          f"The crystal is surrounded by a {box_size}x{box_size}x{box_size} nm box with {num_sol} water molecules and {num_insert} dispersed glycine.\n"
          f"The concentration of glycine is {gly_conc} mol/L.")

    if add_ions:
        print(f"Added {ions} ions (Na+ and Cl-)")



#PARAMETERS:
rx, ry, rz = 3, 3, 3
morph_type = "alpha"
box_size = 6.0
distance = 0.8
system = f"{morph_type}_glycine_crystal_{rx}_{ry}_{rz}_box_{box_size}"
num_insert = 312
ions =  0

build_spherical_system(system, distance, num_insert=num_insert, add_ions=False, box_size=box_size, ions=ions)


