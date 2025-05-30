# Glycine-Nucleation-MD-Simulation

Master Thesis Project for Computational Science at the University of Amsterdam.  
We study the nucleation event of the glycine amino acid through Molecular Dynamics Simulation.  
We find that Glycine follows a two-step nucleation process:  
During the first step the molecules form an amorphous liquid-like cluster.  
In the second step this cluster crystallizes from within until the entire system is in a solid crystal state.  

<table>
  <tr>
    <td align="center">
      <img src="md_images/droplet.png" alt="Amorphous Cluster" width="500"><br>
      <span>Step 1: Amorphous cluster formation</span>
    </td>
    <td align="center">
      <img src="md_images/crystal.png" alt="Crystallization" width="500"><br>
      <span>Step 2: Crystallization within cluster</span>
    </td>
  </tr>
</table>

This repository contains code for creating coordinate files for varying molecular configurations of glycine in water.  
These coordinate files can be used as a starting configuration for Molecular Dynamics Simulation.  
The files are constructed to work with the [GROMACS](https://www.gromacs.org) MD engine.

---

## Directory Structure
```md
Glycine-Nucleation_MD-Simulation
├── Analysis
│   ├── Figures                                  
│   ├── combine_trajectories.ipynb
│   ├── crystal_analyis.ipynb
│   └── plots.ipynb
├── Molecular_System
│   ├── Data
│   ├── glycine_crystal_topology.py
│   ├── glycine_monomer_topology.py
│   ├── glycine_sphere_topology.py
│   └── position_restraints.py
└── md_images


-- Figures/ - Folder where figures and plots are stored
-- combine_trajectories.ipynb - Concatenates two molecular dynamics trajectories to produce a longer trajectory
-- crystal_analysis.ipynb - Notebook for computing and visualization the degree of crystallinity for a molecular system
-- plots.ipynb - Notebook for visualizing measures of the molecular systems
-- Data/ - Contains atomic positions of the glycine molecule and crystals, obtained from experimental studies
-- glycine_crystal_topology.py - Produces coordinates for a perfect glycine crystal in water
-- glycine_monomer_topology.py - Produces coordinates for dispersed glycine monomers in water
-- glycine_sphere_topology.py - Produces coordinates for a spherical glycine crystal with dispersed glycine molecules around it


```


