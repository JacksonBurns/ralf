import os
import math

import scoria
import numpy as np
from scipy.spatial.distance import cdist

# input args
filename = os.path.join(os.getcwd(), "test", "data", "segphos.pdb")
axisatom1 = 5  # this is why Python 0-indexing
axisatom2 = 22
cutoff = 2

# outputs
maxrotation = 0

# load the molecule
mol = scoria.Molecule(filename)

# get the coordinates
coords = mol.get_coordinates()
axisatom1coords = coords[axisatom1]
axisatom2coords = coords[axisatom2]

# break central bond
mol.delete_bond(
    axisatom1,
    axisatom2,
)

# get coords for each half
nextatoms = mol.select_all_atoms_bound_to_selection(
    np.array([axisatom1]),
)
tophalfcoords = mol.select_branch(axisatom1, nextatoms[0])

nextatoms = mol.select_all_atoms_bound_to_selection(
    np.array([axisatom2]),
)
bottomhalfcoords = mol.select_branch(axisatom2, nextatoms[0])

# split molecule into halves
bottomhalf = mol.get_molecule_from_selection(
    tophalfcoords,
)

tophalf = mol.get_molecule_from_selection(
    bottomhalfcoords,
)

# save initial orientation
copy_tophalf = tophalf.copy()

# rotate in one direction until steric clash
clash = False
rot = 0
while not clash and rot < 180:
    tophalf.rotate_molecule_around_a_line_between_points(
        axisatom1coords,
        axisatom2coords,
        math.radians(1),
    )

    # distance of top half atoms to bottom half
    res = cdist(
        tophalf.get_coordinates(),
        bottomhalf.get_coordinates(),
    )

    # ignore axis atoms
    res[0][0] = np.inf

    # check for close atoms
    clash = np.any(res < cutoff)

    rot += 1

outmol = tophalf.merge_with_another_molecule(bottomhalf)

outmol.save_pdb(
    r"output_temp.pdb",
)

print(f"Steric clash {clash} at {rot} degrees.")
maxrotation += rot

# rotate in the other direction until steric clash
clash = False
rot = 0
while not clash and rot > -180:
    copy_tophalf.rotate_molecule_around_a_line_between_points(
        axisatom1coords,
        axisatom2coords,
        math.radians(-1),
    )

    res = cdist(
        copy_tophalf.get_coordinates(),
        bottomhalf.get_coordinates(),
    )

    res[0][0] = np.inf

    clash = np.any(res < cutoff)

    rot -= 1

print(f"Steric clash {clash} at {rot} degrees.")

maxrotation -= rot

# reattach halves
outmol = copy_tophalf.merge_with_another_molecule(bottomhalf)

outmol.save_pdb(
    r"output.pdb",
)

print(f"Maximum rotation of {maxrotation} degrees.")

"""
Load molecule
split into two
rotate half a degree at a time
check for steric clashes
recombine and save
"""
