import os
import math

import scoria

# input args
filename = os.path.join(os.getcwd(), "test", "data", "segphos.pdb")
axisatom1 = 5
axisatom2 = 22
dihedralatom1 = 7
dihedralatom2 = 24

# load the molecule
mol = scoria.Molecule(filename)

# get all the coordinates
coords = mol.get_coordinates()
dihedralatom1coords = coords[dihedralatom1]
axisatom1coords = coords[axisatom1]
axisatom2coords = coords[axisatom2]
dihedralatom2coords = coords[dihedralatom2]

# baseline dihedral angle
angle = mol.get_dihedral_angle(
    dihedralatom1coords,
    axisatom1coords,
    axisatom2coords,
    dihedralatom2coords,
)

print(angle, "radians")
print(angle * 180 / math.pi, "degrees")

# break central bond
mol.delete_bond(
    axisatom1coords,
    axisatom2coords,
)

# get a half of the molecule
molhalf = mol.get_molecule_from_selection()
otherhalf

steric_clash_with_another_molecule(other_mol, cutoff, pairwise_comparison=True)


# add bond back
molhalf.mermerge_with_another_molecule(otherhalf)

mol.save_pdb(
    r"output.pdb",
)
