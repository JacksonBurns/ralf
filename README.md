<h1 align="center">ralf</h1> 
<h3 align="center">Rotation Angle Limit Finder</h3>

<p align="center">  
  <img alt="ralflogo" src="https://github.com/JacksonBurns/ralf/blob/main/ralf_logo.png">
</p> 
<p align="center">
  <img alt="GitHub Repo Stars" src="https://img.shields.io/github/stars/JacksonBurns/ralf?style=social">
  <img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/ralf">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/ralf">
  <img alt="PyPI - License" src="https://img.shields.io/github/license/JacksonBurns/ralf">
</p>

## Online Documentation
[Click here to read the documentation](https://JacksonBurns.github.io/ralf/)

## Usage
`ralf` includes a single function `get_rotation_limit`, which can be imported and used as shown below:
```
from ralf import get_rotation_limit

limit = get_rotation_limit(
  r"test/data/segphos.pdb",
  23,
  6,
  cutoff_distance=1.8,
)
```

Arguments are as follows:
 - pdb_path (str): PDB file for molecule.
 - axis_atom_1 (int): ID of first atom bound in chiral axis.
 - axis_atom_2 (int): ID of second atom bound in chiral axis.
 - cutoff_distance (float, optional): Distance considered overlapping, in Angstrom. Defaults to 1.6.

`get_rotation_limit` returns the maximum rotation angle in degrees.