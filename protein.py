import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


STANDARD_BOND_ANGLE = np.arccos(-1 / 3)

class Protein:
    def __init__(self, pdb_filename):
        """
        Create a Protein object from a PDB file.

        Parameters
        ----------
        pdb_filename: str
            The PDB filename

        Returns
        -------
        """
        # Load protein from PDB file
        self.mol = Chem.MolFromPDBFile(pdb_filename)
        # Sort atoms by residue
        self.atoms_by_residue = dict()
        for atom in self.mol.GetAtoms():
            residue_number = atom.GetPDBResidueInfo().GetResidueNumber()
            name = atom.GetPDBResidueInfo().GetName()
            if residue_number not in self.atoms_by_residue:
                self.atoms_by_residue[residue_number] = dict()
            self.atoms_by_residue[residue_number][name] = atom
        self.num_residues = max(self.atoms_by_residue)
        # Initialize conformer
        self.conf = self.mol.GetConformer()

    def get_dihedral_angle_atoms(self, residue, angle):
        """
        Get the list of atom indices defining a certain dihedral angle.

        Parameters
        ----------
        residue: int
            The residue number.
        angle: str
            A string among 'omega', 'phi', 'psi' specifying the dihedral angle in the residue.

        Returns
        -------
        tuple
            4-tuple listing the indices of the atoms defining the dihedral angle.
        """
        if angle == "omega":
            if residue > 1:
                return self.atoms_by_residue[residue - 1][" CA "].GetIdx(), self.atoms_by_residue[residue - 1][" C  "].GetIdx(), self.atoms_by_residue[residue][" N  "].GetIdx(), self.atoms_by_residue[residue][" CA "].GetIdx()
            else:
                return None
        elif angle == "phi":
            if residue > 1:
                return self.atoms_by_residue[residue - 1][" C  "].GetIdx(), self.atoms_by_residue[residue][" N  "].GetIdx(), self.atoms_by_residue[residue][" CA "].GetIdx(), self.atoms_by_residue[residue][" C  "].GetIdx()
            else:
                return None
        elif angle == "psi":
            if residue < self.num_residues:
                return self.atoms_by_residue[residue][" N  "].GetIdx(), self.atoms_by_residue[residue][" CA "].GetIdx(), self.atoms_by_residue[residue][" C  "].GetIdx(), self.atoms_by_residue[residue + 1][" N  "].GetIdx()
            else:
                return None
        else:
            raise ValueError("unknown angle '{}'".format(angle))

    def get_dihedral_angle(self, residue, angle):
        """
        Get the value of a dihedral angle.

        Parameters
        ----------
        residue: int
            The residue number.
        angle: str
            A string among 'omega', 'phi', 'psi' specifying the dihedral angle in the residue.

        Returns
        -------
        float
            The dihedral angle (in radians).
        """
        atoms = self.get_dihedral_angle_atoms(residue, angle)
        if atoms is None:
            return None
        return rdMolTransforms.GetDihedralRad(self.conf, *atoms)

    def set_dihedral_angle(self, residue, angle, value):
        """
        Set the value of a dihedral angle.

        Parameters
        ----------
        residue: int
            The residue number.
        angle: str
            A string among 'omega', 'phi', 'psi' specifying the dihedral angle in the residue.
        value: float
            The value of the dihedral angle (in radians).

        Returns
        -------
        """
        atoms = self.get_dihedral_angle_atoms(residue, angle)
        if atoms is not None:
            rdMolTransforms.SetDihedralRad(self.conf, *atoms, value)

    def set_bond_lengths(self, length=1.5):
        """
        Set the length of all bonds to a common value.

        Parameters
        ----------
        length: float
            The length of the bond (in Angstrom).

        Returns
        -------
        """
        for residue in range(1, self.num_residues + 1):
            rdMolTransforms.SetBondLength(self.conf, self.atoms_by_residue[residue][" N  "].GetIdx(), self.atoms_by_residue[residue][" CA "].GetIdx(), length)
            rdMolTransforms.SetBondLength(self.conf, self.atoms_by_residue[residue][" CA "].GetIdx(), self.atoms_by_residue[residue][" C  "].GetIdx(), length)
            if residue < self.num_residues:
                rdMolTransforms.SetBondLength(self.conf, self.atoms_by_residue[residue][" C  "].GetIdx(), self.atoms_by_residue[residue + 1][" N  "].GetIdx(), length)

    def set_bond_angles(self, angle=STANDARD_BOND_ANGLE):
        """
        Set all interbond angles to a common value.

        Parameters
        ----------
        angle: float
            The common interbond angle (in radians).

        Returns
        -------
        """
        for residue in range(1, self.num_residues + 1):
            rdMolTransforms.SetAngleRad(self.conf, self.atoms_by_residue[residue][" N  "].GetIdx(), self.atoms_by_residue[residue][" CA "].GetIdx(), self.atoms_by_residue[residue][" C  "].GetIdx(), angle)
            if residue < self.num_residues:
                rdMolTransforms.SetAngleRad(self.conf, self.atoms_by_residue[residue][" CA "].GetIdx(), self.atoms_by_residue[residue][" C  "].GetIdx(), self.atoms_by_residue[residue + 1][" N  "].GetIdx(), angle)
                rdMolTransforms.SetAngleRad(self.conf, self.atoms_by_residue[residue][" C  "].GetIdx(), self.atoms_by_residue[residue + 1][" N  "].GetIdx(), self.atoms_by_residue[residue + 1][" CA "].GetIdx(), angle)

    def reset_dihedral_angles(self):
        """
        Set all dihedral angles to pi rad.

        Parameters
        ----------

        Returns
        -------
        """
        for residue in range(1, self.num_residues + 1):
            self.set_dihedral_angle(residue, "omega", np.pi)
            self.set_dihedral_angle(residue, "phi", np.pi)
            self.set_dihedral_angle(residue, "psi", np.pi)

    def energy(self):
        """
        Evaluate the energy of the conformation from an ab-initio calculation using Psi4.

        Parameters
        ----------

        Returns
        -------
        float
            The evaluated energy.
        """
        # Convert protein to Psi4 molecule format
        psi4_geometry_strings = []
        for atom_idx, atom in enumerate(self.mol.GetAtoms()):
            if atom_idx != atom.GetIdx():
                raise Exception("wrong assumption on order of atoms returned by GetAtoms()")
            position = self.conf.GetAtomPosition(atom.GetIdx())
            psi4_geometry_strings.append("{} {} {} {}".format(atom.GetSymbol(), position.x, position.y, position.z))
        psi4_molecule = psi4.geometry("\n".join(psi4_geometry_strings))
        # Compute energy
        return psi4.energy("mp2/6-31G", molecule=psi4_molecule)

    def save(self, pdb_filename):
        """
        Save the conformation to a PDB file.

        Parameters
        ----------
        pdb_filename: str
            The PDB filename.

        Returns
        -------
        """
        Chem.MolToPDBFile(self.mol, pdb_filename)
