"""
Berta Bori-Bru
IQAC - CSIC
Winter 2026

Code extracted from AF3 Repo. 
Only the necessary to make a CCDmmCIF file from a smiles or mol2 or whatever
Added support to generate the mmCIF file
"""


import collections
from collections.abc import Mapping, Sequence

import numpy as np
import rdkit.Chem as rd_chem
from rdkit.Chem import AllChem

import sys



class UnsupportedMolBondError(Exception):
  """Raised when we try to handle unsupported RDKit bonds."""


_RDKIT_MMCIF_TO_BOND_TYPE: Mapping[str, rd_chem.BondType] = {
    'SING': rd_chem.BondType.SINGLE,
    'DOUB': rd_chem.BondType.DOUBLE,
    'TRIP': rd_chem.BondType.TRIPLE,
}

_RDKIT_BOND_TYPE_TO_MMCIF: Mapping[rd_chem.BondType, str] = {
    v: k for k, v in _RDKIT_MMCIF_TO_BOND_TYPE.items()
}

_RDKIT_BOND_STEREO_TO_MMCIF: Mapping[rd_chem.BondStereo, str] = {
    rd_chem.BondStereo.STEREONONE: 'N',
    rd_chem.BondStereo.STEREOE: 'E',
    rd_chem.BondStereo.STEREOZ: 'Z',
    rd_chem.BondStereo.STEREOCIS: 'Z',
    rd_chem.BondStereo.STEREOTRANS: 'E',
}



def assign_atom_names_from_graph(
    mol: rd_chem.Mol,
    keep_existing_names: bool = False,
) -> rd_chem.Mol:
  """Assigns atom names from the molecular graph.

  The atom name is stored as an atom property 'atom_name', accessible
  with atom.GetProp('atom_name'). If the property is already specified, and
  keep_existing_names is True we keep the original name.

  We traverse the graph in the order of the rdkit atom index and give each atom
  a name equal to '{ELEMENT_TYPE}{INDEX}'. E.g. C5 is the name for the fifth
  unnamed carbon encountered.

  NOTE: A new mol is returned, the original is not changed in place.

  Args:
    mol: Mol object.
    keep_existing_names: If True, atoms that already have the atom_name property
      will keep their assigned names.

  Returns:
    A new mol, with potentially new 'atom_name' properties.
  """
  mol = rd_chem.Mol(mol)

  specified_atom_names = {
      atom.GetProp('atom_name')
      for atom in mol.GetAtoms()
      if atom.HasProp('atom_name') and keep_existing_names
  }

  element_counts = collections.Counter()
  for atom in mol.GetAtoms():
    if not atom.HasProp('atom_name') or not keep_existing_names:
      element = atom.GetSymbol()
      while True:
        element_counts[element] += 1
        # Standardize names by using uppercase element type, as in CCD. Only
        # effects elements with more than one letter, e.g. 'Cl' becomes 'CL'.
        new_name = f'{element.upper()}{element_counts[element]}'
        if new_name not in specified_atom_names:
          break
      atom.SetProp('atom_name', new_name)

  return mol

def mol_to_ccd_cif(
    mol: rd_chem.Mol,
    component_id: str,
    pdbx_smiles: str | None = None,
    include_hydrogens: bool = True,
): #-> cif_dict.CifDict:
  
  """Creates a CCD-like mmcif data block from an rdkit Mol object.

  Only a subset of associated mmcif fields is populated, but that is
  sufficient for further usage, e.g. in featurization code.

  Atom names can be specified via `atom_name` property. For atoms with
  unspecified value of that property, the name is assigned based on element type
  and the order in the Mol object.

  If the Mol object has associated conformers, atom positions from the first of
  them will be populated in the resulting mmcif file.

  Args:
     mol: An rdkit molecule.
     component_id: Name of the molecule to use in the resulting mmcif. That is
       equivalent to CCD code.
     pdbx_smiles: If specified, the value will be used to populate
       `_chem_comp.pdbx_smiles`.
     include_hydrogens: Whether to include atom and bond data involving
       hydrogens.

  Returns:
     An mmcif data block corresponding for the given rdkit molecule.

  Raises:
    UnsupportedMolBond: When a molecule contains a bond that can't be
      represented with mmcif.
  """
  mol = rd_chem.Mol(mol)
  if include_hydrogens:
    mol = rd_chem.AddHs(mol)
  rd_chem.Kekulize(mol)

  # BERTA: make conformer in 3D and minimize
  AllChem.EmbedMolecule(mol, AllChem.ETKDG())
  AllChem.UFFOptimizeMolecule(mol)

  if mol.GetNumConformers() > 0:
    ideal_conformer = mol.GetConformer(0).GetPositions()
    ideal_conformer = np.vectorize(lambda x: f'{x:.3f}')(ideal_conformer)
  else:
    # No data will be populated in the resulting mmcif if the molecule doesn't
    # have any conformers attached to it.
    ideal_conformer = None

  mol_cif = collections.defaultdict(list)
  mol_cif['data_'] = component_id
  mol_cif['_chem_comp.id'] = component_id
  
  # BERTA: Adding other mandatory fields
  mol_cif["_chem_comp.name"]="?"
  mol_cif["_chem_comp.type"]="non-polymer"
  mol_cif['_chem_comp.formula']=f"'{rd_chem.rdMolDescriptors.CalcMolFormula(mol)}'"
  mol_cif['_chem_comp.mon_nstd_parent_comp_id']="?"
  mol_cif['_chem_comp.pdbx_synonyms']="?"
  mol_cif['_chem_comp.formula_weight']=f"{rd_chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}"

  if pdbx_smiles:
    mol_cif['_chem_comp.pdbx_smiles'] = f"'{pdbx_smiles}'"

  mol = assign_atom_names_from_graph(mol, keep_existing_names=True)

  for atom_idx, atom in enumerate(mol.GetAtoms()):
    element = atom.GetSymbol()
    if not include_hydrogens and element in ('H', 'D'):
      continue

    mol_cif['_chem_comp_atom.comp_id'].append(component_id)
    mol_cif['_chem_comp_atom.atom_id'].append(atom.GetProp('atom_name'))
    mol_cif['_chem_comp_atom.type_symbol'].append(atom.GetSymbol().upper())
    mol_cif['_chem_comp_atom.charge'].append(str(atom.GetFormalCharge()))
    mol_cif['_chem_comp_atom.pdbx_leaving_atom_flag'].append("N") # BERTA: hard coding because i don't think there will be leaving groups in predictions
    if ideal_conformer is not None:
      coords = ideal_conformer[atom_idx]
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_x_ideal'].append(coords[0])
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_y_ideal'].append(coords[1])
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_z_ideal'].append(coords[2])

  for bond in mol.GetBonds():
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    if not include_hydrogens and (
        atom1.GetSymbol() in ('H', 'D') or atom2.GetSymbol() in ('H', 'D')
    ):
      continue
    mol_cif['_chem_comp_bond.comp_id'].append(component_id)
    mol_cif['_chem_comp_bond.atom_id_1'].append(
        bond.GetBeginAtom().GetProp('atom_name')
    )
    mol_cif['_chem_comp_bond.atom_id_2'].append(
        bond.GetEndAtom().GetProp('atom_name')
    )
    try:
      bond_type = bond.GetBondType()
      # Older versions of RDKit did not have a DATIVE bond type. Convert it to
      # SINGLE to match the AF3 training setup.
      if bond_type == rd_chem.BondType.DATIVE:
        bond_type = rd_chem.BondType.SINGLE
      mol_cif['_chem_comp_bond.value_order'].append(
          _RDKIT_BOND_TYPE_TO_MMCIF[bond_type]
      )
      mol_cif['_chem_comp_bond.pdbx_stereo_config'].append(
          _RDKIT_BOND_STEREO_TO_MMCIF[bond.GetStereo()]
      )
    except KeyError as e:
      raise UnsupportedMolBondError from e
    mol_cif['_chem_comp_bond.pdbx_aromatic_flag'].append(
        'Y' if bond.GetIsAromatic() else 'N'
    )

  return mol_cif


def make_CCDmmCIF(properties_dict):
  # Classify flags
  SINGLE_FLAGS    = []
  ATOM_FLAGS      = []
  ATOM_DATA       = []
  BOND_FLAGS      = []
  BOND_DATA       = []

  for flag in properties_dict.keys():
      if "_chem_comp_atom." in flag:
          ATOM_FLAGS.append(flag)
          ATOM_DATA.append(properties_dict[flag])
      elif "_chem_comp_bond." in flag:
          BOND_FLAGS.append(flag)
          BOND_DATA.append(properties_dict[flag])
      elif "_chem_comp." in flag:
          SINGLE_FLAGS.append(flag)

  lines = []

  # Write singles
  lines.append(f"data_{properties_dict["data_"]}\n")
  lines.append(f"#\n")

  for flag in SINGLE_FLAGS:
      lines.append(f"{flag} {properties_dict[flag]}\n")

  lines.append(f"#\n")

  # Write atom data
  lines.append(f"loop_\n")
  lines.extend([f"{flag}\n" for flag in ATOM_FLAGS])

  for fields in zip(*ATOM_DATA): #itera
      lines.append(" ".join(fields)+"\n")

  lines.append(f"#\n")

  # Write bond data
  lines.append(f"loop_\n")
  lines.extend([f"{flag}\n" for flag in BOND_FLAGS])

  for fields in zip(*BOND_DATA): #itera
      lines.append(" ".join(fields)+"\n")

  lines.append(f"#\n")

  # write output
  with open(f"{properties_dict["data_"]}.cif","w") as ff:
      ff.writelines(lines)

def main(component_id:str,smiles:str):
  mol             = rd_chem.MolFromSmiles(smiles)

  properties_dict = mol_to_ccd_cif(mol=mol,
                  component_id=component_id,
                  pdbx_smiles=smiles,
                  )

  make_CCDmmCIF(properties_dict)


if __name__=="__main__":
    component_id    = sys.argv[1]
    smiles          = sys.argv[2]

    main(component_id,smiles)

