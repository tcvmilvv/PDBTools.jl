```@meta
CollapsedDocStrings = true
```

# [Secondary structure](@id secondary-structure)

These functions provide an interface to compute the secondary structure assignment of proteins
using the STRIDE and DSSP algorithms, using as inputs vectors of `PDBTools.Atom`s.

```@docs
stride_run
dssp_run
```

The `stride_run` and `dssp_run` functions return a vector of `SSData` objects, each containing the secondary structure assignment and backbone dihedral angles for each residue.

These functions return a vector of `SSData` objects, as defined in [ProteinSecondaryStructures.jl](https://biojulia.dev/ProteinSecondaryStructures.jl). The secondary structure assignment codes are available 
in the [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/overview/#Secondary-structure-classes).

!!! note
    Non-standard residue or atom names may lead to incorrect secondary structure assignments. To ensure accurate results, it is recommended to replace non-standard names with their standard three-letter codes before running these functions, and to verify that all backbone atoms (N, CA, C, O) are present in each residue. Errors or warnings will be issued if residue names are not recognized or backbone atoms are missing.

## Example using STRIDE

The `stride_run` function runs the STRIDE algorithm on the provided array of atoms.
```@example secondary-structure
using PDBTools, LazyArtifacts
atoms = read_pdb(PDBTools.test_file("structure.pdb"), "protein)
ss = stride_run(atoms)
```

To run with DSSP just use the `dssp_run` function instead. 

## Coarse-grained structures (SIRAH)

`dssp_run` and `stride_run` need the atoms of the backbone, so they cannot be used with
coarse-grained structures. For proteins represented with the
[SIRAH](https://sirahff.github.io) force field, [`sirah_ss`](@ref) computes the secondary
structure from the backbone beads `GN`, `GC` and `GO` (at the positions of the N, Cα and O
atoms), following the definition of Darré et al. (J. Chem. Theory Comput. 11, 723, 2015).
Each residue is classified as helix (`"H"`), extended (`"E"`) or coil (`"C"`) from the
dihedral angles between the backbone beads and from hydrogen bond-like interactions
between them.

The output has the same format as that of `dssp_run` and `stride_run`, so it can be used
in the same way, e.g. in secondary structure diagrams:

```@example sirah-ss
using PDBTools, LazyArtifacts
custom_protein_residues!(SIRAH) # one-letter codes of the SIRAH residue names
pdb = read_pdb(PDBTools.test_file("sirah.pdb"))
ss = sirah_ss(select(pdb, "protein"))
```

### Comparison with atomistic methods

The classification has only three classes, so to compare it with DSSP or STRIDE their
codes must be grouped: H, G and I as `"H"`; E and B as `"E"`; all others as `"C"`. Since the
SIRAH beads of the backbone are at the positions of the N, Cα and O atoms, an atomistic
structure can be compared directly with DSSP and STRIDE, by renaming these atoms (the
agreement between DSSP and STRIDE is shown as a reference):

```@example sirah-ss
prot = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
cg = [deepcopy(at) for at in prot if name(at) in ("N", "CA", "O")]
for at in cg
    at.name = Dict("N" => "GN", "CA" => "GC", "O" => "GO")[name(at)]
end
three_class(code) = code in ("H", "G", "I") ? "H" : (code in ("E", "B") ? "E" : "C")
dssp = Dict((s.chain, s.resnum) => three_class(s.sscode) for s in dssp_run(prot))
stride = Dict((s.chain, s.resnum) => three_class(s.sscode) for s in stride_run(prot))
ss_cg = sirah_ss(cg)
# residues present in the three results
common = [s for s in ss_cg if haskey(dssp, (s.chain, s.resnum)) && haskey(stride, (s.chain, s.resnum))]
agreement(a, b) = round(100 * count(s -> a(s) == b(s), common) / length(common); digits = 1)
cg_code(s) = s.sscode
dssp_code(s) = dssp[(s.chain, s.resnum)]
stride_code(s) = stride[(s.chain, s.resnum)]
(sirah_vs_dssp = agreement(cg_code, dssp_code),
 sirah_vs_stride = agreement(cg_code, stride_code),
 dssp_vs_stride = agreement(dssp_code, stride_code))
```

The overall assignment is similar to that of DSSP and STRIDE (the original authors report
~90% of agreement with STRIDE on a large set of proteins), with some systematic
differences that follow from the definition of the method:

- β strands tend to be slightly longer, and may bridge loops of one or two residues;
- helices shorter than about four residues are usually not detected, as the H class
  requires hydrogen bond-like interactions with the residues four positions apart;
- isolated residues may be assigned to H, as there is no minimum length of helices;
- in low-resolution structures (e.g. cryo-EM models), helices with distorted hydrogen
  bonds may be detected by `sirah_ss` and DSSP but assigned as turns by STRIDE.

### Requirements

- Consecutive residues are considered connected if they belong to the same chain and
  their residue numbers differ by one. Check the chain identifiers and the numbering of
  structures written by simulation packages (e.g. chains missing or residue numbers
  restarting within a chain).
- The first and last residues of each chain, and residues missing any of the backbone
  beads, are assigned to `"C"`. Results may thus differ from those of the original
  `sirah_ss` command of SIRAH Tools at chain termini.
- Select the protein (e.g. `select(pdb, "protein")`, after `custom_protein_residues!(SIRAH)`)
  to exclude water, ions and lipids. For very large systems (e.g. whole viral capsids),
  the identification of extended residues compares all pairs of candidates, which may be
  slow: compute the secondary structure of one or a few chains at a time.

```@docs
sirah_ss
```

## Utility functions 

`PDBTools` also reexports the `ss_composition` , `ss_name`, `ss_code`, `ss_number` of `ProteinSecondaryStructures.jl`, which can be used to analyze the secondary structure assignment results. For example, to compute the secondary structure composition from the `ss` vector obtained above:

```@example secondary-structure
ss_composition(ss)
```

```@example secondary-structure
ss_name(ss[1]) # name of the secondary structure of the first residue
```

For further information, refer to the [ProteinSecondaryStructures.jl documentation](https://biojulia.dev/ProteinSecondaryStructures.jl/stable/user_guide/#Retrieving-names,-codes,-and-numeric-codes).




