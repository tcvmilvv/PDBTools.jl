export sirah_ss

#
# Secondary structure of SIRAH coarse-grained proteins, as defined in
#
#   Darré L, Machado MR, Brandner AF, González HC, Ferreira S, Pantano S.
#   SIRAH: a structurally unbiased coarse-grained force field for proteins with
#   aqueous solvation and long-range electrostatics.
#   J. Chem. Theory Comput. 11, 723-739 (2015). DOI: 10.1021/ct5007746
#
# Implemented from the description of the method in the paper (section "Calculated
# Properties"). The backbone of each residue is represented by the beads GN, GC and GO,
# at the positions of the N, Cα and O atoms.
#

# H-bonds: beads of opposite charge within 0.4 nm. In the backbone, GO is negative
# (-0.225e) and GN and GC are positive (+0.125e and +0.10e) (Figure 2 of the paper).
const _SIRAH_HBOND_CUTOFF = 4.0 # Å

# Corrections of the CG dihedrals to the φ and ψ of the Ramachandran plot (degrees)
const _SIRAH_PHI_SHIFT = -20.0
const _SIRAH_PSI_SHIFT = 13.0

# Angle in the interval (-180, 180]
_wrap_angle(a) = a > 180 ? a - 360 : (a <= -180 ? a + 360 : a)

# Dihedral angle (degrees, IUPAC sign convention) defined by four points
function _dihedral_degrees(p1, p2, p3, p4)
    b1 = p2 .- p1
    b2 = p3 .- p2
    b3 = p4 .- p3
    cross(u, v) = (u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1])
    dot(u, v) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    y = sqrt(dot(b2, b2)) * dot(b1, n2)
    x = dot(n1, n2)
    return atan(y, x) * 180 / π
end

_bead_position(at) = (Float64(at.x), Float64(at.y), Float64(at.z))
_bead_distance(a, b) = sqrt(sum(abs2, a .- b))

# Classification from the CG Ramachandran angles (degrees), as in the paper
function _sirah_ramachandran_class(phi, psi)
    if -180 <= phi <= 10 && -120 <= psi <= 45
        return "H"
    elseif (-180 <= phi <= 0 || phi > 135) && (-180 <= psi <= -120 || 45 <= psi <= 180)
        return "E"
    else
        return "C"
    end
end

"""
    sirah_ss(atoms::AbstractVector{<:Atom})

Compute the secondary structure of a protein represented with the SIRAH coarse-grained
force field, following the definition of Darré et al. (2015). Returns a `Vector{SSData}`,
in the same format as the output of [`dssp_run`](@ref) and [`stride_run`](@ref), with one
element per protein residue and the codes `"H"` (helix), `"E"` (extended) or `"C"` (coil).

The protein residues are those with a `GC` backbone bead. The backbone beads `GN`, `GC`
and `GO`, at the positions of the N, Cα and O atoms, define the angles

    φ = dihedral(GOᵢ₋₁, GNᵢ, GCᵢ, GOᵢ) - 20°
    ψ = dihedral(GNᵢ, GCᵢ, GOᵢ, GNᵢ₊₁) + 13°

and each residue is first classified by these angles:

- H: -180° ≤ φ ≤ 10° and -120° ≤ ψ ≤ 45°;
- E: (-180° ≤ φ ≤ 0° or φ > 135°) and (-180° ≤ ψ ≤ -120° or 45° ≤ ψ ≤ 180°);
- C: otherwise.

A hydrogen bond-like interaction is a `GO` bead (negative) within 4 Å of a `GN` or `GC`
bead (positive). A residue in the H region is assigned to H only if its `GO` or `GN` bead
forms such an interaction with the residue four positions after or before it in the chain.
A residue in the E region is assigned to E only if one of its backbone beads forms such an
interaction with a residue that is also in the E region and is more than two positions
away in the chain (or in another chain). All other residues are assigned to C.

Consecutive residues are connected if they belong to the same chain and their residue
numbers differ by one. The first and last residues of each connected segment (for which φ
or ψ is not defined), and residues missing any of the backbone beads, are assigned to C;
their undefined angles are set to `360.0`.

To label the residues by their one-letter codes (e.g. in `secondary_structure_diagram`),
register the SIRAH residue names with `custom_protein_residues!(SIRAH)`.

The overall assignment is similar to that of DSSP and STRIDE (with their codes grouped
in three classes), with systematic differences that follow from the definition of the
method: β strands tend to be slightly longer, helices shorter than about four residues
are usually not detected, and isolated residues may be assigned to H. See the
documentation for a comparison with atomistic methods.

!!! note
    This is an independent implementation of the published definition of the method. The
    original implementation is the `sirah_ss` command of SIRAH Tools (Machado & Pantano,
    Bioinformatics 32, 1568, 2016), whose results may differ at chain termini.

# Example

```julia-repl
julia> using PDBTools, LazyArtifacts

julia> custom_protein_residues!(SIRAH);

julia> pdb = read_pdb(PDBTools.test_file("sirah.pdb"));

julia> ss = sirah_ss(pdb);

julia> [s.sscode for s in ss]
```
"""
function sirah_ss(atoms::AbstractVector{<:Atom})
    # Protein residues: residues with a GC bead
    residues = [r for r in eachresidue(atoms) if any(at -> name(at) == "GC", r)]
    n = length(residues)
    n > 0 || throw(ArgumentError("No residue with SIRAH backbone beads (GN, GC, GO) found."))

    function bead(r, beadname)
        for at in r
            name(at) == beadname && return _bead_position(at)
        end
        return nothing
    end
    GN = [bead(r, "GN") for r in residues]
    GC = [bead(r, "GC") for r in residues]
    GO = [bead(r, "GO") for r in residues]
    complete = [!isnothing(GN[k]) && !isnothing(GC[k]) && !isnothing(GO[k]) for k in 1:n]
    chains = [string(chain(first(r))) for r in residues]
    numbers = [Int(resnum(first(r))) for r in residues]

    # Connected segments: same chain, consecutive residue numbers, all backbone beads present
    connected_to_next = [k < n && complete[k] && complete[k+1] && chains[k] == chains[k+1] &&
                         numbers[k+1] - numbers[k] == 1 for k in 1:n]
    segment = zeros(Int, n)  # segment of each residue
    position = zeros(Int, n) # position of each residue in its segment
    for k in 1:n
        if k > 1 && connected_to_next[k-1]
            segment[k] = segment[k-1]
            position[k] = position[k-1] + 1
        else
            segment[k] = (k > 1 ? segment[k-1] : 0) + 1
            position[k] = 1
        end
    end
    # Residue at a given distance along the same segment (or nothing)
    function partner(k, offset)
        j = k + offset
        (1 <= j <= n && segment[j] == segment[k] && position[j] == position[k] + offset) ? j : nothing
    end

    # Ramachandran angles and first classification
    phi = fill(360.0, n)
    psi = fill(360.0, n)
    codes = fill("C", n)
    for k in 1:n
        complete[k] || continue
        has_previous = !isnothing(partner(k, -1))
        has_next = !isnothing(partner(k, 1))
        if has_previous
            phi[k] = _wrap_angle(_dihedral_degrees(GO[k-1], GN[k], GC[k], GO[k]) + _SIRAH_PHI_SHIFT)
        end
        if has_next
            psi[k] = _wrap_angle(_dihedral_degrees(GN[k], GC[k], GO[k], GN[k+1]) + _SIRAH_PSI_SHIFT)
        end
        if has_previous && has_next
            codes[k] = _sirah_ramachandran_class(phi[k], psi[k])
        end
    end

    hbond(acceptor, donor) = _bead_distance(acceptor, donor) <= _SIRAH_HBOND_CUTOFF

    # Helix: GO(i)⋯GN(i+4) or GN(i)⋯GO(i-4)
    helix = falses(n)
    for k in 1:n
        codes[k] == "H" || continue
        after = partner(k, 4)
        before = partner(k, -4)
        helix[k] = (!isnothing(after) && hbond(GO[k], GN[after])) ||
                   (!isnothing(before) && hbond(GO[before], GN[k]))
    end

    # Extended: GO of one residue within the cutoff of GN or GC of another, both in the
    # E region, more than two positions apart in the same segment (or in different segments)
    extended_candidates = findall(==("E"), codes)
    strand = falses(n)
    for (a, i) in enumerate(extended_candidates), j in extended_candidates[a+1:end]
        segment[i] == segment[j] && abs(position[i] - position[j]) <= 2 && continue
        if hbond(GO[i], GN[j]) || hbond(GO[i], GC[j]) || hbond(GO[j], GN[i]) || hbond(GO[j], GC[i])
            strand[i] = true
            strand[j] = true
        end
    end

    return [SSData(
        string(resname(first(residues[k]))),
        chains[k],
        numbers[k],
        helix[k] ? "H" : (strand[k] ? "E" : "C"),
        phi[k],
        psi[k],
    ) for k in 1:n]
end

@testitem "sirah_ss: helpers" begin
    using PDBTools
    @test PDBTools._wrap_angle(190.0) ≈ -170.0
    @test PDBTools._wrap_angle(-190.0) ≈ 170.0
    @test PDBTools._wrap_angle(180.0) ≈ 180.0
    # trans and cis planar arrangements
    @test abs(PDBTools._dihedral_degrees((1.0, 1.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, -1.0, 0.0))) ≈ 180.0
    @test PDBTools._dihedral_degrees((1.0, 1.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0)) ≈ 0.0 atol = 1e-10
    # sign: +90 and -90
    d1 = PDBTools._dihedral_degrees((0.0, 1.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 0.0, 1.0))
    d2 = PDBTools._dihedral_degrees((0.0, 1.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 0.0, -1.0))
    @test abs(d1) ≈ 90.0
    @test d1 ≈ -d2
    # Ramachandran regions
    @test PDBTools._sirah_ramachandran_class(-60.0, -45.0) == "H"
    @test PDBTools._sirah_ramachandran_class(-120.0, 130.0) == "E"
    @test PDBTools._sirah_ramachandran_class(150.0, -150.0) == "E"
    @test PDBTools._sirah_ramachandran_class(60.0, 30.0) == "C"
end

@testitem "sirah_ss" begin
    using PDBTools, LazyArtifacts

    # A SIRAH structure of the test files: one element per residue with a GC bead
    sirah = read_pdb(PDBTools.test_file("sirah.pdb"))
    ss = sirah_ss(sirah)
    @test ss isa Vector{SSData}
    @test length(ss) == count(at -> name(at) == "GC", sirah)
    @test all(s -> s.sscode in ("H", "E", "C"), ss)
    # termini: undefined angles, coil
    @test ss[1].phi == 360.0
    @test ss[end].psi == 360.0
    @test ss[1].sscode == "C" && ss[end].sscode == "C"
    @test_throws ArgumentError sirah_ss(sirah[1:0])

    # Pseudo-SIRAH backbone built from an atomistic structure: the GN, GC and GO beads are
    # at the positions of the N, CA and O atoms. The assignment must agree with STRIDE for
    # most residues (the paper reports ~90% of agreement).
    prot = read_pdb(PDBTools.test_file("structure.pdb"), "protein")
    cg = [deepcopy(at) for at in prot if name(at) in ("N", "CA", "O")]
    for at in cg
        at.name = Dict("N" => "GN", "CA" => "GC", "O" => "GO")[name(at)]
    end
    ss_cg = sirah_ss(cg)
    @test length(ss_cg) == length(eachresidue(prot))
    ss_stride = stride_run(prot)
    stride_class(code) = code in ("H", "G", "I") ? "H" : (code in ("E", "B") ? "E" : "C")
    stride_codes = Dict((s.chain, s.resnum) => stride_class(s.sscode) for s in ss_stride)
    matched = [stride_codes[(s.chain, s.resnum)] == s.sscode for s in ss_cg if haskey(stride_codes, (s.chain, s.resnum))]
    @test count(matched) / length(matched) > 0.75
    @test any(s -> s.sscode == "H", ss_cg)
end
