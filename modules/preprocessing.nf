/*
 * Preprocessing module
 *
 * Sequence extraction from input PDBs, plus a dummy-mapping helper for
 * the BUILD_CONTIGS step.  Used in the PDB-only input flow (single
 * complex via --pdb_file, or two separate PDBs via --receptor_input +
 * --effector_input).
 */


/*
 * EXTRACT_SEQUENCES
 * -----------------
 * Given a PDB file and chain IDs, extract the amino-acid sequences for
 * the receptor and effector chains.  Outputs a JSON with keys:
 *   receptor_seq, effector_seq, receptor_len, effector_len
 */
process EXTRACT_SEQUENCES {
    tag "extract_seqs"
    label 'cpu'

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path pdb_file
    val  receptor_chain
    val  effector_chain

    output:
    path "sequences.json", emit: sequences_json

    script:
    """
    cat << 'PYEOF' > extract_seqs.py
import json, sys

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEC':'U','PYL':'O',
}

def extract_chain_seq(pdb_path, chain_id):
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                ch = line[21]
                if ch != chain_id:
                    continue
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                if resnum not in residues:
                    residues[resnum] = resname
    seq = ''
    for rn in sorted(residues.keys()):
        aa = THREE_TO_ONE.get(residues[rn], 'X')
        seq += aa
    return seq, sorted(residues.keys())

rec_chain = "${receptor_chain}"
eff_chain = "${effector_chain}"

rec_seq, rec_resnums = extract_chain_seq("${pdb_file}", rec_chain)
eff_seq, eff_resnums = extract_chain_seq("${pdb_file}", eff_chain)

if not rec_seq:
    print(f"ERROR: No residues found for receptor chain {rec_chain}", file=sys.stderr)
    sys.exit(1)
if not eff_seq:
    print(f"ERROR: No residues found for effector chain {eff_chain}", file=sys.stderr)
    sys.exit(1)

result = {
    "receptor_seq": rec_seq,
    "effector_seq": eff_seq,
    "receptor_len": len(rec_seq),
    "effector_len": len(eff_seq),
    "receptor_resnums": rec_resnums,
    "effector_resnums": eff_resnums,
    "receptor_start_pdb": rec_resnums[0] if rec_resnums else 1,
}

with open("sequences.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"Receptor chain {rec_chain}: {len(rec_seq)} residues ({rec_resnums[0]}-{rec_resnums[-1]})")
print(f"Effector chain {eff_chain}: {len(eff_seq)} residues ({eff_resnums[0]}-{eff_resnums[-1]})")
PYEOF

    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \
        python extract_seqs.py
    """
}


/*
 * RESOLVE_CONTIGS
 * ---------------
 * Resolve a raw user contig string to PDB-resolved coordinates.  Used by
 * Branch B (pre-docked complex PDB), where the input PDB is the same one
 * the contig spec refers to and no HADDOCK-driven coordinate remapping
 * is needed.
 *
 * Branch A handles its own resolution + remapping inside BUILD_CONTIGS
 * because HADDOCK can renumber residues; Branch B is simpler and just
 * needs the bare contig_utils.resolve_contigs() pass.
 *
 * Wraps bin/rfdiffusion_contigs.py, which already exposes this exact
 * functionality as a CLI.
 *
 * Output: a single-line text file containing the resolved contig string,
 * consumed downstream as a value channel via .map { it.text.trim() }.
 */
process RESOLVE_CONTIGS {
    tag "resolve_contigs"
    label 'cpu'

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path pdb_file
    val  raw_contigs
    path contigs_script

    output:
    path "processed_contigs.txt", emit: resolved_contigs

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${contigs_script} \\
            --contigs "${raw_contigs}" \\
            --pdb ${pdb_file} \\
            --output processed_contigs.txt
    """
}


/*
 * WRITE_DUMMY_MAPPING
 * -------------------
 * Creates an empty JSON mapping file ({}) so BUILD_CONTIGS gets a real
 * file even when no residue remapping is needed.  load_mapping() in the
 * Python side treats files with <=2 bytes or empty JSON as "no mapping".
 */
process WRITE_DUMMY_MAPPING {
    tag "dummy_${label}"
    label 'cpu'

    input:
    val label

    output:
    path "${label}_trim_mapping.json", emit: mapping

    script:
    """
    printf '{}' > ${label}_trim_mapping.json
    """
}
