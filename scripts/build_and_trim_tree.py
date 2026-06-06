#!/usr/bin/env python3
"""
build_and_trim_tree.py
----------------------
Replicates the KBase kb_gtdbtk post-classify_wf tree-processing steps locally:

  1. Builds the `*.id_to_name-with_proximal_sp_reps.map` file from:
       - id_to_name.map  (GTDB-Tk output)
       - gtdbtk.bac120.summary.tsv (or ar53) (GTDB-Tk output)
  2. Runs trim_tree_to_target_leaves.py twice per tree:
       - First pass  → *-proximals.tree  (target leaves + proximal sp reps, no sisters)
       - Second pass → *-trimmed.tree    (+ sister context branches, + lineage file)

Usage
-----
  python3 build_and_trim_tree.py \\
      --tree          gtdbtk.bac120.classify.tree.1.tree \\
      --id_map        id_to_name.map \\
      --summary       gtdbtk.bac120.summary.tsv \\
      --trim_script   trim_tree_to_target_leaves.py \\
      --archaea_meta  /path/to/ar53_metadata_r214.tsv \\
      --bacteria_meta /path/to/bac120_metadata_r214.tsv \\
      --outdir        ./trimmed_trees

Multiple trees can be processed in one run by repeating --tree:
  python3 build_and_trim_tree.py \\
      --tree gtdbtk.backbone.bac120.classify.tree \\
      --tree gtdbtk.bac120.classify.tree.1.tree   \\
      --tree gtdbtk.bac120.classify.tree.3.tree   \\
      ...

Or use --tree_dir to auto-discover all *.tree files in a directory:
  python3 build_and_trim_tree.py \\
      --tree_dir      ./gtdbtk_output \\
      --id_map        ./gtdbtk_output/id_to_name.map \\
      --summary       ./gtdbtk_output/gtdbtk.bac120.summary.tsv \\
      --trim_script   trim_tree_to_target_leaves.py \\
      --archaea_meta  /path/to/ar53_metadata_r214.tsv \\
      --bacteria_meta /path/to/bac120_metadata_r214.tsv \\
      --outdir        ./trimmed_trees
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Build proximal-sp-rep leaf maps and trim GTDB-Tk trees locally.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input trees — at least one of --tree or --tree_dir is required
    tree_group = p.add_mutually_exclusive_group(required=True)
    tree_group.add_argument(
        "--tree", metavar="TREE", action="append", dest="trees",
        help="Newick tree file(s) from GTDB-Tk classify_wf output. "
             "Repeat the flag for multiple trees.",
    )
    tree_group.add_argument(
        "--tree_dir", metavar="DIR",
        help="Directory containing *.tree files produced by GTDB-Tk. "
             "All .tree files found will be processed.",
    )

    p.add_argument(
        "--id_map", required=True, metavar="FILE",
        help="id_to_name.map file produced by GTDB-Tk "
             "(maps internal leaf IDs like 'id0' to assembly names).",
    )
    p.add_argument(
        "--summary", required=True, metavar="FILE", action="append", dest="summaries",
        help="GTDB-Tk summary TSV file (gtdbtk.bac120.summary.tsv and/or "
             "gtdbtk.ar53.summary.tsv). Repeat for multiple files.",
    )
    p.add_argument(
        "--trim_script", required=True, metavar="FILE",
        help="Path to trim_tree_to_target_leaves.py.",
    )
    p.add_argument(
        "--archaea_meta", required=True, metavar="FILE",
        help="GTDB archaea metadata TSV (e.g. ar53_metadata_r214.tsv).",
    )
    p.add_argument(
        "--bacteria_meta", required=True, metavar="FILE",
        help="GTDB bacteria metadata TSV (e.g. bac120_metadata_r214.tsv).",
    )
    p.add_argument(
        "--outdir", default=".", metavar="DIR",
        help="Output directory for all generated files (default: current dir).",
    )
    p.add_argument(
        "--dry_run", action="store_true",
        help="Print the commands that would be run without executing them.",
    )

    args = p.parse_args()

    # Validate required files
    for label, path in [
        ("--id_map",        args.id_map),
        ("--trim_script",   args.trim_script),
        ("--archaea_meta",  args.archaea_meta),
        ("--bacteria_meta", args.bacteria_meta),
    ]:
        if not os.path.isfile(path):
            p.error(f"{label}: file not found: {path}")

    for s in args.summaries:
        if not os.path.isfile(s):
            p.error(f"--summary: file not found: {s}")

    return args


# ---------------------------------------------------------------------------
# Step 1 — Build id_to_name-with_proximal_sp_reps.map
# ---------------------------------------------------------------------------

def load_id_map(id_map_file: str) -> dict:
    """
    Read id_to_name.map → {leaf_id: assembly_name}

    The file may have 2 or 3 tab-separated columns:
      leaf_id  assembly_name  [lineage]
    Only the first two are used here.
    """
    id_map = {}
    print(f"  Loading id_to_name map from: {id_map_file}")
    with open(id_map_file) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                id_map[parts[0]] = parts[1]
    print(f"  -> {len(id_map)} query entries loaded.")
    return id_map


def load_proximal_sp_reps(summary_files: list) -> dict:
    """
    Parse GTDB-Tk summary TSV(s) and return a dict:
      {assembly_name: [sp_rep_id, ...]}

    Sources for sp rep IDs (same logic as KBase genome_obj_update.py):
      - fastani_reference
      - closest_placement_reference
      - other_related_references(genome_id,species_name,radius,ANI,AF)
    """
    sp_reps_by_query = {}

    single_fields = ["fastani_reference", "closest_placement_reference"]
    multi_field = "other_related_references(genome_id,species_name,radius,ANI,AF)"

    for summary_file in summary_files:
        print(f"  Parsing summary: {summary_file}")
        with open(summary_file) as fh:
            header = fh.readline().rstrip().split("\t")
            for line in fh:
                row = dict(zip(header, line.rstrip().split("\t")))

                # Determine the assembly/query ID column (varies between versions)
                query_id = (
                    row.get("name")
                    or row.get("user_genome")
                    or row.get("Name")
                )
                if not query_id:
                    continue

                sp_reps_by_query.setdefault(query_id, [])

                for field in single_fields:
                    val = row.get(field, "-").strip()
                    if val and val != "-":
                        if val not in sp_reps_by_query[query_id]:
                            sp_reps_by_query[query_id].append(val)

                multi_val = row.get(multi_field, "-").strip()
                if multi_val and multi_val != "-":
                    for hit in multi_val.split(";"):
                        sp_rep_id = hit.split(",")[0].strip()
                        if sp_rep_id and sp_rep_id not in sp_reps_by_query[query_id]:
                            sp_reps_by_query[query_id].append(sp_rep_id)

    total_reps = sum(len(v) for v in sp_reps_by_query.values())
    print(f"  -> {len(sp_reps_by_query)} queries with {total_reps} proximal sp rep hits.")
    return sp_reps_by_query


def get_query_ids_from_tree(tree_path: str, id_map: dict) -> list:
    """
    Return the assembly names of query leaves in the tree.
    Query leaves have IDs starting with 'id' (e.g. id0, id1...).
    """
    try:
        import ete3
    except ImportError:
        # Fallback: read the newick as text and extract id* tokens
        print("  WARNING: ete3 not importable, falling back to regex leaf extraction.")
        with open(tree_path) as fh:
            content = fh.read()
        leaf_ids = re.findall(r'\bid\d+\b', content)
        return [id_map[lid] for lid in leaf_ids if lid in id_map]

    tree = ete3.Tree(tree_path, quoted_node_names=True, format=1)
    query_names = []
    for leaf_name in tree.get_leaf_names():
        if leaf_name.startswith("id") and leaf_name in id_map:
            query_names.append(id_map[leaf_name])
    return query_names


def load_valid_metadata_accessions(bacteria_meta_path: str, archaea_meta_path: str) -> set:
    """
    Legge rapidamente la prima colonna dei file dei metadati di GTDB
    per sapere esattamente quali ID genoma sono disponibili in sicurezza.
    """
    valid_ids = set()
    for path in [bacteria_meta_path, archaea_meta_path]:
        if not path or not os.path.isfile(path):
            continue
        with open(path, 'r') as fh:
            fh.readline() # Salta l'intestazione
            for line in fh:
                parts = line.split('\t', 1)
                if parts:
                    valid_ids.add(parts[0].strip())
    return valid_ids


def build_leaflist_map(tree_path: str,
                       id_map: dict,
                       sp_reps_by_query: dict,
                       outdir: str,
                       bacteria_meta: str,  # <-- Nuovi argomenti aggiunti
                       archaea_meta: str) -> str:
    """
    Costruisce il file temporaneo .map escludendo preventivamente i genomi
    di riferimento che mancano nei metadati per evitare KeyError a valle.
    """
    print("  Pre-scansione dei file metadati per prevenire crash successivi...")
    valid_metadata_ids = load_valid_metadata_accessions(bacteria_meta, archaea_meta)

    # Determina gli sp reps rilevanti per questo specifico albero
    query_names_in_tree = get_query_ids_from_tree(tree_path, id_map)
    sp_rep_ids = []
    
    skipped_count = 0
    for qname in query_names_in_tree:
        for sp_rep in sp_reps_by_query.get(qname, []):
            if sp_rep not in sp_rep_ids:
                clean_rep = sp_rep.replace("GB_", "").replace("RS_", "")
                
                # Controllo di sicurezza stile KBase:
                if (sp_rep in valid_metadata_ids) or (clean_rep in valid_metadata_ids):
                    sp_rep_ids.append(sp_rep)
                else:
                    skipped_count += 1

    if skipped_count > 0:
        print(f"  [Protezione] Rimossi {skipped_count} genomi di riferimento assenti nei metadati locali.")

    # Genera il percorso di output
    tree_stem = re.sub(r'\.tree$', '', os.path.basename(tree_path))
    out_path = os.path.join(outdir, tree_stem + ".id_to_name-with_proximal_sp_reps.map")

    lines = []
        # Your MAG entries (leaf_id → assembly name)
    for leaf_id, assembly_name in id_map.items():
        lines.append(f"{leaf_id}\t{assembly_name}")
        # Proximal GTDB reference entries
    for sp_rep_id in sorted(sp_rep_ids):
        lines.append(f"{sp_rep_id}\t{sp_rep_id}")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"  -> Creato file leaflist map ({len(lines)} righe): {out_path}")
    return out_path

# ---------------------------------------------------------------------------
# Step 2 — Run trim_tree_to_target_leaves.py
# ---------------------------------------------------------------------------

def run_cmd(cmd: list, dry_run: bool):
    print("RUNNING: " + " ".join(cmd))
    if not dry_run:
        result = subprocess.run(cmd, check=True)
        return result.returncode
    return 0


def trim_tree(tree_path: str,
              leaflist_map: str,
              trim_script: str,
              archaea_meta: str,
              bacteria_meta: str,
              outdir: str,
              dry_run: bool):
    """
    Run trim_tree_to_target_leaves.py twice, mirroring what KBase does:
      Pass 1 → -proximals.tree  (no --sisters)
      Pass 2 → -trimmed.tree    (with --sisters + --gtdblineageoutfile)
    """
    tree_stem = re.sub(r'\.tree$', '', os.path.basename(tree_path))

    proximals_tree  = os.path.join(outdir, tree_stem + "-proximals.tree")
    trimmed_tree    = os.path.join(outdir, tree_stem + "-trimmed.tree")
    newleafnames    = os.path.join(outdir, tree_stem + "-newleafnames.map")
    lineage_out     = os.path.join(outdir, tree_stem + "-lineages.map")

    base_cmd = [
        sys.executable, trim_script,
        "--intree",               tree_path,
        "--leaflist",             leaflist_map,
        "--targetleafoutfile",    newleafnames,
        "--archaea_metadata_file", archaea_meta,
        "--bacteria_metadata_file", bacteria_meta,
    ]

    # --- Pass 1: proximals only ---
    print(f"\n  [Pass 1] proximals tree → {proximals_tree}")
    run_cmd(base_cmd + ["--outtree", proximals_tree], dry_run)

    # --- Pass 2: trimmed + sisters + lineages ---
    print(f"\n  [Pass 2] trimmed tree  → {trimmed_tree}")
    run_cmd(base_cmd + [
        "--outtree",           trimmed_tree,
        "--gtdblineageoutfile", lineage_out,
        "--sisters",
    ], dry_run)

    return {
        "proximals_tree": proximals_tree,
        "trimmed_tree":   trimmed_tree,
        "newleafnames":   newleafnames,
        "lineage":        lineage_out,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Collect tree files
    if args.tree_dir:
        tree_files = sorted(Path(args.tree_dir).glob("*.tree"))
        if not tree_files:
            sys.exit(f"ERROR: no .tree files found in {args.tree_dir}")
        tree_files = [str(t) for t in tree_files]
    else:
        tree_files = args.trees

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    print("\n=== Step 1: Loading id_map and summary tables ===")
    id_map = load_id_map(args.id_map)
    sp_reps_by_query = load_proximal_sp_reps(args.summaries)

    print(f"\n=== Processing {len(tree_files)} tree(s) ===")
    results = {}
    for tree_path in tree_files:
        if not os.path.isfile(tree_path):
            print(f"  SKIP (not found): {tree_path}")
            continue

        print(f"\n--- Tree: {tree_path} ---")

        print("\n  Step 1b: Building leaflist map ...")
        leaflist_map = build_leaflist_map(tree_path, id_map, sp_reps_by_query, args.outdir, args.bacteria_meta, args.archaea_meta)

        print("\n  Step 2: Running trim_tree_to_target_leaves.py ...")
        out_files = trim_tree(
            tree_path    = tree_path,
            leaflist_map = leaflist_map,
            trim_script  = args.trim_script,
            archaea_meta = args.archaea_meta,
            bacteria_meta= args.bacteria_meta,
            outdir       = args.outdir,
            dry_run      = args.dry_run,
        )
        results[tree_path] = out_files

    print("\n=== Done ===")
    print(f"Output files written to: {os.path.abspath(args.outdir)}")
    for tree_path, files in results.items():
        print(f"\n  {os.path.basename(tree_path)}")
        for label, path in files.items():
            exists = "✓" if os.path.isfile(path) else "✗ (dry-run)"
            print(f"    [{exists}] {label}: {os.path.basename(path)}")


if __name__ == "__main__":
    main()
