from collections import defaultdict
from pathlib import Path


def get_smg_unitigs(hmmout, mg_frac):
    """
    Get unitigs containing bacterial single-copy marker genes
    """

    # Commands
    # run_FragGeneScan.pl -genome=edges.fasta -out=edges.fasta.frag -complete=0 -train=complete -thread=8 1>edges.fasta.frag.out 2>edges.fasta.frag.err
    # hmmsearch --domtblout edges.fasta.hmmout --cut_tc --cpu 8 /home/mall0133/software/MetaCoAG/metacoag_utils/auxiliary/marker.hmm edges.fasta.frag.faa 1>edges.fasta.hmmout.out 2> edges.fasta.hmmout.err

    smg_unitigs = set()

    unitig_smgs = {}

    with open(hmmout, "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                unitig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = unitig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # unitig name
                unitig_name = "_".join(name_strings)

                if mapped_marker_length > marker_gene_length * mg_frac:
                    smg_unitigs.add(unitig_name)

                    if unitig_name not in unitig_smgs:
                        unitig_smgs[unitig_name] = set()
                        unitig_smgs[unitig_name].add(marker_gene)
                    else:
                        unitig_smgs[unitig_name].add(marker_gene)

    return smg_unitigs


def get_phrog_unitigs(phrogs, e_value, seq_identity):
    """
    Get unitigs containing PHROGs
    """

    # Read phrogs table and get annotations and categories
    phrog_table_file = Path(__file__).parent / "phrogs" / "phrog_annot.tsv"

    phrog_dict = defaultdict(str)

    with open(phrog_table_file, "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("phrog"):
                strings = line.strip().split("\t")
                phrog_dict[f"phrog_{strings[0]}"] = f"{strings[2]} {strings [3]}"

    # Get unitigs containing phrogs
    unitig_phrogs = {}

    with open(phrogs, "r") as myfile:
        for line in myfile.readlines():
            # if "edge_" in line:

            strings = line.strip().split("\t")

            name = strings[0][1:-1]
            phrog_id = strings[1][1:-1].split()[0]
            seqIdentity = float(strings[3])
            evalue = float(strings[4])

            if evalue < e_value and seqIdentity > seq_identity:
                if name not in unitig_phrogs:
                    unitig_phrogs[name] = set([phrog_id])
                else:
                    unitig_phrogs[name].add(phrog_id)

    return unitig_phrogs, phrog_dict
