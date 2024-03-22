import logging
import os
import subprocess

FASTA_LINE_LEN = 60

# Create logger
logger = logging.getLogger("phables 1.3.2")


def write_unitigs(nodes, unitig_names, graph_unitigs, filename, output):
    """
    Write unitigs to FASTA file
    """

    with open(f"{output}/{filename}.fasta", "w+") as myfile:
        for node in nodes:
            unitig_name = unitig_names[node]
            edge_seq = str(graph_unitigs[unitig_name])
            myfile.write(f">{unitig_name}\n")

            chunks = [
                edge_seq[i : i + FASTA_LINE_LEN]
                for i in range(0, len(edge_seq), FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_component_info(all_components, output):
    """
    Write component information to file
    """

    with open(f"{output}/resolved_component_info.txt", "w") as myfile:
        myfile.write(f"Component\t")
        myfile.write(f"Number of nodes\t")
        myfile.write(f"Number of paths\t")
        myfile.write(f"Fraction of unitigs recovered\t")
        myfile.write(f"Maximum degree\t")
        myfile.write(f"Minimum degree\t")
        myfile.write(f"Maximum in degree\t")
        myfile.write(f"Maximum out degree\t")
        myfile.write(f"Average degree\t")
        myfile.write(f"Average in degree\t")
        myfile.write(f"Average out degree\t")
        myfile.write(f"Density\t")
        myfile.write(f"Maximum path length\t")
        myfile.write(f"Minimum path length\t")
        myfile.write(f"Length ratio (long/short)\t")
        myfile.write(f"Maximum coverage path length\t")
        myfile.write(f"Minimum coverage path length\t")
        myfile.write(f"Length ratio (highest cov/lowest cov)\t")
        myfile.write(f"Maximum coverage\t")
        myfile.write(f"Minimum coverage\t")
        myfile.write(f"Coverage ratio (highest/lowest)\n")

        if len(all_components) > 0:
            for component in all_components:
                myfile.write(f"{component.id}\t")
                myfile.write(f"{component.n_nodes}\t")
                myfile.write(f"{component.n_paths}\t")
                myfile.write(f"{component.frac_unitigs}\t")
                myfile.write(f"{component.max_degree}\t")
                myfile.write(f"{component.min_degree}\t")
                myfile.write(f"{component.max_in_degree}\t")
                myfile.write(f"{component.max_out_degree}\t")
                myfile.write(f"{component.avg_degree}\t")
                myfile.write(f"{component.avg_in_degree}\t")
                myfile.write(f"{component.avg_out_degree}\t")
                myfile.write(f"{component.density}\t")
                myfile.write(f"{component.max_path_length}\t")
                myfile.write(f"{component.min_path_length}\t")
                myfile.write(f"{component.min_max_len_ratio}\t")
                myfile.write(f"{component.max_cov_path_length}\t")
                myfile.write(f"{component.min_cov_path_length}\t")
                myfile.write(f"{component.min_max_cov_len_ratio}\t")
                myfile.write(f"{component.max_cov}\t")
                myfile.write(f"{component.min_cov}\t")
                myfile.write(f"{component.min_max_cov_ratio}\n")
        else:
            myfile.write(f"No complex components were resolved.")

    return "resolved_component_info.txt"


def write_res_genome_info(all_resolved_paths, output):
    """
    Write resolved genome information to file
    """

    with open(f"{output}/resolved_genome_info.txt", "w") as myfile:
        myfile.write(f"Path\tCase\tCoverage\tLength\tGC content\tNode order (gfa link format)\tNode order (human readable)\n")
        for genomic_path in all_resolved_paths:
            myfile.write(
                f"{genomic_path.id}\t{genomic_path.bubble_case}\t{genomic_path.coverage}\t{genomic_path.length}\t{genomic_path.gc}\t{genomic_path.node_order}\t{genomic_path.node_order_human}\n"
            )

    return "resolved_genome_info.txt"


def write_path(final_genomic_paths, output):
    """
    Write genomic paths to a single FASTA file
    """

    with open(f"{output}/resolved_paths.fasta", "a+") as myfile:
        for genomic_path in final_genomic_paths:
            myfile.write(f">{genomic_path.id}\n")

            chunks = [
                genomic_path.path[i : i + FASTA_LINE_LEN]
                for i in range(0, genomic_path.length, FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_path_fasta(final_genomic_paths, output_genomes_path):
    """
    Write genomic paths to individual FASTA files
    """

    if not os.path.isdir(f"{output_genomes_path}"):
        subprocess.run("mkdir -p " + output_genomes_path, shell=True)

    for genomic_path in final_genomic_paths:
        with open(f"{output_genomes_path}/{genomic_path.id}.fasta", "w+") as myfile:
            myfile.write(f">{genomic_path.id}\n")

            chunks = [
                genomic_path.path[i : i + FASTA_LINE_LEN]
                for i in range(0, genomic_path.length, FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_component_phrog_info(resolved_components, comp_phrogs, output):
    """
    Write PHROGs found in resolved components
    """

    with open(f"{output}/component_phrogs.txt", "w") as myfile:
        myfile.write(f"Phage component\tPHROG\n")
        for comp in resolved_components:
            myfile.write(f"phage_{comp}\t{comp_phrogs[comp]}\n")

    return "component_phrogs.txt"


def init_files(output):
    """
    Initialise files and folders
    """

    open(f"{output}/resolved_edges.fasta", "a").close()
    open(f"{output}/resolved_paths.fasta", "a").close()
    open(f"{output}/resolved_genome_info.txt", "a").close()
    open(f"{output}/resolved_component_info.txt", "a").close()
    open(f"{output}/component_phrogs.txt", "a").close()

    if not os.path.isdir(f"{output}/resolved_phages"):
        subprocess.run(f"mkdir -p {output}/resolved_phages", shell=True)
