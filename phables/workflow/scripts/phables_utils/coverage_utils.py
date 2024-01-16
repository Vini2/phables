import gc
import glob
import os
import pickle
from collections import defaultdict

import pysam


def get_unitig_coverage(coverage):
    """
    Get coverage values of unitigs
    """

    unitig_coverages = {}

    with open(coverage, "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("Contig"):
                strings = line.strip().split()

                unitig_name = strings[0]

                coverage_sum = sum([float(x) for x in strings[1:]])

                unitig_coverages[unitig_name] = coverage_sum

    return unitig_coverages


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])

    for read in bam.fetch(region=region_string):
        if (
            read.is_secondary
            or read.is_supplementary
            or not read.is_paired
            or read.mapping_quality <= 30
        ):
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

    return read_dict


def get_junction_pe_coverage(bam_path, output):
    """
    Get number of paired end reads supporting a junction
    """

    link_counts = defaultdict(int)

    if os.path.isfile(f"{output}/junction_pe_coverage.pickle"):
        with open(f"{output}/junction_pe_coverage.pickle", "rb") as handle:
            link_counts = pickle.load(handle)

    else:
        bam_files = glob.glob(bam_path + "/*.bam")

        for bam_file in bam_files:
            bam = pysam.AlignmentFile(bam_file, "rb")

            read_pairs = read_pair_generator(bam)

            for read1, read2 in read_pairs:
                if read1.reference_name != read2.reference_name:
                    link_counts[(read1.reference_name, read2.reference_name)] += 1

        with open(f"{output}/junction_pe_coverage.pickle", "wb") as handle:
            pickle.dump(link_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return link_counts


def get_sub_path_coverage(sub_path_cov, bam_path, output):
    """
    Get coverage values of sub paths using long reads
    """

    if os.path.isfile(f"{output}/sub_path_coverage.pickle"):
        sub_path_cov = defaultdict(int)
        with open(f"{output}/sub_path_coverage.pickle", "rb") as handle:
            sub_path_cov = pickle.load(handle)

    else:
        bam_files = glob.glob(bam_path + "/*.bam")

        for bam_file in bam_files:
            unitig_reads = defaultdict(set)

            bam = pysam.AlignmentFile(bam_file, "rb")

            for read in bam:
                if not read.is_unmapped:  # Only consider mapped reads
                    query_id = read.query_name
                    target_id = read.reference_name
                    unitig_reads[target_id].add(query_id)

            for sub_path in sub_path_cov.keys():
                if len(sub_path) == 3:
                    node1 = sub_path[0]
                    node2 = sub_path[1]
                    node3 = sub_path[2]

                    intersection_set = unitig_reads[node1].intersection(
                        unitig_reads[node2], unitig_reads[node3]
                    )
                    sub_path_cov[sub_path] += len(intersection_set)

                elif len(sub_path) == 2:
                    node1 = sub_path[0]
                    node2 = sub_path[1]

                    intersection_set = unitig_reads[node1].intersection(
                        unitig_reads[node2]
                    )
                    sub_path_cov[sub_path] += len(intersection_set)

            del unitig_reads
            del bam
            gc.collect()

    return sub_path_cov


def get_graph_spanning_reads(gaf_path, output):
    """
    Get number of reads spanning across a junction
    """

    junction_reads = defaultdict(int)

    if os.path.isfile(f"{output}/graph_spanning_reads.pickle"):
        with open(f"{output}/graph_spanning_reads.pickle", "rb") as handle:
            junction_reads = pickle.load(handle)

    else:
        gaf_files = glob.glob(gaf_path + "/*.gaf")

        for gaf_file in gaf_files:
            with open(gaf_file, "r") as myfile:
                for line in myfile.readlines():
                    strings = line.strip().split("\t")

                    if strings[5].count(">") == 2:
                        edges = strings[5].split(">")[1:]
                        junction_reads[(edges[0], edges[1])] += 1

                    elif strings[5].count("<") == 2:
                        edges = strings[5].split("<")[1:]
                        junction_reads[(edges[1], edges[0])] += 1

        with open(f"{output}/graph_spanning_reads.pickle", "wb") as handle:
            pickle.dump(junction_reads, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return junction_reads
