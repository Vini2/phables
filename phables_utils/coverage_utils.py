from collections import defaultdict
import os
import pickle
import pysam
import glob

def get_edge_coverage(coverage):

    contig_coverages = {}

    with open(coverage, "r") as myfile:
        for line in myfile.readlines():

            if not line.startswith("Contig"):
                strings = line.strip().split()

                contig_name = strings[0].replace("contig", "edge")

                coverage_sum = sum([float(x) for x in strings[1:]])

                contig_coverages[contig_name] = coverage_sum

    return contig_coverages


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    
    for read in bam.fetch(region=region_string):
        if read.is_secondary or read.is_supplementary:
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

    link_counts = defaultdict(list)

    if os.path.isfile(f"{output}/junction_pe_coverage.pickle"):
        with open(f"{output}/junction_pe_coverage.pickle", "rb") as handle:
            link_counts = pickle.load(handle)

    else:

        bam_files = glob.glob(bam_path+"/*.bam")

        for bam_file in bam_files:

            bam = pysam.AlignmentFile(bam_file, 'rb')

            read_pairs = read_pair_generator(bam)

            for read1, read2 in read_pairs:
                if read1.reference_name != read2.reference_name:
                    if read1.query_alignment_length/read1.qlen >= 0.9 and read2.query_alignment_length/read2.qlen >= 0.9:
                        if read1.reference_name=='edge_91047' and read2.reference_name=='edge_5628':
                            print(read1.reference_start, read1.reference_end, read2.reference_start, read2.reference_end)
                        link_counts[(read1.reference_name, read2.reference_name)].append(read1.query_name)

        with open(f"{output}/junction_pe_coverage.pickle", "wb") as handle:
            pickle.dump(link_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return link_counts


def get_graph_spanning_reads(gaf_path, output):

    junction_reads = defaultdict(int)

    if os.path.isfile(f"{output}/graph_spanning_reads.pickle"):
        with open(f"{output}/graph_spanning_reads.pickle", "rb") as handle:
            junction_reads = pickle.load(handle)

    else:
        
        gaf_files = glob.glob(gaf_path+"/*.gaf")

        for gaf_file in gaf_files:

            with open(gaf_file, "r") as myfile:

                for line in myfile.readlines():

                    strings = line.strip().split("\t")

                    # print(strings[0], strings[5])

                    if strings[5].count(">") == 2:
                        edges = strings[5].split(">")[1:]
                        junction_reads[(edges[0], edges[1])] += 1

                    elif strings[5].count("<") == 2:
                        edges = strings[5].split("<")[1:]
                        junction_reads[(edges[1], edges[0])] += 1

        with open(f"{output}/graph_spanning_reads.pickle", "wb") as handle:
            pickle.dump(junction_reads, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    return junction_reads