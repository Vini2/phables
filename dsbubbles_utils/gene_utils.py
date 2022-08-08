def get_smg_contigs(hmmout, mg_frac):

    # Commands
    # run_FragGeneScan.pl -genome=edges.fasta -out=edges.fasta.frag -complete=0 -train=complete -thread=8 1>edges.fasta.frag.out 2>edges.fasta.frag.err
    # hmmsearch --domtblout edges.fasta.hmmout --cut_tc --cpu 8 /home/mall0133/software/MetaCoAG/auxiliary/marker.hmm edges.fasta.frag.faa 1>edges.fasta.hmmout.out 2> edges.fasta.hmmout.err

    smg_contigs = set()

    contig_smgs = {}

    with open(hmmout, "r") as myfile:

        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = contig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # Contig name
                contig_name = "_".join(name_strings)

                if mapped_marker_length > marker_gene_length * mg_frac:
                    smg_contigs.add(contig_name)

                    if contig_name not in contig_smgs:
                        contig_smgs[contig_name] = set()
                        contig_smgs[contig_name].add(marker_gene)
                    else:
                        contig_smgs[contig_name].add(marker_gene)

    return smg_contigs, contig_smgs


def get_phrog_contigs(phrogs, align_score, seq_identity):

    contig_phrogs = {}

    with open(phrogs, "r") as myfile:

        for line in myfile.readlines():

            if line.startswith("contig_"):

                strings = line.strip().split()

                name = strings[0].replace("contig", "edge")
                phrog = " ".join(strings[12:])
                alnScore = float(strings[2])
                seqIdentity = float(strings[3])

                if alnScore > align_score and seqIdentity > seq_identity:

                    if name not in contig_phrogs:
                        contig_phrogs[name] = set([phrog])
                    else:
                        contig_phrogs[name].add(phrog)

    return contig_phrogs
