def get_circular(paths):

    circular = {}

    with open(paths, "r") as myfile:

        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                if strings[3] == "Y":
                    contig_name = strings[0].replace("contig", "edge")
                    contig_length = int(strings[1])
                    circular[contig_name] = contig_length

    return circular
