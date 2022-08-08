class GenomePath:
    def __init__(self, id, node_order, path, coverage, length):
        self.id = id
        self.path = path
        self.coverage = coverage
        self.length = length
        self.node_order = node_order


class GenomeComponent:
    def __init__(
        self,
        id,
        n_nodes,
        max_degree,
        max_in_degree,
        max_out_degree,
        avg_degree,
        avg_in_degree,
        avg_out_degree,
        density,
    ):
        self.id = id
        self.n_nodes = n_nodes
        self.max_degree = max_degree
        self.max_in_degree = max_in_degree
        self.max_out_degree = max_out_degree
        self.avg_degree = avg_degree
        self.avg_in_degree = avg_in_degree
        self.avg_out_degree = avg_out_degree
        self.density = density
