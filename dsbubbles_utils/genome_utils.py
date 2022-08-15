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
        max_path_length = 0,
        min_path_length = 0,
        min_max_len_ratio = 0,
        max_cov_path_length = 0,
        min_cov_path_length = 0,
        min_max_cov_len_ratio = 0,
        max_cov = 0,
        min_cov = 0,
        min_max_cov_ratio = 0
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
        self.max_path_length = max_path_length,
        self.min_path_length = min_path_length,
        self.min_max_len_ratio = min_max_len_ratio,
        self.max_cov_path_length = max_cov_path_length,
        self.min_cov_path_length = min_cov_path_length,
        self.min_max_cov_len_ratio = min_max_cov_len_ratio,
        self.max_cov = max_cov,
        self.min_cov = min_cov,
        self.min_max_cov_ratio = min_max_cov_ratio

