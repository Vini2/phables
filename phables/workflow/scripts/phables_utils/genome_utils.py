# Class for genome path
class GenomePath:
    def __init__(
        self, id, bubble_case, node_order, node_order_human, node_id_order, path, coverage, length, gc
    ):
        self.id = id
        self.bubble_case = bubble_case
        self.path = path
        self.coverage = coverage
        self.length = length
        self.node_order = node_order
        self.node_order_human = node_order_human
        self.node_id_order = node_id_order
        self.gc = gc


# Class for genome component
class GenomeComponent:
    def __init__(
        self,
        id,
        n_nodes,
        n_paths,
        max_degree,
        min_degree,
        max_in_degree,
        max_out_degree,
        avg_degree,
        avg_in_degree,
        avg_out_degree,
        density,
        max_path_length,
        min_path_length,
        min_max_len_ratio,
        max_cov_path_length,
        min_cov_path_length,
        min_max_cov_len_ratio,
        max_cov,
        min_cov,
        min_max_cov_ratio,
        frac_unitigs,
    ):
        self.id = id
        self.n_nodes = n_nodes
        self.n_paths = n_paths
        self.max_degree = max_degree
        self.min_degree = min_degree
        self.max_in_degree = max_in_degree
        self.max_out_degree = max_out_degree
        self.avg_degree = avg_degree
        self.avg_in_degree = avg_in_degree
        self.avg_out_degree = avg_out_degree
        self.density = density
        self.max_path_length = max_path_length
        self.min_path_length = min_path_length
        self.min_max_len_ratio = min_max_len_ratio
        self.max_cov_path_length = max_cov_path_length
        self.min_cov_path_length = min_cov_path_length
        self.min_max_cov_len_ratio = min_max_cov_len_ratio
        self.max_cov = max_cov
        self.min_cov = min_cov
        self.min_max_cov_ratio = min_max_cov_ratio
        self.frac_unitigs = frac_unitigs
