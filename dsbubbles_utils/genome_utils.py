class GenomePath:
    def __init__(self, id, node_order, path, coverage, length):
        self.id = id
        self.path = path
        self.coverage = coverage
        self.length = length
        self.node_order = node_order
