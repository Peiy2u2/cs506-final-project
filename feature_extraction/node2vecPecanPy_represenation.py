import argparse
import networkx as nx
import pandas as pd
from pecanpy import pecanpy as node2vec


class GeneExpressionGraph:
    """A class to represent a gene expression graph where nodes are cells and genes,
    and edges are weighted by gene expression counts."""

    def __init__(self, csv_path: str):
        """Initializes the GeneExpressionGraph.

        Parameters:
        - csv_path (str): Path to the preprocessed gene expression CSV file.
        """
        self.csv_path = csv_path
        self.graph = nx.Graph()
        self._load_data()

    def _load_data(self):
        """Loads gene expression data from CSV and constructs a bipartite graph.
        - Nodes: Cells and Genes.
        - Edges: Between cells and genes, weighted by gene expression counts.
        - Edges with weight 0 are not added.
        """
        # Load data: rows = cells, columns = genes
        self.data = pd.read_csv(self.csv_path, index_col=0)

        # Add nodes (cells and genes)
        self.graph.add_nodes_from(self.data.index, type="cell")
        self.graph.add_nodes_from(self.data.columns, type="gene")

        # Add edges with weights
        for cell, row in self.data.iterrows():
            for gene, weight in row.items():
                if weight > 0:  # Ignore zero weights
                    self.graph.add_edge(cell, gene, weight=weight)

    def get_graph(self):
        """Returns the constructed graph."""
        return self.graph


class PecanPyEmbedding:
    """A class to generate node2vec embeddings using PecanPy library."""

    def __init__(
        self,
        graph: nx.Graph,
        p: float = 1.0,
        q: float = 1.0,
        dimensions: int = 128,
        walk_length: int = 40,
        num_walks: int = 10,
        workers: int = 4,
        window: int = 10,
        min_count: int = 1,
        verbose: bool = True,
    ):
        """Initializes the PecanPyEmbedding class.

        Parameters:
        - graph (nx.Graph): The gene expression graph.
        - dimensions (int): Number of embedding dimensions. Default is 128.
        - walk_length (int): Length of each random walk.
        - num_walks (int): Number of walks per node.
        - workers (int): Number of CPU cores for parallel computation.
        - window (int): Window size for Word2Vec model.
        - min_count (int): Minimum count for Word2Vec model.
        - p (float): Return parameter for Node2Vec.
        - q (float): In-out parameter for Node2Vec.
        - verbose (bool): Whether to show progress messages.
        """
        self.graph = graph
        self.dimensions = dimensions
        self.walk_length = walk_length
        self.num_walks = num_walks
        self.workers = workers
        self.window = window
        self.min_count = min_count
        self.p = p
        self.q = q
        self.verbose = verbose
        self.model = None
        self.embeddings = None

    def generate_embeddings(self):
        """Runs the PecanPy algorithm to generate embeddings."""
        # Save graph to temporary edgelist file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.edg') as f:
            # Write edges with weights
            for u, v, data in self.graph.edges(data=True):
                weight = data.get('weight', 1.0)
                f.write(f"{u}\t{v}\t{weight}\n")
            f.flush()

            # Initialize PecanPy graph object
            g = node2vec.DenseOTF(p=self.p, q=self.q, workers=self.workers, verbose=self.verbose)
            g.read_edg(f.name, weighted=True, directed=False)

            # Generate embeddings
            self.model = g.embed(
                dim=self.dimensions,
                num_walks=self.num_walks,
                walk_length=self.walk_length,
                window_size=self.window,
            )

        # Convert embeddings to DataFrame for cell nodes only
        cell_nodes = [
            str(node)
            for node, attr in self.graph.nodes(data=True)
            if attr.get("type") == "cell"
        ]
        
        # Create mapping from node names to indices
        node_id_map = {node: idx for idx, node in enumerate(g.nodes)}
        
        embeddings = {}
        for node in cell_nodes:
            if node in node_id_map:
                embeddings[node] = self.model[node_id_map[node]]

        self.embeddings = pd.DataFrame.from_dict(embeddings, orient='index')
        self.embeddings.columns = [f"embedding_dim_{i+1}" for i in range(self.dimensions)]

    def save_embeddings(self, output_path: str):
        """Saves the generated embeddings to a CSV file."""
        if self.embeddings is not None:
            self.embeddings.to_csv(output_path)
        else:
            raise ValueError("Embeddings have not been generated. Run generate_embeddings() first.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Node2Vec for gene expression graph embedding using PecanPy."
    )
    parser.add_argument(
        "--input_path",
        type=str,
        required=True,
        help="Path to input gene expression CSV file",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Path to save embeddings CSV file",
    )
    parser.add_argument(
        "--p", type=float, default=1.0, help="Return parameter for Node2Vec"
    )
    parser.add_argument(
        "--q", type=float, default=1.0, help="In-out parameter for Node2Vec"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of workers for parallel processing",
    )
    parser.add_argument(
        "--walk_length", type=int, default=40, help="Length of each random walk"
    )
    parser.add_argument(
        "--num_walks", type=int, default=10, help="Number of walks per node"
    )
    args = parser.parse_args()

    # Create the gene expression graph
    gene_graph = GeneExpressionGraph(args.input_path)
    graph = gene_graph.get_graph()

    # Generate embeddings using PecanPy
    pecanpy_embedding = PecanPyEmbedding(
        graph,
        p=args.p,
        q=args.q,
        workers=args.workers,
        walk_length=args.walk_length,
        num_walks=args.num_walks,
    )
    pecanpy_embedding.generate_embeddings()

    # Save the cell embeddings
    pecanpy_embedding.save_embeddings(args.output_path)
    print(f"Cell embeddings saved to {args.output_path}")