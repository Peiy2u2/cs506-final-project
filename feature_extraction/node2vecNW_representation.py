import argparse

import networkx as nx
import pandas as pd
from node2vec import Node2Vec


class GeneExpressionGraph:
    """
    A class to represent a gene expression graph where nodes are cells and genes,
    and edges are weighted by gene expression counts.
    """

    def __init__(self, csv_path: str):
        """
        Initializes the GeneExpressionGraph.

        Parameters:
        - csv_path (str): Path to the preprocessed gene expression CSV file.
        """
        self.csv_path = csv_path
        self.graph = nx.Graph()
        self._load_data()

    def _load_data(self):
        """
        Loads gene expression data from CSV and constructs a bipartite graph.
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


class Node2vecEmbedding:
    """
    A class to generate 128-dimensional Node2vec embeddings for cell nodes in a gene expression graph.
    """

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
        batch_words: int = 4,
    ):
        """
        Initializes the Node2vecEmbedding class.

        Parameters:
        - graph (nx.Graph): The gene expression graph.
        - dimensions (int): Number of embedding dimensions. Default is 128.
        - walk_length (int): Length of each random walk.
        - num_walks (int): Number of walks per node.
        - workers (int): Number of CPU cores for parallel computation.
        - window (int): Window size for Word2Vec model.
        - min_count (int): Minimum count for Word2Vec model.
        - batch_words (int): Batch size for Word2Vec model.
        - p (float): Return parameter for Node2Vec (controls likelihood of revisiting nodes).
        - q (float): In-out parameter for Node2Vec (controls exploration vs. homophily).
        """
        self.graph = graph
        self.dimensions = dimensions
        self.walk_length = walk_length
        self.num_walks = num_walks
        self.workers = workers
        self.window = window
        self.min_count = min_count
        self.batch_words = batch_words
        self.p = p
        self.q = q
        self.model = None
        self.embeddings = None

    def generate_embeddings(self):
        """
        Runs the  algorithm to generate embeddings for cell nodes.
        """
        # Run Node2Vec
        node2vec = Node2Vec(
            self.graph,
            dimensions=self.dimensions,
            walk_length=self.walk_length,
            num_walks=self.num_walks,
            workers=self.workers,
            p=self.p,
            q=self.q,
        )
        self.model = node2vec.fit(
            window=self.window, min_count=self.min_count, batch_words=self.batch_words
        )

        # Extract embeddings for cell nodes only
        cell_nodes = [
            node
            for node, attr in self.graph.nodes(data=True)
            if attr.get("type") == "cell"
        ]
        embeddings = {cell: self.model.wv[cell] for cell in cell_nodes}

        # Convert to DataFrame
        self.embeddings = pd.DataFrame.from_dict(embeddings, orient="index")
        self.embeddings.columns = [
            f"embedding_dimension_{i+1}" for i in range(self.dimensions)
        ]

    def save_embeddings(self, output_path: str):
        """
        Saves the generated embeddings to a CSV file.

        Parameters:
        - output_path (str): Path to save the embeddings CSV file.
        """
        if self.embeddings is not None:
            self.embeddings.to_csv(output_path)
        else:
            raise ValueError(
                "Embeddings have not been generated. Run generate_embeddings() first."
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Node2Vec for gene expression graph embedding."
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

    # Step 2: Run Node2Vec to get embeddings
    deepwalk = Node2vecEmbedding(
        graph,
        p=args.p,
        q=args.q,
        workers=args.workers,
        walk_length=args.walk_length,
        num_walks=args.num_walks,
    )
    deepwalk.generate_embeddings()

    # Step 3: Save the cell embeddings
    deepwalk.save_embeddings(args.output_path)
    print(f"Cell embeddings saved to {args.output_path}")
