import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class PCAProcessor:
    """
    A class to perform PCA on gene expression data.
    """

    def __init__(self, input_path, variance_threshold, output_path, output_plot):
        """
        Initializes the PCAProcessor.

        Parameters:
        - input_path (str): Path to the preprocessed gene expression CSV file.
        - variance_threshold (float): Minimum variance to be explained by selected principal components.
        - output_path (str): Path to save PCA results.
        - output_plot (str): Path to save cumulative variance plot.
        """
        self.input_path = input_path
        self.variance_threshold = variance_threshold
        self.output_path = output_path
        self.output_plot = output_plot
        self.df = None
        self.pca_result = None
        self.k = None

    def load_data(self):
        """Loads gene expression data from a CSV file."""
        self.df = pd.read_csv(self.input_path, index_col=0)

    def apply_pca(self):
        """Applies PCA to the standardized gene expression data."""
        scaler = StandardScaler()
        df_scaled = scaler.fit_transform(self.df.iloc[:, :-1])

        pca = PCA()
        self.pca_result = pca.fit_transform(df_scaled)

        cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
        self.k = np.argmax(cumulative_variance >= self.variance_threshold) + 1
        print(
            f"Smallest k covering at least {self.variance_threshold * 100}% variance: {self.k}"
        )

    def save_results(self):
        """Saves PCA results with selected k components."""
        df_pca = pd.DataFrame(
            self.pca_result[:, : self.k],
            index=self.df.index,
            columns=[f"PC{i+1}" for i in range(self.k)],
        )
        df_pca.to_csv(self.output_path)

    def plot_cumulative_variance(self):
        """Plots cumulative explained variance and saves the figure."""
        pca = PCA()
        pca.fit(StandardScaler().fit_transform(self.df.iloc[:, :-1]))
        cumulative_variance = np.cumsum(pca.explained_variance_ratio_)

        plt.figure(figsize=(8, 6))
        plt.plot(
            range(1, len(cumulative_variance) + 1),
            cumulative_variance,
            marker="o",
            linestyle="--",
        )
        plt.axhline(
            y=self.variance_threshold,
            color="g",
            linestyle="--",
            label=f"{self.variance_threshold * 100}% Variance",
        )
        plt.axvline(x=self.k, color="r", linestyle="--", label=f"k={self.k}")
        plt.xlabel("Number of Principal Components")
        plt.ylabel("Cumulative Explained Variance")
        plt.title("PCA Cumulative Variance Plot")
        plt.legend()
        plt.savefig(self.output_plot)
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform PCA on gene expression data.")
    parser.add_argument(
        "--input_path", type=str, required=True, help="Path to input CSV file."
    )
    parser.add_argument(
        "--variance_threshold",
        type=float,
        default=0.8,
        help="Variance threshold for selecting PCs.",
    )
    parser.add_argument(
        "--output_path", type=str, required=True, help="Path to save PCA results."
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to save cumulative variance plot.",
    )

    args = parser.parse_args()

    pca_processor = PCAProcessor(
        args.input_path, args.variance_threshold, args.output_path, args.output_plot
    )
    pca_processor.load_data()
    pca_processor.apply_pca()
    pca_processor.save_results()
    pca_processor.plot_cumulative_variance()

    print(f"PCA results saved to {args.output_path}")
    print(f"Cumulative variance plot saved to {args.output_plot}")
