"""Approximate kNN transformer (pynndescent) returning a CSR distance matrix in the same
format as sklearn's KNeighborsTransformer(mode='distance'). Works around a numba typing
error in scanpy's built-in 'pynndescent' path with the current umap-learn/numba versions."""
import numpy as np, scipy.sparse as sp
from pynndescent import NNDescent
from sklearn.base import BaseEstimator, TransformerMixin

class NNDTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, n_neighbors=20, metric="euclidean", random_state=0, n_jobs=-1):
        self.n_neighbors, self.metric, self.random_state, self.n_jobs = n_neighbors, metric, random_state, n_jobs
    def fit(self, X, y=None):
        self.index_ = NNDescent(X, n_neighbors=self.n_neighbors, metric=self.metric,
                                random_state=self.random_state, n_jobs=self.n_jobs)
        self.n_ = X.shape[0]; return self
    def transform(self, X):
        idx, d = self.index_.neighbor_graph
        idx, d = np.ascontiguousarray(idx), np.ascontiguousarray(d, dtype=np.float32)
        n, k = idx.shape
        return sp.csr_matrix((d.ravel(), idx.ravel(), np.arange(0, n * k + 1, k)), shape=(n, self.n_))
    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)
