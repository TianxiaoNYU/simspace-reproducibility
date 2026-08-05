"""Import and minimal execution tests for the shared benchmark environment."""

from __future__ import annotations

import importlib

import numpy as np
import torch
from torch_geometric.data import Data
from torch_sparse import SparseTensor


EXPECTED_IMPORTS = {
    "GraphST": "GraphST",
    "STAGATE": "STAGATE_pyG",
    "SpaGCN": "SpaGCN",
    "spCLUE": "spCLUE",
}


for label, module_name in EXPECTED_IMPORTS.items():
    module = importlib.import_module(module_name)
    print(f"PASS import {label}: {module.__file__}")

# Exercise the compiled PyG extensions required by STAGATE.
edge_index = torch.tensor([[0, 1], [1, 0]], dtype=torch.long)
graph = Data(x=torch.ones((2, 3)), edge_index=edge_index)
sparse = SparseTensor.from_edge_index(graph.edge_index, sparse_sizes=(2, 2))
assert sparse.nnz() == 2
assert np.isfinite(graph.x.numpy()).all()
print("PASS PyTorch/PyG/torch_sparse minimal execution")

# Verify that rpy2 is bound to the external R library and can load mclust.
from rpy2.robjects.packages import importr

mclust = importr("mclust")
print(f"PASS host R mclust: {mclust.__version__}")
