import numpy as np

z_edges = np.linspace(0, 6, 60)
z_mid = np.array([(z_edges[i] + z_edges[i+1])/2 for i in range(len(z_edges)-1)])