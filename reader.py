import pandas as pd
import numpy as np

class Reader():
    def __init__(self, path, cascade=True):
        if cascade:
            df = pd.read_csv(path, skiprows=1)
        else:
            df = pd.read_csv(path)
        self.geo = np.vstack([list(df.columns), df.to_numpy()])

    def n_pt_x(self, col):
        lst = col.strip("Point (").split(")")
        lst.append(lst[1].strip(" = {"))
        lst.remove(lst[1])
        n = int(lst[0])
        x = float(lst[1])
        return n, x

    def pt_y(self, col):
        return float(col)

    def pt_z(self, col):
        return float(col)

    def pt_size(self, col):
        col.strip(" ")
        col = col[:-2]
        return float(col)

    def collate_all(self, col_n, g=None):
        if g is None:
            g = self.geo
        return self.n_pt_x(g[col_n][0])[0], self.n_pt_x(g[col_n][0])[1], self.pt_y(g[col_n][1]), self.pt_z(g[col_n][2]), self.pt_size(g[col_n][3])

    def collate_pt(self, col_n, g=None):
        if g is None:
            g = self.geo
        return self.n_pt_x(g[col_n][0])[0], self.n_pt_x(g[col_n][0])[1], self.pt_y(g[col_n][1]), self.pt_z(g[col_n][2])
    
    def collate_3d(self, col_n, g=None):
        if g is None:
            g = self.geo
        return self.n_pt_x(g[col_n][0])[1], self.pt_y(g[col_n][1]), self.pt_z(g[col_n][2])
        
    def collate_2d(self, col_n, g=None):
        if g is None:
            g = self.geo
        return self.n_pt_x(g[col_n][0])[1], self.pt_y(g[col_n][1])
    
    def size(self, g=None):
        if g is None:
            g = self.geo
        return len(g)