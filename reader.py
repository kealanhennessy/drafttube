import pandas as pd
import numpy as np

def n_pt_x(col):
    lst = col.strip("Point (").split(")")
    lst.append(lst[1].strip(" = {"))
    lst.remove(lst[1])
    n = int(lst[0])
    x = float(lst[1])
    return n, x

def pt_y(col):
    return float(col)

def pt_z(col):
    return float(col)

def pt_size(col):
    col.strip(" ")
    col = col[:-2]
    return float(col)

def collate_all(col_n, g):
    return n_pt_x(g[col_n][0])[0], n_pt_x(g[col_n][0])[1], pt_y(g[col_n][1]), pt_z(g[col_n][2]), pt_size(g[col_n][3])

def collate_3d(col_n, g):
    return n_pt_x(g[col_n][0])[1], pt_y(g[col_n][1]), pt_z(g[col_n][2])
        
def collate_2d(col_n, g):
    return n_pt_x(g[col_n][0])[1], pt_y(g[col_n][1])