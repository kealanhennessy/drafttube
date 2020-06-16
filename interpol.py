import numpy as np
import sys
sys.path.append("/Users/kealan/Desktop/drafttube/")
from reader import *
import matplotlib.pyplot as plt
import math
from scipy import interpolate

class Interpolate():
    def __init__(self, points, dir):
        self.points = points # Specify the number of points we expect to represent a single planar figure (S) in our .geo file
        self.dir = dir # Specify the directory of our .geo file. Should be a single portion of the draft tube

    def cart2pol(self, x, y):
        """Convert a set of cartesian coordinates x, y to polar coordinates r, theta (default to radians)"""
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return r, theta

    def closest(self, lst, K):
        """Return the nearest value to integer K in a list of floating point values lst"""
        ind = min(range(len(lst)), key = lambda i: abs(lst[i] - K))
        return [lst[ind - 1], lst[ind]] # Has to be mutable

    def read_polar(self, x, y):
        """Take the our i = 10000 interpolated points in Cartesian (x, y) space and return each 
        (r, theta) vector from [1, 360] degrees (360 = 0). Note that these vectors are not ordered.
        FIX: The indexes are messed up - [90] is 270deg, [180] is 360deg, [270] is 90deg, [0] is 180deg
        """
        pol = np.asarray([cart2pol(y[i], x[i]) for i in range(len(x))]) # Convert to an (r, theta) numpy array, default in radians
        pol_f = np.ndarray.flatten(pol) # Flatten to be able to manipulate and still be able to re-pair later on
        for i in range(1, len(pol_f), 2): # Convert radians to degrees
            deg = np.rad2deg(pol_f[i])
            pol_f[i] = deg
        degs = [] # List to append re-pair
        for i in range(-180, 180): # Find the closest discreet deg value for 360deg out of 10,000 points
            degs.append(closest(pol_f, i))
        assert np.shape(degs)[0] is 360
        for i in range(len(degs)):
            r = int(round(degs[i][1])) # Convert the new degree values to integers
            degs[i][1] = r
            if degs[i][1] <= 0: # Convert -180:180 to 0:360
                degs[i][1] += 360
        return degs

    def read_points(self):
        """Reads the .geo file and returns the number of points as specified in a single list"""
        geo = Reader(self.dir, cascade=False)
        pts = []
        for i in range(0, self.points):
            pts.append(np.array(geo.collate_xy(i)))
        points = np.asarray(pts)
        return points
        
    def interpolate_spline(self, points):
        """Interpolate the single list of points given by read_points."""
        # Split the points into respective lists
        x = np.array([points[i][0] for i in range(len(points))])
        y = np.array([points[i][1] for i in range(len(points))])
        # Append the starting x,y coordinates
        x = np.r_[x, x[0]]
        y = np.r_[y, y[0]]
        # Fit splines to x = f(u) and y = g(u), treating both as periodic. Also note that the input parameter 
        # s = 0 is needed in order to force the spline fit to pass through all the input points.
        tck, u = interpolate.splprep([x, y], s=0, per=True)
        # Evaluate the spline fits for n = 10000 evenly spaced distance values
        # 10000 points should be good enough to be able to approximate 360 to a high degree of accurately
        xi, yi = interpolate.splev(np.linspace(0, 1, 10000), tck)
        return read_polar(xi, yi)