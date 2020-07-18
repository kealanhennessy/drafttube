import numpy as np
from reader import *
import matplotlib.pyplot as plt
import math
from scipy import interpolate
import sys

class planeInterpolator():
    def __init__(self, baselineNum):
        self.num = baselineNum
        
    def cone_interp(self):
        """Read the points from the cone and return all 360 points per layer"""
        geo = Reader("drafttube_shapes/baseline" + self.num +"/geo/cone.geo", cascade=False)
        allPoints = []
        t = [i for i in range(0, 360)]
        t_rad = [np.deg2rad(t[i]) for i in range(len(t))]
        for k in range(1, 129):
            pts = []
            points_final = []
            if k == 1:
                index = list(range(256))
            else:
                index = list(range((k-1)*128, (k-1)*128 + 128))
            for i in index:
                pts.append(np.array(geo.collate_xyz(i)))
            points = np.asarray(pts)
            points, x_mid, z_mid = self.shift2origin(points)
            for i in range(len(index)):
                r , theta = self.cart2pol(points[i][0], points[i][1])
                points[i] = [r, theta, 0]
            f = self.interp(points)
            r_values = f(t_rad)
            for i in range(len(t)):
                point = [r_values[i]*np.cos(t_rad[i]), r_values[i]*np.sin(t_rad[i]), 0]
                xf = point[0] + x_mid
                zf = point[2] + z_mid
                points_final.append([xf, point[1], zf])
            #points_final = np.array(points_final)
            allPoints.append(points_final)
        return allPoints

    def cmd_interp(self, shape):
        """takes in a string that is the drafttube_shape
        Read the points from the cone and return all 360 points per layer"""
        geo = Reader("drafttube_shapes/baseline" + self.num + "/geo/" + shape + ".geo", cascade=False)
        allPoints = []
        t = [i for i in range(0, 360)]
        t_rad = [np.deg2rad(t[i]) for i in range(len(t))]
        k = 130
        for k in range(65):
            pts = []
            points_final = []
            for i in range(k, 8320, 65):
                pts.append(np.array(geo.collate_xyz(i)))
            points = np.asarray(pts)
            points, x_mid, z_mid = self.shift2origin(points)
            ind = self.closest_index(points[:, 1], 0)
            phi = math.pi/2 # angle of the plane
            y_axis = [0,1,0]
            for i in range(len(index)):
                pt = points[i]
                dot_product = np.dot(y_axis, pt)
                y = np.array(y_axis)
                p = np.array(pt)

                # account for specific quadrant, is there a function that gives these values for arccos?
                theta_original = np.arccos(dot_product/(np.linalg.norm(p)*np.linalg.norm(y)))
                if p[0] > 0 and p[1] >= 0:
                    theta = (math.pi/2) - theta_original
                elif p[1] >= 0:
                    theta = (math.pi/2) + theta_original
                elif p[0] < 0 and p[1] < 0:
                    theta = (math.pi/2) + theta_original
                else:
                    theta = (5*math.pi/2) - theta_original
                    r = np.linalg.norm(p)
                    points[i] = [r, theta, 0]
            f = self.interp(points)
            r_values = f(t_rad)
            for i in range(len(t)):
                point = [r_values[i]*np.cos(t_rad[i]), r_values[i]*np.sin(t_rad[i]), 0]
                x, y, z = self.plane_rotation(phi,point)
                xf = x + x_mid
                zf = z + z_mid
                points_final.append([xf, y, zf])
            points_final = np.array(points_final)
            allPoints.append(points_final)
        return allPoints

    def elbow_interp(self):
        "Read points from elbow and returns all of the 360 points per layer"
        geo = Reader("drafttube_shapes/baseline"+ self.num + "/geo/elbow.geo", cascade=False)
        allPoints = []
        t = [i for i in range(0, 360)]
        t_rad = [np.deg2rad(t[i]) for i in range(len(t))]
        for k in range(1, 66):
            pts = []
            points_final = []
            if k == 1:
                index = list(range(256))
            elif k <= 30 or (k>=34 and k<=63):
                index = list(range((k-1)*128, (k-1)*128 + 128))
            elif k == 32:
                index = list(range(3969, 4051, 2))+list(range(4050, 4098, 2))+list(range(4099, 4147, 2))+list(range(4146, 4223, 2))
            elif k == 33:
                index = list(range(3968, 4050, 2))+list(range(4051, 4099, 2))+list(range(4098, 4146, 2))+list(range(4147, 4225, 2))
            elif k == 64:
                index = list(range(8064,8093)) + list(range(8126, 8175)) + list(range(8224, 8248))+ list(range(8272, 8296)) + [8106, 8121]
            else:
                index = list(range(8093,8106)) + list(range(8107, 8121)) + list(range(8122, 8126))+ list(range(8175, 8224))+list(range(8248,8272))+list(range(8296, 8320))
            for i in index:
                pts.append(np.array(geo.collate_xyz(i)))
            points = np.asarray(pts)
            points, x_mid, z_mid = self.shift2origin(points)

            ind = self.closest_index(points[:, 1], 0)
            phi = np.arctan(np.array(abs(points[ind][2]/points[ind][0]))) # angle of the plane
            y_axis = [0,1,0]
            for i in range(len(index)):
                pt = points[i]
                dot_product = np.dot(y_axis, pt)
                y = np.array(y_axis)
                p = np.array(pt)

                # account for specific quadrant, is there a function that gives these values for arccos?
                theta_original = np.arccos(dot_product/(np.linalg.norm(p)*np.linalg.norm(y)))
                if p[0] > 0 and p[1] >= 0:
                    theta = (math.pi/2) - theta_original
                elif p[1] >= 0:
                    theta = (math.pi/2) + theta_original
                elif p[0] < 0 and p[1] < 0:
                    theta = (math.pi/2) + theta_original
                else:
                    theta = (5*math.pi/2) - theta_original
                    r = np.linalg.norm(p)
                    points[i] = [r, theta, 0]
            f = self.interp(points)
            r_values = f(t_rad)
            for i in range(len(t)):
                point = [r_values[i]*np.cos(t_rad[i]), r_values[i]*np.sin(t_rad[i]), 0]
                x, y, z = self.plane_rotation(phi,point)
                xf = x + x_mid
                zf = z + z_mid
                points_final.append([xf, y, zf])
            points_final = np.array(points_final)
            allPoints.append(points_final)
        return allPoints

    def interp(self, points):
        """takes in the unorganized points at (0,0,0) and in r, theta and returns the interpolated functions"""
        r_theta = np.array([[points[i][0], points[i][1]] for i in range(len(points))]) # get r and theta from points
        sort = np.argsort(r_theta[:,1])
        r_theta_sorted = r_theta[sort,:] # sorted pairs
        r_plot = [r_theta_sorted[i][0] for i in range(len(r_theta_sorted))] # get r values form sorted
        theta_plot = [r_theta_sorted[i][1] for i in range(len(r_theta_sorted))] # get theta values from sorted
        f = interpolate.PchipInterpolator(theta_plot, r_plot, extrapolate = True)
        return f

    def shift2origin(self, points):
        """takes in the points and shifts them to the origin, return the points, and
         x_mid and z_mid"""
        x_mid = np.average(points[:,0])
        z_mid = np.average(points[:,2])
        for i in range(len(points)):
            points[i][0] = points[i][0] - x_mid
            points[i][2] = points[i][2] - z_mid
        return points, x_mid, z_mid

    def plane_rotation(self,phi, point):
        """rotate in 3d space takes in a point in cartesian coordinates 3d space and rotates it"""
        x = point[0]*np.cos(-phi) + point[2]*np.sin(-phi)
        y = point[1]
        z = -point[0]*np.sin(-phi) + point[2]*np.cos(-phi)
        return x, y, z

    def cart2pol(self,x, y):
        """Convert a set of cartesian coordinates x, y to polar coordinates r, theta (default to radians)"""
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return r, theta

    def closest_index(self,lst, K):
        """Return the nearest value to integer K in a list of floating point values lst"""
        ind = min(range(len(lst)), key = lambda i: abs(lst[i] - K))
        return ind