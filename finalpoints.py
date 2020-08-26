import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
import scipy.integrate as integrate
from planeInterp import *
import scipy
import numpy.fft
import scipy.fft
from copy import copy, deepcopy

class FinalPoints():
    def __init__(self):
        self.cone1, self.elbow1, self.connector1, self.middle1, self.diffuser1 = self.line_interp_set("1")
        self.cone4, self.elbow4, self.connector4, self.middle4, self.diffuser4 = self.line_interp_set("4")
        self.cone5, self.elbow5, self.connector5, self.middle5, self.diffuser5 = self.line_interp_set("5")
        self.cone6, self.elbow6, self.connector6, self.middle6, self.diffuser6 = self.line_interp_set("6")


        self.numberPointsTube = 100
        x0, y0, z0 = self.line_interp(self.cone1, self.elbow1, self.connector1, self.middle1, self.diffuser1, 0, num = self.numberPointsTube)
        x179, y179, z179= self.line_interp(self.cone1, self.elbow1, self.connector1, self.middle1, self.diffuser1, 179, num = self.numberPointsTube)
        x04, y04, z04 = self.line_interp(self.cone4, self.elbow4, self.connector4, self.middle4, self.diffuser4, 0, num = self.numberPointsTube)
        x1794, y1794, z1794= self.line_interp(self.cone4, self.elbow4, self.connector4, self.middle4, self.diffuser4, 179, num = self.numberPointsTube)
        x05, y05, z05= self.line_interp(self.cone5, self.elbow5, self.connector5, self.middle5, self.diffuser5, 0, num = self.numberPointsTube)
        x1795, y1795, z1795= self.line_interp(self.cone5, self.elbow5, self.connector5, self.middle5, self.diffuser5, 179, num = self.numberPointsTube)
        x06, y06, z06= self.line_interp(self.cone6, self.elbow6, self.connector6, self.middle6, self.diffuser6, 0, num = self.numberPointsTube)

        self.xorigin = [(2*x1794[i]+x06[i]+ x05[i])/4 for i in range(0, 55)]
        self.zorigin = [(2*z1794[i]+z06[i] +z05[i])/4 for i in range(0, 55)]
        for i in range(17):
            self.xorigin[i] = -0.1

        self.xorigin.extend([x1794[i] - 0.12 for i in range(63, 79)])
        self.zorigin.extend([z1794[i] + 0.025 for i in range(63, 79)])


        final = 3.4639191327609593
        distx = self.xorigin[-1] - self.xorigin[-2]
        distz = self.zorigin[-1] - self.zorigin[-2]

        while self.xorigin[len(self.xorigin)-1]  < final:
            self.xorigin.append(self.xorigin[-1]+distx)
            self.zorigin.append(self.zorigin[-1]+distz)
        self.xorigin[-1] = final
        self.zorigin[-1] = self.zorigin[-2]
        self.r_theta = []

    def points(self, tubeNum, spacing):
        """ take the lines of points and returns the x, y, z coordinates of the points and angles of the planes
        spacing is either equal, e, or chebyshev, c. These xyz need to be smoothed, then converted back to r, theta"""
        tck, u = interpolate.splprep((self.xorigin, self.zorigin), s=0)
        self.numPlanes = 20;
        if spacing =="e":
            points = np.linspace(0, 1, self.numPlanes) #equal spacing
            #print(points)
        elif spacing == "c":
            points = [(1-math.cos(math.pi*(2*(k+1)-2)/(2*self.numPlanes-2)))/2 for k in range(self.numPlanes)] #chebyshev spacing
            #print(points)
        self.xPo, self.zPo= interpolate.splev(points, tck, der = 0)
        xDer, zDer= interpolate.splev(points, tck, der = 1)

        allLines = []
        for i in range(360):
            if tubeNum == 1:
                xf, yf, zf = self.line_interp(self.cone1, self.elbow1, self.connector1, self.middle1, self.diffuser1, i, num = 100)
            elif tubeNum == 4:
                xf, yf, zf = self.line_interp(self.cone4, self.elbow4, self.connector4, self.middle4, self.diffuser4, i, num = 100)
            elif tubeNum == 5:
                xf, yf, zf = self.line_interp(self.cone5, self.elbow5, self.connector5, self.middle5, self.diffuser5, i, num = 100)
            elif tubeNum == 6:
                xf, yf, zf = self.line_interp(self.cone6, self.elbow6, self.connector6, self.middle6, self.diffuser6, i, num = 100)
            linePoints = [[xf[j], yf[j], zf[j]] for j in range(len(xf))]
            allLines.append(linePoints)
        everyPoint = []
        r_theta = []
        angles = []
        for w in range(self.numPlanes):
            newPlane = []
            planeNum = w
            point = [self.xPo[planeNum], 0, self.zPo[planeNum]] # location on origin line
            d = np.dot(point, [xDer[planeNum], 0, zDer[planeNum]])
            p_2 = [xDer[planeNum], 0, zDer[planeNum], d]
            a = np.sqrt((xDer[planeNum])**2 + (zDer[planeNum])**2)
            normal_vector = [xDer[planeNum]/a, 0, zDer[planeNum]/a]
            for j in range(360):
                line = allLines[j]
                point1, point2 = self.closest2index(line, normal_vector, point, planeNum/self.numPlanes)
                if np.sign(point1) == -1:
                    newPlane.append([line[99][0], line[99][1], line[99][2]])
                    continue
                elif np.sign(point2) == -1:
                    newPlane.append([line[0][0], line[0][1], line[0][2]])
                    continue
                newPoints = np.linspace(point2/100, point1/100, 15)
                if tubeNum == 1:
                    x2, y2, z2 = self.line_interp(self.cone1, self.elbow1, self.connector1, self.middle1, self.diffuser1, j, unew = newPoints)
                elif tubeNum == 4:
                    x2, y2, z2 = self.line_interp(self.cone4, self.elbow4, self.connector4, self.middle4, self.diffuser4, j, unew = newPoints)
                elif tubeNum == 5:
                    x2, y2, z2 = self.line_interp(self.cone5, self.elbow5, self.connector5, self.middle5, self.diffuser5, j, unew = newPoints)
                elif tubeNum == 6:
                    x2, y2, z2 = self.line_interp(self.cone6, self.elbow6, self.connector6, self.middle6, self.diffuser6, j, unew = newPoints)
                point1, point2 = self.closest2index2([x2, y2, z2], normal_vector, point)
                line = [x2[point1], y2[point1], z2[point1]] +[x2[point2]-x2[point1], y2[point2]- y2[point1],z2[point2]- z2[point1]]
                t = (p_2[3] - (p_2[0]*line[0] +p_2[1]*line[1] +p_2[2]*line[2]))/(p_2[0]*line[3] +p_2[1]*line[4] + p_2[2]*line[5])
                finalPoint = [(line[0]+t*line[3]), (line[1]+t*line[4]), (line[2]+t*line[5])]
                newPlane.append([finalPoint[0], finalPoint[1], finalPoint[2]])
            for j in range(len(newPlane)):
                newPlane[j][0] = newPlane[j][0] - point[0]
                newPlane[j][2] = newPlane[j][2] - point[2]
            if np.sign(xDer[planeNum]*zDer[planeNum]) == -1:
                phi =  np.arctan(np.abs(xDer[planeNum]/zDer[planeNum])) # angle of the plane
            elif np.sign(xDer[planeNum]*zDer[planeNum]) == 1:
                phi = math.pi/2 + np.arctan(zDer[planeNum]/xDer[planeNum])
            elif np.sign(xDer[planeNum]) == 0:
                phi = 0;
            else:
                phi = math.pi/2
            y_axis = [0,1,0]
            for i in range(360):
                if phi == math.pi/2:
                    newPlane[i] = [newPlane[i][2], newPlane[i][1], newPlane[i][0]]
                if phi > math.pi/2:
                    newPlane[i] = [newPlane[i][2], newPlane[i][1], -newPlane[i][0]]
                pt = newPlane[i]
                dot_product = np.dot(y_axis, pt)
                y = np.array(y_axis)
                p = np.array(pt)
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
                newPlane[i] = [r, theta, 0]
            t = [i for i in range(360)]
            t_rad = [np.deg2rad(t[i]) for i in range(len(t))]
            print(w)
            f = self.interp(newPlane)
            r_values = f(t_rad)
            for i in range(len(t)):
                point_temp = [r_values[i]*np.cos(t_rad[i]), r_values[i]*np.sin(t_rad[i]), 0]
                x, y, z = self.plane_rotation(phi,point_temp)
                xf = x + point[0]
                zf = z + point[2]
                if tubeNum == 1 and zf<-1.14297783e+00:
                    zf = -1.14297783e+00
                newPlane[i] = [xf, y, zf] # 360 points
            everyPoint.append(newPlane)
            angles.append(phi)
            r_theta.append(r_values)
        r_theta = self.smoothing(r_theta)
        return r_theta, everyPoint, angles

    def pol2cart(self, points, angles):
        point_cart = []
        for w in range(self.numPlanes):
            r_values = points[w]
            newPlane = []
            t = [i for i in range(360)]
            t_rad = [np.deg2rad(t[i]) for i in range(len(t))]
            point = [self.xPo[w], 0, self.zPo[w]]
            phi = angles[w]
            for i in range(len(t)):
                point_temp = [r_values[i]*np.cos(t_rad[i]), r_values[i]*np.sin(t_rad[i]), 0]
                x, y, z = self.plane_rotation(phi,point_temp)
                xf = x + point[0]
                zf = z + point[2]
                newPlane.append([xf, y, zf]) # 360 points
            point_cart.append(newPlane)
        return point_cart


    def line_interp_set(self, num):
        """sets up the baseline interpolator"""
        baseline1 = planeInterpolator(num)
        cone = baseline1.cone_interp()
        elbow = baseline1.elbow_interp()
        connector = baseline1.cmd_interp("connector")
        middle = baseline1.cmd_interp("middle")
        diffuser = baseline1.cmd_interp("diffuser")
        return cone, elbow, connector, middle, diffuser

    def line_interp(self, cone, elbow, connector, middle, diffuser, k, num = 50, unew = []):
        """gives the kth line 0<=k<=359 num is a string and k is an interger"""
        """num if the number of points you want for the draft tube lines"""
        """the u values can also be defined specifically with unew"""
        x = []
        y = []
        z = []
        dist = 0.001
        i = 0
        for j in range(len(cone)):
            if k<180:
                if i == 0:
                    x.append(cone[i][k+180][0])
                    y.append(cone[i][k+180][1])
                    z.append(cone[i][k+180][2]+0.10011979192495346)
                    i = i+1
                else:
                    xnext = cone[j][k+180][0]
                    ynext = cone[j][k+180][1]
                    znext = cone[j][k+180][2]+0.10011979192495346
                    distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                    #print(distance)
                    if distance > dist:
                        x.append(xnext)
                        y.append(ynext)
                        z.append(znext)
                        i = i+1
            else:
                if i == 0:
                    x.append(cone[i][k-180][0])
                    y.append(cone[i][k-180][1])
                    z.append(cone[i][k-180][2]+0.10011979192495346)
                    i = i+1
                else:
                    xnext = cone[j][k-180][0]
                    ynext = cone[j][k-180][1]
                    znext = cone[j][k-180][2]+0.10011979192495346
                    distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                    if distance > dist:
                        x.append(xnext)
                        y.append(ynext)
                        z.append(znext)
                        i = i+1
        a = i
        for j in list(range(0, 31)) + [63, 64, 32, 31] + list(range(33, 63)):
            if i-a == 0:
                x.append(elbow[i-a][k][0])
                y.append(elbow[i-a][k][1])
                z.append(elbow[i-a][k][2]+0.10011979192495346)
                i = i+1
            else:
                xnext = elbow[j][k][0]
                ynext = elbow[j][k][1]
                znext = elbow[j][k][2]+0.10011979192495346
                distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                if distance > dist:
                    x.append(xnext)
                    y.append(ynext)
                    z.append(znext)
                    i = i+1
        b = i
        for j in range(65):
            if i-b == 0:
                x.append(connector[i-b][k][0])
                y.append(connector[i-b][k][1])
                z.append(connector[i-b][k][2]+0.10011979192495346)
                i = i+1
            else:
                xnext = connector[j][k][0]
                ynext = connector[j][k][1]
                znext = connector[j][k][2]+0.10011979192495346
                distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                if distance > 0.001:
                    x.append(xnext)
                    y.append(ynext)
                    z.append(znext)
                    i = i+1
        b = i
        for j in range(65):
            if i-b == 0:
                x.append(middle[i-b][k][0])
                y.append(middle[i-b][k][1])
                z.append(middle[i-b][k][2]+0.10011979192495346)
                i = i+1
            else:
                xnext = middle[j][k][0]
                ynext = middle[j][k][1]
                znext = middle[j][k][2]+0.10011979192495346
                distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                if distance > 0.001:
                    x.append(xnext)
                    y.append(ynext)
                    z.append(znext)
                    i = i+1
        b = i
        for j in range(65):
            if i-b == 0:
                x.append(diffuser[i-b][k][0])
                y.append(diffuser[i-b][k][1])
                z.append(diffuser[i-b][k][2]+0.10011979192495346)
                i = i+1
            else:
                xnext = diffuser[j][k][0]
                ynext = diffuser[j][k][1]
                znext = diffuser[j][k][2]+0.10011979192495346
                distance = (((xnext-x[i-1])**2) + ((ynext-y[i-1])**2) + ((znext-z[i-1])**2))**0.5
                if distance > 0.001:
                    x.append(xnext)
                    y.append(ynext)
                    z.append(znext)
                    i = i+1
        if len(unew) == 0:
            unew = np.linspace(0, 1, num, endpoint=False)
        tck, u = interpolate.splprep([x, y,z], s=0)
        xf, yf, zf = interpolate.splev(unew, tck, der = 0)
        return xf, yf, zf


    def closest2index(self, linePoints, normal_vector, point, a):
        """points is a list of [x, y, z] points and plane is by [a, b, c, d], where ax+by+cz+d = 0"""
        # a is how far along the line the plane is from 0 to 1
        dist = []
        index = math.floor(a*len(linePoints)) # gives the index along the line that is near the point
        b = math.floor(0.22*len(linePoints)) # 0.1 is arbitrary, need to choose a number that will limit search
        c = index - b
        if c < 0 :
            c = 0
        d = index+b
        if d>len(linePoints):
            d = len(linePoints)
        testPoints = np.array(np.linspace(c, d , d - c, endpoint=False))
        for i in testPoints:
            i = int(i)
            vector = [linePoints[i][0]-point[0], linePoints[i][1]-point[1], linePoints[i][2]-point[2]]
            normal_distance = np.dot(vector, normal_vector)
            dist.append([i, normal_distance])
            distance = (((point[0]-linePoints[i][0])**2) + ((point[1]-linePoints[i][1])**2) + ((point[2]-linePoints[i][2])**2))**0.5
            if distance > 1 and point[0]< 1:
                dist[len(dist)-1][1] += 15
            elif distance > 2:
                dist[len(dist)-1][1] += 15
        closest = min(range(len(dist)), key = lambda i: (np.abs(dist[i][1])))
        if dist[closest][0] == 0 and np.sign(dist[closest][1]) == np.sign(dist[closest+1][1]):
            return 0, -1
        elif dist[closest][0] == 99 and np.sign(dist[closest][1]) == np.sign(dist[closest-1][1]):  # change later
            return -1, 0
        dist1 = dist[closest][1]
        secondvalue = 1000
        for i in range(len(dist)):
            if np.abs(dist[i][1]) < np.abs(secondvalue) and np.abs(dist[i][1]) > np.abs(dist1) and not np.sign(dist1) == np.sign(dist[i][1]):
                second = i
                secondvalue = dist[i][1]
        return dist[closest][0], dist[second][0]

    def closest2index2(self, linePoints, normal_vector, point):
        """points is a list of [x, y, z] points and plane is by [a, b, c, d], where ax+by+cz+d = 0"""
        dist = []
        for i in range(len(linePoints[0])):
            vector = [linePoints[0][i]-point[0], linePoints[1][i]-point[1], linePoints[2][i]-point[2]]
            normal_distance = np.dot(vector, normal_vector)
            dist.append([i, normal_distance])
        closest = min(range(len(dist)), key = lambda i: (np.abs(dist[i][1])))
        dist1 = dist[closest][1]
        secondvalue = 1000
        for i in range(len(dist)):
            if np.abs(dist[i][1]) < np.abs(secondvalue) and np.abs(dist[i][1]) > np.abs(dist1) and not np.sign(dist1) == np.sign(dist[i][1]):
                second = i
                secondvalue = dist[i][1]
        return dist[closest][0], dist[second][0]

    def interp(self, points):
        """takes in the unorganized points at (0,0,0) and in r, theta and returns the interpolated functions"""
        r_theta = np.array([[points[i][0], points[i][1]] for i in range(len(points))]) # get r and theta from points
        sort = np.argsort(r_theta[:,1])
        r_theta_sorted = r_theta[sort,:] # sorted pairs
        r_plot = [r_theta_sorted[i][0] for i in range(len(r_theta_sorted))] # get r values form sorted
        theta_plot = [r_theta_sorted[i][1] for i in range(len(r_theta_sorted))] # get theta values from sorted
        for i in range(359):
            if theta_plot[i] == theta_plot[i+1]:
                print(theta_plot[i])
        f = interpolate.PchipInterpolator(theta_plot, r_plot)
        return f

    def plane_rotation(self,phi, point):
        """rotate in 3d space takes in a point in cartesian coordinates 3d space and rotates it"""
        x = point[0]*np.cos(-phi) + point[2]*np.sin(-phi)
        y = point[1]
        z = -point[0]*np.sin(-phi) + point[2]*np.cos(-phi)
        return x, y, z

    def smoothing(self, r):
        # r is a 2D matrix that contain the radius r(z,theta) values in the matrix in this format:
        # The columns represent different theta values at any given row
        # The rows represent z along the origin line
        # I.E. Each row is at a specific z on the origin line and the columns on that row represent r values at thetas from 0-2*pi

        # Transforms
        r_f1=scipy.fft.dct(r,axis=0)
        r_f=numpy.fft.rfft(r_f1,axis=1)

        # Smoothing
        cutoff_z=5
        cutoff_theta=120
        r_f_smooth = deepcopy(r_f)

        for k_z in range(r_f.shape[0]):
            for k_theta in range(r_f.shape[1]):
                if k_theta>cutoff_theta and k_z>cutoff_z:
                    dist=-1*(k_theta**2+k_z**2)**0.5
                    r_f_smooth[k_z,k_theta]=np.exp(dist)*r_f_smooth[k_z,k_theta]

        r_if1=numpy.fft.irfft(r_f_smooth,axis=1)
        r_smooth=scipy.fft.idct(r_if1,axis=0)

        # r is the original matrix. r_smooth is the smoothed out matrix

        return r_smooth
