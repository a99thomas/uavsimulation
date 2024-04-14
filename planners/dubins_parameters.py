# dubins_parameters
#   - Dubins parameters that define path between two configurations
#
# mavsim_matlab 
#     - Beard & McLain, PUP, 2012
#     - Update history:  
#         3/26/2019 - RWB
#         4/2/2020 - RWB
#         3/30/2022 - RWB

import numpy as np


class DubinsParameters:
    '''
    Class that contains parameters for a Dubin's car path

    Attributes
    ----------
        p_s : np.ndarray (3x1)
            inertial position of start position, in meters
        chi_s : float
            course of start position in radians, measured from North
        p_e : np.ndarray (3x1)
            inertial position of end position, in meters
        chi_e : float
            course of end position in radians, measured from North
        R : float
            radius of start and end circles, from north
        center_s : np.ndarray (3x1)
            inertial center of start circle
        dir_s : int 
            direction of start circle: +1 CW, -1 CCW
        center_e : np.ndarray (3x1)
            inertial center of end circle
        dir_e : int 
            direction of end circle: +1 CW, -1 CCW
        length : float
            length of straight line segment
        r1 : np.ndarray (3x1)
            position on half plane for transition from start circle to straight-line
        n1 : np.ndarray (3x1)
            unit vector defining half plane for transition from start circle to straight-line, and from straight line to end circle
        r2 : np.ndarray (3x1)
            position on half plane for transition from straight line to end circle
        r3 : np.ndarray (3x1)
            position on half plane for end of dubins path
        n3 : np.ndarray (3x1)
            unit vector defining half plane for end of dubins path

    Methods
    ----------
    update(ps, chis, pe, chie, R)
        : create new Dubins path from start to end poses, with specified radius
    compute_parameters()
        : construct four dubins paths and pick the shortest and define all associated parameters.
    compute_points()
        : find equally spaced points along dubins path - for plotting and collision checking
    '''
    def __init__(self):
        self.radius = 150
        self.center_s = np.array([[0],[0],[-100]])
        self.center_e = np.zeros((3,1))
        self.q1 = np.zeros((3,1))
        self.dir_s = 0
        self.dir_m = 0
        self.dir_e = 0
        self.r1 = np.zeros((3,1))
        self.n1 = np.zeros((3,1))
        self.r2 = np.zeros((3,1))
        self.r3 = np.zeros((3,1))
        self.n3 = np.zeros((3,1))
        self.RLR = 0
        self.LRL = 0
        self.CCC = 0

    def update(self, ps, chis, pe, chie, R):
        self.p_s = ps
        self.chi_s = chis
        self.p_e = pe
        self.chi_e = chie
        self.radius = R
        self.compute_parameters()

    def compute_parameters(self):
        ps = self.p_s.reshape(3,)
        pe = self.p_e.reshape(3,)
        chis = self.chi_s
        chie = self.chi_e
        R = self.radius
        ell = np.linalg.norm(pe[0:2] - ps[0:2])

        ##### TODO #####
        if ell < 2 * R:
            print('Error in Dubins Parameters: The distance between nodes must be larger than 2R.')
            
        else:
            # print("chie",chie)
            # compute start and end circles
            crs = ps + R * rotz(np.pi/2) @ np.array([np.cos(chis), np.sin(chis), 0]).T
            cls = ps + R * rotz(-np.pi/2) @ np.array([np.cos(chis), np.sin(chis), 0]).T
            cre = pe + R * rotz(np.pi/2) @ np.array([np.cos(chie), np.sin(chie), 0]).T
            cle = pe + R * rotz(-np.pi/2) @ np.array([np.cos(chie), np.sin(chie), 0]).T
            
            if np.linalg.norm(pe-ps) <= 4*self.radius:
                if np.linalg.norm(cre-crs) <= 4*self.radius:
                    self.RLR = 1
                if np.linalg.norm(cle-cls) <= 4*self.radius:
                    self.LRL = 1
            # compute L1
            theta = np.arctan2(cre.item(1) - crs.item(1), cre.item(0) - crs.item(0))
            
            L1 = np.linalg.norm(crs-cre) + R * mod(2*np.pi + mod(theta - np.pi/2) - mod(chis - np.pi/2)) + R*mod(2*np.pi + mod(chie - np.pi/2) - mod(theta - np.pi/2))

            # compute L2
            theta2 = theta - np.pi/2 + np.arcsin(2*R/ell)

            d_theta_s = mod(2 * np.pi + mod(theta2) - mod(chis - np.pi/2))
            d_theta_e = mod(2 * np.pi + mod(theta2 + np.pi) - mod(chie + np.pi/2))
            L2 = np.sqrt(ell**2 - 4 * R **2) + R * d_theta_s + R * d_theta_e

            # compute L3
            theta3 = np.arccos(2*R/ell)
            d_theta_s = mod(2 * np.pi + mod(chis + np.pi/2) - mod(theta + theta3))
            d_theta_e = mod(2 * np.pi + mod(chie - np.pi/2) - mod(theta + theta3 - np.pi))
            L3 = np.sqrt(ell**2 - 4 * R **2) + R * d_theta_s + R * d_theta_e
            L4 = np.linalg.norm(cls-cle) + R * mod(np.pi * 2 + mod(chis + np.pi/2) - mod(theta + np.pi/2)) + R * mod(2 * np.pi + mod(theta + np.pi/2) - mod(chie + np.pi/2))
            L5 = np.inf
            L6 = np.inf
            if self.RLR:
                ell = np.linalg.norm(cre - crs)
                theta1 = mod(np.arccos(ell/(4*R)))
                cm = crs + 2*R*rotz(theta1) * (cre - crs)/ell
                L5 = R * mod(2*np.pi + mod(np.pi/2 - theta2) - mod(chis - np.pi/2)) + \
                R * mod(2*np.pi + mod(np.pi/2 - theta3) + mod(3*np.pi/2 - theta2)) + \
                R * mod(2*np.pi + mod(chie - np.pi/2) - mod(3*np.pi/2 - theta3))

            if self.LRL:
                ell = np.linalg.norm(cle - cls)
                theta1 = mod(np.arccos(ell/(4*R)))
                cm = crs + 2*R*rotz(theta1) * (cre - crs)/ell
                L6 = R * mod(2*np.pi - mod(np.pi/2 - theta2) + mod(chis - np.pi/2)) + \
                R * mod(2*np.pi + mod(np.pi/2 - theta3) - mod(3*np.pi/2 - theta2)) + \
                R * mod(2*np.pi - mod(chie - np.pi/2) + mod(3*np.pi/2 - theta3))


            # L is the minimum distance
            L = np.min([L1, L2, L3, L4, L5, L6])
            # print(L1, L2, L3, L4, L5, L6)
            min_idx = np.argmin([L1, L2, L3, L4, L5, L6])
            e1 = np.array([1, 0, 0])
            z1 = np.zeros((3,1))
            z2 = np.zeros((3,1))
            q1 = np.zeros((3,1))

            if min_idx == 0:
                self.center_s = crs
                self.dir_s = 1
                self.dir_e = 1
                self.center_e = cre
                q1 = (self.center_e - self.center_s) / np.linalg.norm(self.center_e - self.center_s)
                z1 = self.center_s + R * rotz(-np.pi/2) @ q1
                z2 = self.center_e + R * rotz(-np.pi/2) @ q1
                # print(self.center_e)
                # print(self.center_s)
                print("PATH 1!")
                pass

            elif min_idx == 1:
                self.center_s = crs
                self.dir_s = 1
                self.dir_e = -1
                self.center_e = cle
                ell = np.linalg.norm(self.center_e - self.center_s)
                theta = np.arctan2(self.center_e.item(1) - self.center_s.item(1), self.center_e.item(0) - self.center_s.item(0))
                theta2 = theta - np.pi/2 + np.arcsin(2 * R / ell)
                q1 = rotz(theta2 + np.pi/2) @ e1
                z1 = self.center_s + R * rotz(theta2) @ e1
                z2 = self.center_e + R * rotz(theta2 + np.pi) @ e1
                print("PATH 2!")
                pass

            elif min_idx == 2:
                self.center_s = cls
                self.dir_s = -1
                self.dir_e = 1
                self.center_e = cre
                ell = np.linalg.norm(self.center_e - self.center_s)
                theta = np.arctan2(self.center_e.item(1) - self.center_s.item(1), self.center_e.item(0) - self.center_s.item(0))
                theta2 = np.arccos(2 * R / ell)
                q1 = rotz(theta + theta2 - np.pi/2) @ e1
                z1 = self.center_s + R * rotz(theta + theta2) @ e1
                z2 = self.center_e + R * rotz(theta + theta2 - np.pi) @ e1
                print("PATH 3!")
                # quit()
                pass

            elif min_idx == 3:
                self.center_s = cls
                self.dir_s = -1
                self.dir_e = -1
                self.center_e = cle
                q1 = (self.center_e - self.center_s) / np.linalg.norm(self.center_e - self.center_s)
                z1 = self.center_s + R * rotz(np.pi/2) @ q1
                z2 = self.center_e + R * rotz(np.pi/2) @ q1
                print("PATH 4!")
                pass

            # elif min_idx == 4:
            #     self.CCC = 1
            #     self.center_s = cls
            #     self.dir_s = -1
            #     self.dir_m = 1
            #     self.dir_e = -1
            #     self.center_e = cle
            #     q1 = (self.center_e - self.center_s) / np.linalg.norm(self.center_e - self.center_s)
            #     z1 = self.center_s + R * rotz(-np.pi/2) @ q1
            #     # zm = self.center_m + R * rotz(np.pi/2) @ q1
            #     z2 = self.center_e + R * rotz(-np.pi/2) @ q1
            #     # print("PATH 5!")
            #     pass

            # elif min_idx == 5:
            #     self.CCC = 1
            #     self.center_s = cls
            #     self.dir_s = -1
            #     self.dir_m = 1
            #     self.dir_e = -1
            #     self.center_e = cle
            #     q1 = (self.center_e - self.center_s) / np.linalg.norm(self.center_e - self.center_s)
            #     z1 = self.center_s + R * rotz(np.pi/2) @ q1
            #     z2 = self.center_e + R * rotz(np.pi/2) @ q1
            #     # print("PATH 6!")
            #     pass
            # quit()
            self.length = L
            self.q1 = q1.reshape(3,1)
            # print(q1)
            # print(self.q1)
            # self.center_s = 0
            # self.dir_s = 0
            # self.center_e = 0
            # self.dir_e = 0
            self.r1 = z1.reshape(3,1)
            self.n1 = q1.reshape(3,1)
            self.r2 = z2.reshape(3,1)
            self.r3 = pe.reshape(3,1)
            self.n3 = (rotz(chie) @ e1).reshape(3,1)

    def compute_points(self):
        ##### TODO ##### - uncomment lines and remove last line
        Del = 0.1  # distance between point

        # points along start circle
        
        th1 =  np.arctan2(self.p_s.item(1) - self.center_s.item(1),
                         self.p_s.item(0) - self.center_s.item(0))
        th1 = mod(th1)
        th2 = np.arctan2(self.r1.item(1) - self.center_s.item(1),
                         self.r1.item(0) - self.center_s.item(0))
        th2 = mod(th2)
        th = th1
        theta_list = [th]
        if self.dir_s > 0:
            if th1 >= th2:
                while th < th2 + 2*np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2*np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)

        points = np.array([[self.center_s.item(0) + self.radius * np.cos(theta_list[0]),
                            self.center_s.item(1) + self.radius * np.sin(theta_list[0]),
                            self.center_s.item(2)]])
        for angle in theta_list:
            new_point = np.array([[self.center_s.item(0) + self.radius * np.cos(angle),
                                   self.center_s.item(1) + self.radius * np.sin(angle),
                                   self.center_s.item(2)]])
            points = np.concatenate((points, new_point), axis=0)

        # points along straight line
        sig = 0
        while sig <= 1:
            new_point = np.array([[(1 - sig) * self.r1.item(0) + sig * self.r2.item(0),
                                   (1 - sig) * self.r1.item(1) + sig * self.r2.item(1),
                                   (1 - sig) * self.r1.item(2) + sig * self.r2.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
            sig += Del

        # points along end circle
        th2 = np.arctan2(self.p_e.item(1) - self.center_e.item(1),
                         self.p_e.item(0) - self.center_e.item(0))
        th2 = mod(th2)
        th1 = np.arctan2(self.r2.item(1) - self.center_e.item(1),
                         self.r2.item(0) - self.center_e.item(0))
        th1 = mod(th1)
        th = th1
        theta_list = [th]
        if self.dir_e > 0:
            if th1 >= th2:
                while th < th2 + 2 * np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2 * np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)
        for angle in theta_list:
            new_point = np.array([[self.center_e.item(0) + self.radius * np.cos(angle),
                                   self.center_e.item(1) + self.radius * np.sin(angle),
                                   self.center_e.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
        
        return points


def rotz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])


def mod(x):
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


