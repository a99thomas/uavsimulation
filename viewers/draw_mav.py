"""
mavsim_python: drawing tools
Aaron Thomas 1/16/24
"""
import numpy as np
import pyqtgraph.opengl as gl
from tools.rotations import euler_to_rotation
from tools.drawing import rotate_points, translate_points, points_to_mesh


class DrawMav:
    def __init__(self, state, window, scale=10):
        """
        Draw the Spacecraft.

        The input to this function is a (message) class with properties that define the state.
        The following properties are assumed:
            state.north  # north position
            state.east  # east position
            state.altitude   # altitude
            state.phi  # roll angle
            state.theta  # pitch angle
            state.psi  # yaw angle
        """
        self.unit_length = scale
        mav_position = np.array([[state.north], [state.east], [-state.altitude]])  # NED coordinates
        # attitude of mav as a rotation matrix R from body to inertial
        R_bi = euler_to_rotation(state.phi, state.theta, state.psi)
        # convert North-East Down to East-North-Up for rendering
        self.R_ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        # get points that define the non-rotated, non-translated spacecraft and the mesh colors
        self.mav_points, self.mav_index, self.mav_meshColors = self.get_mav_points()
        self.mav_body = self.add_object(
            self.mav_points,
            self.mav_index,
            self.mav_meshColors,
            R_bi,
            mav_position)
        window.addItem(self.mav_body)  # add spacecraft to plot     

    def update(self, state):
        mav_position = np.array([[state.north], [state.east], [-state.altitude]])  # NED coordinates
        # attitude of mav as a rotation matrix R from body to inertial
        R_bi = euler_to_rotation(state.phi, state.theta, state.psi)
        self.mav_body = self.update_object(
            self.mav_body,
            self.mav_points,
            self.mav_index,
            self.mav_meshColors,
            R_bi,
            mav_position)

    def add_object(self, points, index, colors, R, position):
        rotated_points = rotate_points(points, R)
        translated_points = translate_points(rotated_points, position)
        translated_points = self.R_ned @ translated_points
        # convert points to triangular mesh defined as array of three 3D points (Nx3x3)
        mesh = points_to_mesh(translated_points, index)
        object = gl.GLMeshItem(
            vertexes=mesh,  # defines the triangular mesh (Nx3x3)
            vertexColors=colors,  # defines mesh colors (Nx1)
            drawEdges=True,  # draw edges between mesh elements
            smooth=False,  # speeds up rendering
            computeNormals=False)  # speeds up rendering
        return object

    def update_object(self, object, points, index, colors, R, position):
        rotated_points = rotate_points(points, R)
        translated_points = translate_points(rotated_points, position)
        translated_points = self.R_ned @ translated_points
        # convert points to triangular mesh defined as array of three 3D points (Nx3x3)
        mesh = points_to_mesh(translated_points, index)
        object.setMeshData(vertexes=mesh, vertexColors=colors)
        return object

    def get_mav_points(self, fuse_l1=2, fuse_l2=1, fuse_l3=4, fuse_w=1, fuse_h=1, wing_l=1.5, wing_w=4, tailwing_l=1, tail_h=1, tailwing_w=2.5):
        """"
            Points that define the spacecraft, and the colors of the triangular mesh
            Define the points on the spacecraft following information in Appendix C.3
        """
        # points are in XYZ coordinates
        #   define the points on the spacecraft according to Appendix C.3
        points = self.unit_length * np.array([
            [fuse_l1, 0, 0],  # point 1 [0]
            [fuse_l2, fuse_w/2, -fuse_h/2],  # point 2 [1]
            [fuse_l2, -fuse_w/2, -fuse_h/2],  # point 3 [2]
            [fuse_l2, -fuse_w/2, fuse_h/2],  # point 4 [3]
            [fuse_l2, fuse_w/2, fuse_h/2],  # point 5 [4]
            [-fuse_l3, 0, 0],  # point 6 [5]
            [0, wing_w/2, 0],  # point 7 [6]
            [-wing_l, wing_w/2, 0],  # point 8 [7]
            [-wing_l, -wing_w/2, 0],  # point 9 [8]
            [0, -wing_w/2, 0],  # point 10 [9]
            [-fuse_l3+tailwing_l, tailwing_w/2, 0],  # point 11 [10]
            [-fuse_l3, tailwing_w/2, 0],  # point 12 [11]
            [-fuse_l3, -tailwing_w/2, 0],
            [-fuse_l3+tailwing_l, -tailwing_w/2, 0],
            [-fuse_l3+tailwing_l, 0, 0],
            [-fuse_l3, 0, -tail_h]
            ]).T
        # point index that defines the mesh
        index = np.array([
            [0, 1, 2],  # front 1
            [0, 2, 3],  # front 2
            [0, 3, 4],  # front 3
            [0, 4, 1],  # front 4
            [1, 5, 2],  # top 1
            [4, 5, 3],  # bottom 1
            [2, 5, 3], # left 1
            [1, 5, 4], # right 1
            [6, 7, 8],  # wing 1
            [8, 9, 6],  # wing 2
            [10, 11, 12],  # tailwing 1
            [12, 13, 10],  # tailwing 2
            [5, 14, 15],  # tail 1
            [14, 15, 5],  # tail 2
            ])
        #   define the colors for each face of triangular mesh
        red = np.array([1., 0., 0., 1])
        green = np.array([0., 1., 0., 1])
        blue = np.array([0., 0., 1., 1])
        yellow = np.array([1., 1., 0., 1])
        purple = np.array([1., 0., 1., 1])
        meshColors = np.empty((13, 3, 4), dtype=np.float32)
        meshColors[0] = yellow  # front 1
        meshColors[1] = yellow  # front 2
        meshColors[2] = yellow  # front 3
        meshColors[3] = yellow  # front 4
        meshColors[4] = blue  # top 1
        meshColors[5] = blue  # bottom 1
        meshColors[6] = blue  # left 1
        meshColors[7] = blue  # right 1
        meshColors[8] = red  # wing 1
        meshColors[9] = red  # wing 2
        meshColors[10] = green  # tailwing 1
        meshColors[11] = green  # tailwing 2
        meshColors[12] = purple # tail 1
        return points, index, meshColors

