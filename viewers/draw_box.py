"""
mavsim_python: drawing tools
Aaron Thomas 1/16/24
"""
import numpy as np
import pyqtgraph.opengl as gl
from tools.rotations import euler_to_rotation
from tools.drawing import rotate_points, translate_points, points_to_mesh


class DrawBox:
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
        box_position = np.array([[state.north], [state.east], [-state.altitude]])  # NED coordinates
        # attitude of box as a rotation matrix R from body to inertial
        R_bi = euler_to_rotation(state.phi, state.theta, state.psi)
        # convert North-East Down to East-North-Up for rendering
        self.R_ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        # get points that define the non-rotated, non-translated spacecraft and the mesh colors
        self.box_points, self.box_index, self.box_meshColors = self.get_box_points()
        self.box_body = self.add_object(
            self.box_points,
            self.box_index,
            self.box_meshColors,
            R_bi,
            box_position)
        window.addItem(self.box_body)  # add spacecraft to plot     

    def update(self, state):
        box_position = np.array([[state.north], [state.east], [-state.altitude]])  # NED coordinates
        # attitude of box as a rotation matrix R from body to inertial
        R_bi = euler_to_rotation(state.phi, state.theta, state.psi)
        self.box_body = self.update_object(
            self.box_body,
            self.box_points,
            self.box_index,
            self.box_meshColors,
            R_bi,
            box_position)

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

    def get_box_points(self, length=2., width=2., height=0.25):
        """"
            Points that define the spacecraft, and the colors of the triangular mesh
            Define the points on the spacecraft following information in Appendix C.3
        """
        # points are in XYZ coordinates
        #   define the points on the spacecraft according to Appendix C.3
        points = self.unit_length * np.array([
            [-length/2, -width/2, -height/2],  # point 1 [0]
            [length/2, -width/2, -height/2],  # point 2 [1]
            [-length/2, width/2, -height/2],  # point 3 [2]
            [length/2, width/2, -height/2],  # point 4 [3]
            [-length/2, -width/2, height/2],  # point 5 [4]
            [length/2, -width/2, height/2],  # point 6 [5]
            [-length/2, width/2, height/2],  # point 7 [6]
            [length/2, width/2, height/2]  # point 8 [7]
            ]).T
        # point index that defines the mesh
        index = np.array([
            [0, 1, 2],  # Bottom 1
            [1, 2, 3],  # Bottom 2
            [0, 1, 5],  # Side 1a
            [0, 4, 5],  # Side 1b
            [0, 2, 6],  # Side 2a
            [0, 4, 6],  # Side 2b
            [1, 3, 5], # Side 3a
            [1, 5, 7], # Side 3b
            [2, 3, 7],  # Side 4a
            [2, 6, 7],  # Side 4b
            [4, 5, 6],  # Top 1
            [5, 6, 7],  # Top 2
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
        meshColors[2] = blue  # front 3
        meshColors[3] = blue  # front 4
        meshColors[4] = blue  # top 1
        meshColors[5] = blue  # bottom 1
        meshColors[6] = blue  # left 1
        meshColors[7] = blue  # right 1
        meshColors[8] = blue  # wing 1
        meshColors[9] = blue  # wing 2
        meshColors[10] = red  # tailwing 1
        meshColors[11] = red  # tailwing 2
        return points, index, meshColors

