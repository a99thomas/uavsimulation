"""
mavsim_python: drawing tools
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - RWB
        3/30/2022 - RWB
        7/13/2023 - RWB
        3/25/2024 - RWB
"""

import numpy as np
from message_types.msg_path import MsgPath
from tools.rotations import euler_to_rotation
import time


class Delivery_Zone:
    def __init__(self, location):
        self._path1 = MsgPath()
        self._path2 = MsgPath()

        drop_location = np.array([location[0],location[1], 0])
        # self._delivery_location = location
        self.construct_line(drop_location)
        # self.rewind = False


    def construct_line(self, location):
        self._path1.plot_updated = False
        self._path1.type = 'line'  
        self._path1.airspeed = 0, 
        self._path1.line_origin = location - np.array([[20, 20, 0]])
        self._path1.line_direction = np.array([[0.04], [0.04], [0]])

        self._path2.plot_updated = False
        self._path2.type = 'line'  
        self._path2.airspeed = 0, 
        self._path2.line_origin = location - np.array([[20, -20, 0]])
        self._path2.line_direction = np.array([[0.04], [-0.04], [0]])
        