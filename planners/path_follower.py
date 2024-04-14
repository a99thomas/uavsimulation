import numpy as np
from math import sin, cos
from message_types.msg_autopilot import MsgAutopilot
from tools.wrap import wrap


class PathFollower:
    def __init__(self):
        ##### TODO #####
        self.chi_inf = np.radians(50)  # approach angle for large distance from straight-line path
        self.k_path = 0.01  # path gain for straight-line path following
        self.k_orbit = 1.  # path gain for orbit following
        self.gravity = 9.81
        self.autopilot_commands = MsgAutopilot()  # message sent to autopilot

    def update(self, path, state):
        if path.type == 'line':
            self._follow_straight_line(path, state)
        elif path.type == 'orbit':
            self._follow_orbit(path, state)
        return self.autopilot_commands

    def _follow_straight_line(self, path, state):
        ##### TODO #####
        #airspeed command
        chi_q = np.arctan2(path.line_direction[1],path.line_direction[0])
        self.autopilot_commands.airspeed_command = path.airspeed
        chi_q = wrap(chi_q, state.chi)[0]
        # print(chi_q)

        # course command
        Ri_P = np.array([[np.cos(chi_q), np.sin(chi_q), 0],
                         [-np.sin(chi_q), np.cos(chi_q), 0],
                         [0, 0, 1]])
        p = np.array([[state.north, state.east, state.altitude]]).T
        r = path.line_origin
        self.e_p = Ri_P @ (p - r)

        # self.chi_inf
        # print(self.e_p)
        # print("hi", chi_q - self.chi_inf * (2/np.pi) * np.arctan(self.k_path * self.e_p[1]))
        chi_d = - self.chi_inf * (2/np.pi) * np.arctan(self.k_path * self.e_p[1])
        # print(chi_d)
        self.autopilot_commands.course_command = chi_q + chi_d[0]

        # altitude command
        k = np.array([[0, 0, 1]])
        q = path.line_direction.T
        n = np.cross(k, q)[0,:] / np.linalg.norm(np.cross(k, q))
        s = self.e_p[:,0] - np.dot(self.e_p[:,0], n) * n
        hd = -r[2] - np.sqrt(s[0]**2 + s[1]**2) * (q[0,2]/np.sqrt(q[0,0]**2 + q[0,1]**2))
        self.autopilot_commands.altitude_command = hd[0]
        # feedforward roll angle for straight line is zero
        self.autopilot_commands.phi_feedforward = 0

    def _follow_orbit(self, path, state):
        ##### TODO #####
        # airspeed command
        if path.orbit_direction == "CW":
            direction = 1.0
        else:
            direction = -1.0
        c = path.orbit_center
        rho = path.orbit_radius
        d = np.sqrt((state.north - c.item(0))**2 + (state.east - c.item(1))**2)
        varphi = np.arctan2(state.east - c.item(1), state.north - c.item(0))
        varphi = wrap(varphi, state.chi)
        self.autopilot_commands.airspeed_command = path.airspeed

        orbit_error = (d - rho)/rho
        # print("orbit_error", orbit_error)
        # course command
        self.autopilot_commands.course_command = varphi + direction * (np.pi/2. + np.arctan(self.k_orbit * orbit_error))

        # altitude command
        self.autopilot_commands.altitude_command = -path.orbit_center.item(2)
        
        # roll feedforward command
        if orbit_error < 0.015: 
        
            self.autopilot_commands.phi_feedforward = direction * np.arctan(state.Vg **2 / self.gravity / path.orbit_radius / np.cos(state.chi - state.psi))
            
        else:
            self.autopilot_commands.phi_feedforward = 0


