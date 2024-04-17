"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
mavsim_python
    - Beard & McLain, PUP, 2012
    - Update history:  
        2/24/2020 - RWB
        7/13/2023 - RWB
        1/17/2024 - RWB
"""
import numpy as np
# load message types
from message_types.msg_state import MsgState
import parameters.aerosonde_parameters as MAV
from tools.rotations import quaternion_to_rotation, quaternion_to_euler, euler_to_rotation
from numpy import sin, cos, tan, arcsin, arccos, arctan

class BoxDynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        self.released = False
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.north0],  # (0)
                               [MAV.east0],   # (1)
                               [MAV.down0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0],    # (12)
                               [0],   # (13)
                               [0],   # (14)
                               ])
        # initialize true_state message
        self.true_state = MsgState()
        self.true_state.Vdrop = 0 
        self.prev_sim_time = 0

    ###################################
    # public functions
    def update(self, box_state, mav_state, waypoints, sim_time):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        if self.released:
            self.calculate_position(sim_time)
        else:
            box_offset = (euler_to_rotation(self.true_state.phi, self.true_state.theta, self.true_state.psi) @ np.array([[0, 0, 8]]).T)

            self._state[0] = mav_state.north + box_offset[0]
            self._state[1] = mav_state.east + box_offset[1]
            self._state[2]= -mav_state.altitude + box_offset[2]
            box_state.north = mav_state.north
            box_state.east = mav_state.east
            box_state.altitude = -mav_state.altitude + 6
            box_state.chi = mav_state.chi
            box_state.alpha = mav_state.alpha
            box_state.Va = mav_state.Va
            self.true_state.chi = mav_state.chi
            self.true_state.phi = mav_state.phi
            self.true_state.psi = mav_state.psi
            self.true_state.theta = mav_state.theta
            self.true_state.beta = mav_state.beta
            self.true_state.alpha = mav_state.alpha
            self.true_state.gamma = mav_state.gamma
            self._determine_release(mav_state, waypoints, sim_time)

            # self_state[3] = mav_state.Va
        # update the message class for the true state
        self._update_true_state()

        

    def _determine_release(self, mav_state, waypoints, sim_time):
        destination = np.reshape(waypoints.ned[0:3,-2],(3,1))
        self.destination = destination
        prev_waypoint = np.reshape(waypoints.ned[0:3,-3],(3,1))
        Va = mav_state.Va
        # time_before_destination = np.roots([4.905, Va * np.sin(self.true_state.theta), -self.true_state.altitude])[1]
        time_before_destination = (-Va * np.sin(self.true_state.theta) + np.sqrt((Va * np.sin(self.true_state.theta))**2 - 4*4.905 * -self.true_state.altitude))/9.81
        self.vel_at_drop = np.reshape([[Va * np.cos(self.true_state.chi)],[Va*np.
                                                                      sin(self.true_state.chi)], [Va*np.sin(self.true_state.theta)]], (3,1))
        dropping_point = destination - time_before_destination * self.vel_at_drop*1.06
        # course_on_drop = waypoints.course.item(-1)
        self.true_state.Vdrop = self.vel_at_drop[2]

        self.halfspace_r = np.reshape(dropping_point,(3,1))
        self.halfspace_n = np.reshape((destination - prev_waypoint ) / np.linalg.norm(destination - prev_waypoint),(3,1))
        mav_pos = np.array([[mav_state.north, mav_state.east, -mav_state.altitude]]).T
        if self._inHalfSpace(mav_pos) and np.linalg.norm(mav_pos - destination) <= 400:
            self.released = True
            self.drop_time = sim_time
            self.prev_sim_time = sim_time
            self.height_at_drop = -self.true_state.altitude
            self.Vdrop = self.vel_at_drop.item(2)
            # print(self.height_at_drop)
            # print(self.vel_at_drop)
            print("Dropping")

    def calculate_position(self, sim_time):
        time_since_drop = sim_time - self.drop_time
        time_step = sim_time - self.prev_sim_time
        altitude = -self._state[2]
        if self._state[2] < 0:
            # print(sim_time)
            # print("Vel", self.vel_at_drop)
            # prin't("step",time_step)
            self._state[0] += self.vel_at_drop[0] * time_step
            self._state[1] += self.vel_at_drop[1] * time_step
            # altitude += self.Vdrop*time_step - 4.905*time_step**2
            # print("height",-self.height_at_drop)
            # print(self.Vdrop*time_since_drop)
            # print(-4.905*time_since_drop**2)
            altitude = -self.height_at_drop + self.Vdrop*time_since_drop - 4.905*time_since_drop**2
            # self.Vdrop += -9.81*time_step
            self._state[2] = -altitude
            # print(self._state[0:3])
            if self._state[2] > 0:
                print("Pizza delivered", np.round(np.linalg.norm(self._state[0:2] - self.destination[0:2]),1), "feet from the target")
        else:
            self._state[2] = 0

        self.prev_sim_time = sim_time
        # print("hi")


    def update2(self, box_state, mav_state, sim_time, trigger):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        if self.released:
            self.calculate_position2(sim_time)
        else:
            box_offset = (euler_to_rotation(self.true_state.phi, self.true_state.theta, self.true_state.psi) @ np.array([[0, 0, 8]]).T)

            self._state[0] = mav_state.north + box_offset[0]
            self._state[1] = mav_state.east + box_offset[1]
            self._state[2]= -mav_state.altitude + box_offset[2]
            box_state.north = mav_state.north
            box_state.east = mav_state.east
            # box_state.altitude = mav_state.altitude + 6
            box_state.chi = mav_state.chi
            box_state.alpha = mav_state.alpha
            box_state.Va = mav_state.Va
            self.true_state.chi = mav_state.chi
            self.true_state.phi = mav_state.phi
            self.true_state.psi = mav_state.psi
            self.true_state.theta = mav_state.theta
            self.true_state.beta = mav_state.beta
            self.true_state.alpha = mav_state.alpha
            self.true_state.gamma = mav_state.gamma
            # print(self.true_state.phi)
            Va = mav_state.Va
            self.vel_at_drop = np.reshape([[Va * np.cos(self.true_state.chi)],[Va*np.
                                                                      sin(self.true_state.chi)], [Va*np.sin(self.true_state.theta)]], (3,1))
            # course_on_drop = waypoints.course.item(-1)
            self.true_state.Vdrop = self.vel_at_drop[2]
            if trigger == True:
                self.released = True
                self.drop_time = sim_time
                self.height_at_drop = -self.true_state.altitude
                self.Vdrop = self.vel_at_drop.item(2)
                self.prev_sim_time = sim_time

            # self_state[3] = mav_state.Va
        # update the message class for the true state
        self._update_true_state()

    def calculate_position2(self, sim_time):
        time_since_drop = sim_time - self.drop_time
        time_step = sim_time - self.prev_sim_time
        altitude = -self._state[2]
        if self._state[2] < 0:
            # print(sim_time)
            # print("Vel", self.vel_at_drop)
            # prin't("step",time_step)
            self._state[0] += self.vel_at_drop[0] * time_step
            self._state[1] += self.vel_at_drop[1] * time_step
            # altitude += self.Vdrop*time_step - 4.905*time_step**2
            # print("height",-self.height_at_drop)
            # print(self.Vdrop*time_since_drop)
            # print(-4.905*time_since_drop**2)
            altitude = -self.height_at_drop + self.Vdrop*time_since_drop - 20*4.905*time_since_drop**2
            # self.Vdrop += -9.81*time_step
            self._state[2] = -altitude
            # print(self._state[0:3])
            if self._state[2] > 0:
                print("Pizza delivered")
        else:
            self._state[2] = 0

        self.prev_sim_time = sim_time

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = quaternion_to_euler(self._state[6:10])

        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = np.linalg.norm(self._state[3:5])
        pdot = quaternion_to_rotation(self._state[6:10]) @ self._state[3:6]
        # self.true_state.alpha = 0
        # self.true_state.beta = 0
        # self.true_state.phi = phi
        # self.true_state.theta = 0
        # self.true_state.psi = psi
        self.true_state.Vg = 0
        # self.true_state.gamma = 0
        # self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = 0
        self.true_state.we = 0
        self.true_state.bx = 0
        self.true_state.by = 0
        self.true_state.bz = 0
        self.true_state.camera_az = 0
        self.true_state.camera_el = 0

    def _inHalfSpace(self, 
                     pos: np.ndarray)->bool:
        
        '''Is pos in the half space defined by r and n?'''
        # print((pos-self._halfspace_r).T @ self._halfspace_n)
        if (pos-self.halfspace_r).T @ self.halfspace_n >= 0:
            return True
        else:
            return False