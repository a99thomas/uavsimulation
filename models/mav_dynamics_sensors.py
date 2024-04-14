"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
mavsim_python
    - Beard & McLain, PUP, 2012
    - Update history:  
        2/24/2020 - RWB
"""
import numpy as np
from message_types.msg_sensors import MsgSensors
import parameters.aerosonde_parameters as MAV
import parameters.sensor_parameters as SENSOR
from models.mav_dynamics_control import MavDynamics as MavDynamicsNoSensors
from tools.rotations import quaternion_to_rotation, quaternion_to_euler, euler_to_rotation
from numpy import sin, cos
from math import radians, degrees
import time

class MavDynamics(MavDynamicsNoSensors):
    def __init__(self, Ts):
        super().__init__(Ts)
        # initialize the sensors message
        self._sensors = MsgSensors()
        # random walk parameters for GPS
        self._gps_eta_n = 0.
        self._gps_eta_e = 0.
        self._gps_eta_h = 0.
        # timer so that gps only updates every ts_gps seconds
        self._t_gps = 999.  # large value ensures gps updates at initial time.

    def sensors(self):
        "Return value of sensors on MAV: gyros, accels, absolute_pressure, dynamic_pressure, GPS"
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        # simulate rate gyros(units are rad / sec)
        self._sensors.gyro_x = self._state.item(10) + SENSOR.gyro_x_bias + np.random.normal(0., SENSOR.gyro_sigma)
        self._sensors.gyro_y = self._state.item(11) + SENSOR.gyro_y_bias + np.random.normal(0., SENSOR.gyro_sigma)
        self._sensors.gyro_z = self._state.item(12) + SENSOR.gyro_z_bias + np.random.normal(0., SENSOR.gyro_sigma)

        # simulate accelerometers(units of g)
        # print(self._forces.item(0))
        self._sensors.accel_x = self._forces.item(0) / MAV.mass + MAV.gravity * np.sin(theta) + np.random.normal(0.,SENSOR.accel_sigma)
        self._sensors.accel_y = self._forces.item(1) / MAV.mass - MAV.gravity * np.cos(theta) * np.sin(phi) + np.random.normal(0.,SENSOR.accel_sigma)
        self._sensors.accel_z = self._forces.item(2) / MAV.mass - MAV.gravity * np.cos(theta) * np.cos(phi) + np.random.normal(0.,SENSOR.accel_sigma)

        # simulate magnetometers
        # magnetic field in provo has magnetic declination of 12.5 degrees
        # and magnetic inclination of 66 degrees
        Rm_i = euler_to_rotation(0, radians(-66), radians(12.5))
        # Rm_i = euler_to_rotation(0,0,0)
        m_i = Rm_i @ [[1],[0],[0]]
        # print(m_i)
        m_b = euler_to_rotation(phi, theta, psi).T @ m_i
        # print(m_b)
        R_b_v1 = np.array([[cos(theta), sin(phi)*sin(theta), sin(theta)*cos(phi)],
                           [0, cos(phi), -sin(phi)],
                           [-sin(theta), cos(theta)*sin(phi), cos(theta)*cos(phi)]])
        
        y = m_b + np.random.normal(0.,SENSOR.mag_sigma) + SENSOR.mag_beta
        self._sensors.mag_x = y[0]
        self._sensors.mag_y = y[1]
        self._sensors.mag_z = y[2]
        # print(y)
        m_v1 = R_b_v1 @ y
        psi_m = -np.arctan2(m_v1[1], m_v1[0])
        true_psi_m = psi_m + radians(12.5)
        # print("Magnetometer magnetic heading:",round(degrees(psi_m),3), "deg -- Magnetometer True heading:", round(degrees(true_psi_m[0]),3), "deg")
        # print(self._sensors.mag_z)

        # simulate pressure sensors
        self._sensors.abs_pressure = MAV.rho * MAV.gravity * (-self._state.item(2)) + np.random.normal(0., SENSOR.abs_pres_sigma)
        self._sensors.diff_pressure = MAV.rho * self._Va**2 / 2 + np.random.normal(0., SENSOR.diff_pres_sigma)
        self._sensors.theta = self.true_state.theta
        self._sensors.phi = self.true_state.phi
        self._sensors.psi = self.true_state.psi
        
        # simulate GPS sensor
        if self._t_gps >= SENSOR.ts_gps:
            self._gps_eta_n = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_n + np.random.normal(0.,SENSOR.gps_n_sigma)
            self._gps_eta_e = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_e + np.random.normal(0.,SENSOR.gps_e_sigma)
            self._gps_eta_h = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_h + np.random.normal(0.,SENSOR.gps_h_sigma)
            self._sensors.gps_n = self._state.item(0) + self._gps_eta_n
            self._sensors.gps_e = self._state.item(1) + self._gps_eta_e
            self._sensors.gps_h = -self._state.item(2) + self._gps_eta_h
            self._sensors.gps_Vg = np.sqrt((self._Va*np.cos(psi) + self._wind.item(0))**2 + (self._Va*np.sin(psi) + self._wind.item(1))**2) + np.random.normal(0.,SENSOR.gps_Vg_sigma)
            self._sensors.gps_course = np.arctan2(self._Va * np.sin(psi) + self._wind.item(1), self._Va*np.cos(psi) + self._wind.item(0)) + np.random.normal(0.,SENSOR.gps_course_sigma)
            self._t_gps = 0#SENSOR.ts_gps
        else: 
            self._t_gps += self._ts_simulation
        return self._sensors

    def external_set_state(self, new_state):
        self._state = new_state

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        pdot = quaternion_to_rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self.true_state.theta * 180/np.pi#self._alpha 
        self.true_state.beta = self.true_state.phi * 180/np.pi #self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta 
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
        self.true_state.bx = SENSOR.gyro_x_bias
        self.true_state.by = SENSOR.gyro_y_bias
        self.true_state.bz = SENSOR.gyro_z_bias
        self.true_state.camera_az = self._state.item(13)
        self.true_state.camera_el = self._state.item(14)