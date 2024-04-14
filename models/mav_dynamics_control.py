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
from models.mav_dynamics import MavDynamics as MavDynamicsForces
from math import exp
# load message types
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta
import parameters.aerosonde_parameters as MAV
from tools.rotations import quaternion_to_rotation, quaternion_to_euler, euler_to_rotation
from numpy import sin, cos, tan


class MavDynamics(MavDynamicsForces):
    def __init__(self, Ts):
        super().__init__(Ts)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # update velocity data and forces and moments
        self._update_velocity_data()
        self._forces_moments(delta=MsgDelta())
        # update the message class for the true state
        self._update_true_state()


    ###################################
    # public functions
    def update(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)
        super()._rk4_step(forces_moments)
        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)
        # update the message class for the true state
        self._update_true_state()

    ###################################
    # private functions
    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]
        # print("Wind", steady_state + gust)

        # phi, theta, psi = quaternion_to_euler(self._state[6:10])

        V = self._state[3:6]

        
        ##### TODO #####
        # convert wind vector from world to body frame (self._wind = ?)
        self._wind = quaternion_to_rotation(self._state[6:10]) @ (steady_state) + gust

        # velocity vector relative to the airmass ([ur , vr, wr]= ?)]
        urwr = V - self._wind
        # print("ur:", ur)
        ur, vr, wr = (urwr.item(i) for i in range(len(urwr)))
        # compute airspeed (self._Va = ?)\
        self._Va = np.sqrt(ur**2 + vr**2 + wr**2) 
        # print(np.sqrt(ur**2 + vr**2 + wr**2))
        # compute angle of attack (self._alpha = ?)
        self._alpha = np.arctan2(wr, ur)
        # print("alpha", self._alpha)

        # compute sideslip angle (self._beta = ?)
        # print(np.sqrt(ur**2 + vr**2 + wr**2))
        # self._beta = np.arctan2(vr, np.sqrt(ur**2 + wr**2))
        # print("var", ur, vr, wr)
        self._beta = np.arcsin(vr/self._Va)
        
        # print(self._beta)
        # print(self._alpha)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        ##### TODO ######
        # extract states (phi, theta, psi, p, q, r)
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # compute gravitational forces ([fg_x, fg_y, fg_z])
        # [fg_x, fg_y, fg_z] = MAV.mass * MAV.gravity * np.array([[-sin(theta)],
        #                                                         [cos(theta)*sin(phi)],
        #                                                         [cos(theta)*cos(phi)]])

        mass = MAV.mass
        g = MAV.gravity
        rho = MAV.rho
        # update the message class for the true state
        self._update_true_state()
        alpha = self._alpha
        beta = self._beta
        c = MAV.c

        # print(alpha)
        # compute Lift and Drag coefficients (CL, CD)
        AR = MAV.b**2 / MAV.S_wing
        # C_L_alpha = np.pi * AR / (1. + np.sqrt(1 + (AR/2.)**2))
        # CL = MAV.C_ell_0 + C_L_alpha * alpha
        M = MAV.M
        sigma = (1 + exp(-M*(alpha-MAV.alpha0)) + exp(M*(alpha+MAV.alpha0))) / \
                (1 + exp(-M*(alpha-MAV.alpha0)) * (1 + exp(M*(alpha+MAV.alpha0))))
        # print("sigma", sigma)
        CL = (1-sigma)*(MAV.C_L_0 + MAV.C_L_alpha*alpha) + sigma * (2*np.sign(alpha)*sin(alpha)**2*cos(alpha))

        e = MAV.e
        CD = MAV.C_D_p + (MAV.C_L_0 + MAV.C_L_alpha * alpha)**2 / (np.pi * e * AR)

        # compute Lift and Drag Forces (F_lift, F_drag)
        # F_lift = 0.5 * MAV.rho * self._Va**2 * MAV.S_wing * (CL + MAV.C)
        # F_lift = 0.5 * MAV.rho * self._Va**2 * MAV.S_wing * CL
        # F_drag = 0.5 * MAV.rho * self._Va**2 * MAV.S_wing * CD

        # propeller thrust and torque
        # print("throttle", delta.throttle)
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)
        # print("Thrust", thrust_prop.item(0))
        thrust_prop = float(thrust_prop)
        torque_prop = float(torque_prop)
        # compute longitudinal forces in body frame (fx, fz)
        Cx = -CD * cos(alpha) + CL * sin(alpha)
        Cx_q = - MAV.C_D_q * cos(alpha) + MAV.C_L_q * sin(alpha)
        Cx_del_e = - MAV.C_D_delta_e * cos(alpha) + MAV.C_L_delta_e * sin(alpha)
        Cz = -CD * sin(alpha) - CL * cos(alpha)
        Cz_q = - MAV.C_D_q * sin(alpha) - MAV.C_L_q * cos(alpha)
        Cz_del_e = - MAV.C_D_delta_e * sin(alpha) - MAV.C_L_delta_e* cos(alpha)
        f = np.array([[-mass * g * sin(theta)],
                      [mass*g*cos(theta)*sin(phi)],
                      [mass*g*cos(theta)*cos(phi)]]) + .5 * rho * self._Va**2 * MAV.S_wing * \
            np.array([[Cx + Cx_q * c * q/(2*self._Va)],
                      [MAV.C_Y_0 + MAV.C_Y_beta * beta + MAV. C_Y_p * beta/(2*self._Va) * p + MAV.C_Y_r*beta/(2*self._Va)*r ],
                      [Cz + Cz_q * c * q / (2*self._Va)]]) + .5 * rho * self._Va**2 * MAV.S_wing * \
            np.array([[Cx_del_e*delta.elevator],
                      [MAV.C_Y_delta_a * delta.aileron + MAV.C_Y_delta_r * delta.rudder],
                      [Cz_del_e * delta.elevator]]) + np.array([[thrust_prop, 0, 0]]).T


        # print("thrust:", thrust_prop)
        # compute logitudinal torque in body frame (My)
        Coeff = .5 * rho * self._Va**2 * MAV.S_wing
        b = MAV.b
        self._forces = f[0:3]

        moments = Coeff * np.array([[b * (MAV.C_ell_0 + MAV.C_ell_beta * beta + MAV.C_ell_p * p * b / (2*self._Va) + MAV.C_ell_r*b*r/(2*self._Va))],
                                    [c * (MAV.C_m_0 + MAV.C_m_alpha*alpha + MAV.C_m_q*c*q/(2*self._Va))],
                                    [b * (MAV.C_n_0 + MAV.C_n_beta * beta + MAV.C_n_p*b*p/(2*self._Va) + MAV.C_n_r*b*r/(2*self._Va))]]) \
                +  Coeff * np.array([[b * (MAV.C_ell_delta_a * delta.aileron + MAV.C_ell_delta_r * delta.rudder)],
                                       [c * MAV.C_m_delta_e * delta.elevator],
                                       [b * (MAV.C_n_delta_a*delta.aileron + MAV.C_n_delta_r*delta.rudder)]]) \
                - np.array([[torque_prop, 0, 0]]).T
        # print("beta", Coeff)
        # compute lateral torques in body frame (Mx, Mz)

        forces_moments = np.array([f[0], f[1], f[2], moments[0], moments[1], moments[2]])

        return forces_moments

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller
        ##### TODO #####
        # map delta_t throttle command(0 to 1) into motor input voltage
        v_in = MAV.V_max * delta_t

        # Angular speed of propeller (omega_p = ?)
        rho = MAV.rho
        D = MAV.D_prop
        a = rho * np.power(D, 5) / ((2*np.pi)**2) * MAV.C_Q0
        b = rho * np.power(D, 4) / (2*np.pi) * MAV.C_Q1 * Va + MAV.KQ*MAV.KV/MAV.R_motor
        c = rho * np.power(D, 3) * MAV.C_Q2 * Va**2 - MAV.KQ * v_in / MAV.R_motor + MAV.KQ * MAV.i0
        omega_p = (-b + np.sqrt(b**2 - 4*a*c)) / (2.*a)
        J = 2 * np.pi * Va / (omega_p * D)
        CT = MAV.C_T2 * J**2 + MAV.C_T1 * J + MAV.C_T0
        CQ = MAV.C_Q2 * J**2 + MAV.C_Q1 * J + MAV.C_Q0

        n = omega_p / (2*np.pi)
        thrust_prop = rho * n**2 * np.power(D, 4) * CT
        torque_prop = rho * n**2 * np.power(D, 5) * CQ

        return thrust_prop, torque_prop

    def _update_true_state(self):
        # rewrite this function because we now have more information
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        pdot = quaternion_to_rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
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
        self.true_state.bx = 0
        self.true_state.by = 0
        self.true_state.bz = 0
        self.true_state.camera_az = 0
        self.true_state.camera_el = 0
