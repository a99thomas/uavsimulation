"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
"""
import numpy as np
from scipy import stats
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
from tools.wrap import wrap
from message_types.msg_state import MsgState
from message_types.msg_sensors import MsgSensors
from numpy import sin, cos, tan 

class Observer:
    def __init__(self, ts_control, initial_measurements = MsgSensors()):
        # initialized estimated state message
        self.estimated_state = MsgState()
        # use alpha filters to low pass filter gyros and accels
        # alpha = Ts/(Ts + tau) where tau is the LPF time constant

        ##### TODO #####
        self.lpf_gyro_x = AlphaFilter(alpha=0.9, y0=initial_measurements.gyro_x)
        self.lpf_gyro_y = AlphaFilter(alpha=0.9, y0=initial_measurements.gyro_y)
        self.lpf_gyro_z = AlphaFilter(alpha=0.9, y0=initial_measurements.gyro_z)
        self.lpf_accel_x = AlphaFilter(alpha=0.9, y0=initial_measurements.accel_x)
        self.lpf_accel_y = AlphaFilter(alpha=0.9, y0=initial_measurements.accel_y)
        self.lpf_accel_z = AlphaFilter(alpha=0.9, y0=initial_measurements.accel_z)
        # use alpha filters to low pass filter absolute and differential pressure
        self.lpf_abs = AlphaFilter(alpha=0.9, y0=initial_measurements.abs_pressure)
        self.lpf_diff = AlphaFilter(alpha=0.9, y0=initial_measurements.diff_pressure)
        # ekf for phi and theta
        self.attitude_ekf = EkfAttitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = EkfPosition()

    def update(self, measurement):
        ##### TODO #####
        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurement.gyro_x)
        self.estimated_state.q = self.lpf_gyro_y.update(measurement.gyro_y)
        self.estimated_state.r = self.lpf_gyro_z.update(measurement.gyro_z)

        # invert sensor model to get altitude and airspeed
        self.estimated_state.altitude = self.lpf_abs.update(measurement.abs_pressure)/(CTRL.rho * CTRL.gravity)
        self.estimated_state.Va = np.sqrt(2 * self.lpf_diff.update(measurement.diff_pressure)/ CTRL.rho)
        
        # estimate phi and theta with simple ekf
        self.attitude_ekf.update(measurement, self.estimated_state)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(measurement, self.estimated_state)
        self.estimated_state.alpha =  self.estimated_state.theta * 180/np.pi
        self.estimated_state.beta = self.estimated_state.phi * 180/np.pi
        self.estimated_state.bx = 0   #measurement.psi * 180/np.pi
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0

        return self.estimated_state


class AlphaFilter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.5, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        ##### TODO #####
        self.y = u * (1 - self.alpha) + self.y * (self.alpha)
        return self.y


class EkfAttitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    def __init__(self):
        ##### TODO #####
        self.Q = np.diag([1, 0.5])
        self.Q_gyro = np.diag([0.1, 0.1, 0.1])
        self.R_accel = np.diag([1., 1., 1.])
        self.N = 3 # number of prediction step per sample
        self.xhat = np.array([[0.0], [0.0]]) # initial state: phi, theta
        self.P = np.diag([0.008, 0.004])
        self.Ts = SIM.ts_control/self.N
        self.gate_threshold = stats.chi2.isf(q=0.01, df = 3)
        # print("grapes", self.gate_threshold)

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        ##### TODO #####
        # xdot = np.zeros((2,1))
        xdot = np.array([[state.p + state.q * sin(state.phi) * sin(state.theta) + state.r * cos(state.phi) * tan(state.theta)],
                      [state.q * cos(state.phi) + state.r * cos(state.phi)]])
        return xdot

    def h(self, x, measurement, state):
        # measurement model y
        ##### TODO #####
        h_ = np.array([[state.q * state.Va * sin(state.theta) + CTRL.gravity * sin(state.theta)],  # x-accel
                        [state.r * state.Va * cos(state.theta) - state.p * state.Va * sin(state.theta) - CTRL.gravity * cos(state.theta) * sin(state.phi)],# y-accel
                        [-state.q * state.Va * cos(state.theta) - CTRL.gravity * cos(state.theta) * cos(state.phi)]])  # z-accel
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        ##### TODO #####
        # x = np.array([state.phi,state.theta])
        self.x_hat = np.zeros(2)
        Tout = self.Ts * 3
        # Tp = self.Ts
        for i in range(0, self.N):
            Tp = Tout/self.N
            # self.P = state.
            self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)
            # A = jacobian(self.f, self.xhat, measurement, state)
            A = np.array([[state.q * cos(state.phi) * tan(state.theta) - state.r * sin(state.phi)*tan(state.theta), state.q * sin(state.phi) - state.r * cos(state.phi)],
                          [-state.q * sin(state.phi) - state.r * cos(state.phi), 0]])
            
            self.P = self.P + (Tp/self.N) * (A @ self.P + self.P @ A.T + self.Q)
            # print("P2:", self.P)
            if np.any(np.isnan(self.P)):
                quit()
            # np.zeros((2,2))

    def measurement_update(self, measurement, state):
        # measurement updates
        h = self.h(self.xhat, measurement, state)
        # C = jacobian(self.h, self.xhat, measurement, state)
        C = np.array([[0, state.q * state.Va * cos(state.theta) + CTRL.gravity*cos(state.theta)],
                      [-CTRL.gravity * cos(state.phi) * cos(state.theta), -state.r * state.Va * sin(state.theta) - state.p * state.Va * cos(state.theta) + CTRL.gravity * sin(state.phi) * sin(state.theta)],
                      [CTRL.gravity * sin(state.phi) * cos(state.theta), (state.q * state.Va + CTRL.gravity * cos(state.phi))*sin(state.theta)]])
        y = np.array([[measurement.accel_x, measurement.accel_y, measurement.accel_z]]).T
       
        ##### TODO #####
        S_inv = np.linalg.inv(self.R_accel + C @ self.P @ C.T)
        
        if (y-h).T @ S_inv @ (y-h) < self.gate_threshold:
            L = self.P @ C.T @ S_inv
            tmp = np.eye(2) - L @ C
            self.P = tmp @ self.P @ tmp.T + L @ self.R_accel @ L.T
            self.xhat =  self.xhat + L @ (y-h)


class EkfPosition:
    # implement continous-discrete EKF to estimate pn, pe, Vg, chi, wn, we, psi

    def __init__(self):
        self.Q = np.diag([
                    0.3**2, # pn
                    0.3**2, # pe
                    0.3**2, # Vg
                    0.1**2, # chi
                    .01**2, # wn
                    .01**2, # we
                    0.1**2, # psi
                    ])
        self.P = np.diag([
                    0.1**2, # pn
                    0.1**2, # pe
                    0.1**2, # Vg
                    0.1**2, # chi
                    0.1**2, # wn
                    0.1**2, # we
                    0.1**2, # psi
                    ])
        self.xhat = np.array([
            [0.0], # pn
            [0.0], # pe
            [0.1], # Vg
            [0.0], # chi
            [0.0], # wn
            [0.0], # we
            [0.0], # psi
            ])
        self.R_gps = np.diag([
                    0.1**2,  # y_gps_n
                    0.1**2,  # y_gps_e
                    0.1**2,  # y_gps_Vg
                    0.1**2,  # y_gps_course
                    ])
        self.R_pseudo = np.diag([
                    1.**2,  # pseudo measurement #1
                    1.**2,  # pseudo measurement #2
                    ])
        self.N = 1  # number of prediction step per sample
        self.Ts = SIM.ts_control / self.N
        self.gps_n_old = 0
        self.gps_e_old = 0
        self.gps_Vg_old = 0
        self.gps_course_old = 0
        self.pseudo_threshold = stats.chi2.isf(q=0.00001, df=2)
        self.gps_threshold = 100000 # don't gate GPS
        

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.north = self.xhat.item(0)
        state.east = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        xdot = np.array([[0],
                       [0],
                       [0],
                       [0],
                       [0.0],
                       [0.0],
                       [0],
                       ])
        Va = state.Va
        q = state.q #measurement.gyro_y
        r = state.r #measurement.gyro_z
        pn, pe, Vg, chi, wn, we, psi = (x.item(i) for i in range(len(x)))
        if Vg == 0:
            Vg = 0.01

        psidot = (q*sin(state.phi) + r * cos(state.phi)) / cos(state.theta)
        # Vgdot = Va * psidot *(we * cos(psi) - wn * sin(psi))/Vg
        Vgdot = ((Va * cos(psi) + wn) * (-psidot * Va * sin(psi)) + (Va * sin(psi) + we) * (psidot * Va * cos(psi))) / Vg
        xdot = np.array([[Vg * cos(chi)],
                         [Vg * sin(chi)],
                         [Vgdot],
                         [(CTRL.gravity / Vg) * tan(state.phi)],
                         [0.0],
                         [0.0],
                         [psidot]])       
        return xdot

    def h_gps(self, x, measurement, state):
        # measurement model for gps measurements y=h(x,u)
        y = np.array([
            [x.item(0)], #pn
            [x.item(1)], #pe
            [x.item(2)], #Vg
            [x.item(3)], #chi
        ])
        return y

    def h_pseudo(self, x, measurement, state):
        # measurement model for wind triangale pseudo measurement y=h(x,u)
        pn, pe, Vg, chi, wn, we, psi = (x.item(i) for i in range(len(x)))
        y = np.array([
            [state.Va * cos(psi) + wn - Vg * cos(chi)],  # wind triangle x
            [state.Va * sin(psi) + wn - Vg * sin(chi)],  # wind triangle y
        ])
        
        return y

    def propagate_model(self, measurement, state):
        # model propagation
        # self.xhat = np.zeros((7,1))
        for i in range(0, self.N):
            # propagate model
            # print("1")
            # print(self.xhat[3,0])
            Tp = self.Ts/self.N
            self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)

            # compute Jacobian
            A = jacobian(self.f, self.xhat, measurement, state)
            # convert to discrete time models
            Ad = np.eye(7) + A * Tp + A @ A * Tp**2
            # print("Ad = ", np.round(Ad,2))

            # update P with discrete time model
            # self.P = np.zeros((7,7))
            self.P = Ad @ self.P @ Ad.T + Tp**2 * self.Q


    def measurement_update(self, measurement, state):
        # always update based on wind triangle pseudo measurement
        yhat = self.h_pseudo(self.xhat, measurement, state)
        C = jacobian(self.h_pseudo, self.xhat, measurement, state)
        y = np.array([[0, 0]]).T
        S_inv = np.linalg.inv(self.R_pseudo + C @ self.P @ C.T)
        # print("Sinv", S_inv)
        # print(stats.chi2.isf(q=0.001, df=2))
        if (y-yhat).T @ S_inv @ (y-yhat) < self.pseudo_threshold:
        # if True:
            L = self.P @ C.T @ S_inv
            LC = L @ C 
            self.P = (np.eye(LC.shape[0]) - LC) @ self.P @ (np.eye(LC.shape[0]) - LC).T + L @ self.R_pseudo @ L.T
            # print("P1", self.P)
            # print(self.xhat[3,0])
            self.xhat = self.xhat + L @ (y-yhat)

        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):

            yhat = self.h_gps(self.xhat, measurement, state)
            C = jacobian(self.h_gps, self.xhat, measurement, state)
            y_chi = wrap(measurement.gps_course, yhat[3, 0]) ###Big problem here with yhat
            y = np.array([[measurement.gps_n,
                           measurement.gps_e,
                           measurement.gps_Vg,
                           y_chi]]).T
            S_inv = np.linalg.inv(self.R_gps + C @ self.P @ C.T)
            if (y-yhat).T @ S_inv @ (y-yhat) < self.gps_threshold:
                L = self.P @ C.T @ S_inv
                # print("P12", self.P)
                self.xhat = self.xhat + L @ (y-yhat)
                LC = L @ C 
                self.P = (np.eye(LC.shape[0]) - LC) @ self.P @ (np.eye(LC.shape[0]) - LC).T + L @ self.R_gps @ L.T
            

            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course








# class EkfPosition:
#     # implement continous-discrete EKF to estimate pn, pe, Vg, chi, wn, we, psi
#     def __init__(self):
#         self.Q = np.diag([
#                     .0,  # pn
#                     .0,  # pe
#                     .1,  # Vg
#                     0.1, # chi
#                     0.0, # wn
#                     0.0, # we
#                     0.1, # psi
#                     ])
#         self.P = np.diag([
#                     0.1**2, # pn
#                     0.1**2, # pe
#                     0.1**2, # Vg
#                     0.1**2, # chi
#                     0.1**2, # wn
#                     0.1**2, # we
#                     0.1**2, # psi
#                     ])
#         self.xhat = np.array([
#             [0.0], # pn
#             [0.0], # pe
#             [0.0], # Vg
#             [0.0], # chi
#             [0.0], # wn
#             [0.0], # we
#             [0.0], # psi
#             ])
#         self.R_gps = np.diag([
#                     1**2,  # y_gps_n
#                     0.1**2,  # y_gps_e
#                     0.1**2,  # y_gps_Vg
#                     0.1**2,  # y_gps_course
#                     ])
#         self.R_pseudo = np.diag([
#                     0.1**2,  # pseudo measurement #1
#                     0.1**2,  # pseudo measurement #2
#                     ])
#         self.N = 1  # number of prediction step per sample
#         self.Ts = (SIM.ts_control / self.N)
#         self.gps_n_old = 0
#         self.gps_e_old = 0
#         self.gps_Vg_old = 0.001
#         self.gps_course_old = 0
#         self.pseudo_threshold = stats.chi2.isf(q=0.01, df=7)
#         self.gps_threshold = 100000 # don't gate GPS

#     def update(self, measurement, state):
#         self.propagate_model(measurement, state)
#         self.measurement_update(measurement, state)
#         state.north = self.xhat.item(0)
#         # print(state.north)
#         state.east = self.xhat.item(1)
#         state.Vg = self.xhat.item(2)
#         state.chi = measurement.gps_course#self.xhat.item(3)
#         state.wn = self.xhat.item(4)
#         state.we = self.xhat.item(5)
#         state.psi = self.xhat.item(6)


#     def f(self, x, measurement, state):
#         # system dynamics for propagation model: xdot = f(x, u)
#         f_ = np.array([[0],
#                        [0],
#                        [0],
#                        [0],
#                        [0.0],
#                        [0.0],
#                        [0],
#                        ])
#         Va = state.Va
#         q = state.q#measurement.gyro_y
#         r = state.r#measurement.gyro_z
#         pn, pe, Vg, chi, wn, we, psi = (x.item(i) for i in range(len(x)))
#         if Vg == 0:
#             Vg = 0.01

#         psidot = (q*sin(state.phi) + r * cos(state.phi)) / cos(state.theta)
#         Vgdot = Va * psidot *(we * cos(psi) - wn * sin(psi))/Vg
#         # ((Va * cos(psi) + wn) * (-psidot * Va * sin(psi)) \
#         #          + (Va * sin(psi) + we) * (psidot * Va * cos(psi))) / Vg
#         xdot = np.array([[Vg * cos(chi)],
#                          [Vg * sin(chi)],
#                          [Vgdot],
#                          [(CTRL.gravity / Vg) * tan(state.phi) * cos(chi - psi)],
#                          [0.0],
#                          [0.0],
#                          [psidot]])
#         f_ = xdot
#         return f_

#     def h_gps(self, x, measurement, state):
#         # measurement model for gps measurements
#         h_ = np.array([
#             # [measurement.gps_n], #pn
#             # [measurement.gps_e], #pe
#             # [measurement.gps_Vg], #Vg
#             # [measurement.gps_course], #chi
#             [x.item(0)],
#             [x.item(1)],
#             [x.item(2)],
#             [x.item(3)]
#         ])
#         return h_

#     def h_pseudo(self, x, measurement, state):
#         # measurement model for wind triangle pseudo measurement
#         pn, pe, Vg, chi, wn, we, psi = (x.item(i) for i in range(len(x)))
#         h_ = np.array([
#             [state.Va * cos(psi) + wn - Vg * cos(chi)],  # wind triangle x
#             [state.Va * sin(psi) + wn - Vg * sin(chi)],  # wind triangle y
#         ])
#         # print("H+ ", h_)
#         return h_

#     def propagate_model(self, measurement, state):
#         # model propagation

#         print("3")
#         self.xhat = np.zeros((7,1)) 
#         Tout = self.Ts 
#             # Tp = self.Ts
#         Tp = Tout/self.N

#         for i in range(0, self.N):
#             Tp = Tout/self.N
#             # self.P = state.
#             print("4")
#             self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)
#             # A = jacobian(self.f, self.xhat, measurement, state)
#             A = jacobian(self.f, self.xhat, measurement, state)
#             y_F = self.f(self.xhat, measurement, state)
#             psidot = y_F[6]
#             chidot = CTRL.gravity * tan(state.psi) * cos(state.chi - state.psi) / state.Vg
#             Vgdot = ((state.Va * cos(state.psi) + state.wn) * (state.Va * psidot * sin(state.psi)) + (state.Va * sin(state.psi) + state.we) * (state.Va * psidot * cos(state.psi)))/state.Vg
#             DVgdot_DPsi =  (-psidot * state.Va * (state.wn * cos(state.psi) + state.we * sin(state.psi)))/state.Vg
#             Dchidot_DVg = -CTRL.gravity * tan(state.phi) * cos(state.chi - state.psi) / (state.Vg**2)
#             Dchidot_Dchi = -CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi) / state.Vg
#             Dchidot_Dpsi = CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi)
            
#             self.A = np.array([[0., 0., cos(state.chi), -state.Vg * sin(state.chi), 0., 0., 0.],
#                             [0., 0., sin(state.chi), state.Vg * cos(state.chi), 0., 0., 0.],
#                             [0., 0., -Vgdot[0] / state.Vg, 0., -psidot[0] * state.Va * sin(state.psi), psidot[0] * state.Va * cos(state.psi), DVgdot_DPsi[0]],
#                             [0., 0., Dchidot_DVg, Dchidot_Dchi, 0., 0., Dchidot_Dpsi],
#                             [0., 0., 0., 0., 0., 0., 0.],
#                             [0., 0., 0., 0., 0., 0., 0.],
#                             [0., 0., 0., 0., 0., 0., 0.]])
#             # np.array([[state.q * cos(state.phi) * tan(state.theta) - state.r * sin(state.phi)*tan(state.theta), state.q * sin(state.phi) - state.r * cos(state.phi)],
#             #               [-state.q * sin(state.phi) - state.r * cos(state.phi), 0]])
            
#             self.P = self.P + (Tp/self.N) * (A @ self.P + self.P @ A.T + self.Q)
       
#         psidot = y_F[6]
#         chidot = CTRL.gravity * tan(state.psi) * cos(state.chi - state.psi) / state.Vg
#         Vgdot = ((state.Va * cos(state.psi) + state.wn) * (state.Va * psidot * sin(state.psi)) + (state.Va * sin(state.psi) + state.we) * (state.Va * psidot * cos(state.psi)))/state.Vg
#         DVgdot_DPsi =  (-psidot * state.Va * (state.wn * cos(state.psi) + state.we * sin(state.psi)))/state.Vg
#         Dchidot_DVg = -CTRL.gravity * tan(state.phi) * cos(state.chi - state.psi) / (state.Vg**2)
#         Dchidot_Dchi = -CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi) / state.Vg
#         Dchidot_Dpsi = CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi)
        
#         self.A = np.array([[0., 0., cos(state.chi), -state.Vg * sin(state.chi), 0., 0., 0.],
#                         [0., 0., sin(state.chi), state.Vg * cos(state.chi), 0., 0., 0.],
#                         [0., 0., -Vgdot[0] / state.Vg, 0., -psidot[0] * state.Va * sin(state.psi), psidot[0] * state.Va * cos(state.psi), DVgdot_DPsi[0]],
#                         [0., 0., Dchidot_DVg, Dchidot_Dchi, 0., 0., Dchidot_Dpsi],
#                         [0., 0., 0., 0., 0., 0., 0.],
#                         [0., 0., 0., 0., 0., 0., 0.],
#                         [0., 0., 0., 0., 0., 0., 0.]])
        
#         self.P = self.P + (Tp/self.N) * (A @ self.P + self.P @ A.T + self.Q)

#         # Ad = np.eye(7) + self.A * Tp + self.A @ self.A * Tp**2
#         # self.P = self.P + (Tp/self.N) * (Ad @ self.P + self.P @ Ad.T + self.Q)
#         # print("P = ", self.P)
#         # convert to discrete time models
#         # self.P = 

#             # update P with discrete time model
#             # self.P = np.zeros((7,7))

#             # 
#             # Tout = self.Ts * 3
#             # # Tp = self.Ts
#             # Tp = Tout/self.N
#             # f = self.f(self.xhat, measurement, state)
#             # psidot = f[6]
#             # chidot = CTRL.gravity * tan(state.psi) * cos(state.chi - state.psi) / state.Vg
#             # Vgdot = ((state.Va * cos(state.psi) + state.wn) * (state.Va * psidot * sin(state.psi)) + (state.Va * sin(state.psi) + state.we) * (state.Va * psidot * cos(state.psi)))/state.Vg
            
#             # for i in range(0, self.N):
#             #     self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state).T
#             #     # A = jacobian(self.f, self.xhat, measurement, state)
#             #     DVgdot_DPsi =  (-psidot * state.Va * (state.wn * cos(state.psi) + state.we * sin(state.psi)))/state.Vg
#             #     Dchidot_DVg = -CTRL.gravity * tan(state.phi) * cos(state.chi - state.psi) / (state.Vg**2)
#             #     Dchidot_Dchi = -CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi) / state.Vg
#             #     Dchidot_Dpsi = CTRL.gravity * tan(state.phi) * sin(state.chi - state.psi)
                
#             #     self.A = np.array([[0., 0., cos(state.chi), -state.Vg * sin(state.chi), 0., 0., 0.],
#             #                   [0., 0., sin(state.chi), state.Vg * cos(state.chi), 0., 0., 0.],
#             #                   [0., 0., -Vgdot[0] / state.Vg, 0., -psidot[0] * state.Va * sin(state.psi), psidot[0] * state.Va * cos(state.psi), DVgdot_DPsi[0]],
#             #                   [0., 0., Dchidot_DVg, Dchidot_Dchi, 0., 0., Dchidot_Dpsi],
#             #                   [0., 0., 0., 0., 0., 0., 0.],
#             #                   [0., 0., 0., 0., 0., 0., 0.],
#             #                   [0., 0., 0., 0., 0., 0., 0.]])
#             #     A = self.A
#             #     # A = jacobian()
#             #     Ad = np.eye(7) + A * Tp + A @ A * Tp**2
#             #     # self.P = Ad @ self.P @ Ad.T + Tp**2 * self.Q
#             #     self.P = self.P + (Tp/self.N) * (A @ self.P + self.P @ A.T + self.Q)
                
#             #     print("P:", self.P)
#             #     if np.any(np.isnan(self.P)):
#             #         quit()


#     def measurement_update(self, measurement, state):
#         # always update based on wind triangle pseudo measurement
       
#         y_p = self.h_pseudo(self.xhat, measurement,state)
#         C = jacobian(self.h_pseudo, self.xhat, measurement, state)
#         y_m = np.array([[0, 0]]).T
#         S_inv = np.linalg.inv(self.R_pseudo + C @ self.P @ C.T)#np.zeros((2,2))
#         if True:
#             L = self.P @ C.T @ S_inv
#             LC = L @ C
#             tmp = np.eye(7) - LC
#             self.P = tmp @ self.P @ tmp.T + L @ self.R_pseudo @ L.T
            
#             print("5") 
#             self.xhat =  self.xhat + L @ (y_m-y_p)
            
#         # only update GPS when one of the signals changes
#         if (measurement.gps_n != self.gps_n_old) \
#             or (measurement.gps_e != self.gps_e_old) \
#             or (measurement.gps_Vg != self.gps_Vg_old) \
#             or (measurement.gps_course != self.gps_course_old):
            
#             # y_p = self.h_gps(self.xhat, measurement, state)
#             yhat = self.h_gps(self.xhat, measurement, state)
#             C = jacobian(self.h_gps, self.xhat, measurement, state)
#             y_chi = wrap(measurement.gps_course, yhat[3, 0])
#             y_m = np.array([[measurement.gps_n,
#                            measurement.gps_e,
#                            measurement.gps_Vg,
#                            y_chi]]).T
            
#             S_inv = np.linalg.inv(self.R_gps + C @ self.P @ C.T)#np.zeros((2,2))
       
#             if (y_m-yhat).T @ S_inv @ (y_m-yhat) < self.gps_threshold:
#             # if True:
#                 L = self.P @ C.T @ S_inv

#                 print("6")
#                 self.xhat = self.xhat + L @ (y_m - yhat)
#                 LC = L @ C
#                 tmp = np.eye(7) - LC
#                 self.P = tmp @ self.P @ tmp.T + L @ self.R_gps @ L.T

#             # update stored GPS signals
#             self.gps_n_old = measurement.gps_n
#             self.gps_e_old = measurement.gps_e
#             self.gps_Vg_old = measurement.gps_Vg
#             self.gps_course_old = measurement.gps_course
            


def jacobian(fun, x, measurement, state):
    # compute jacobian of fun with respect to x
    f = fun(x, measurement, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.0001  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, measurement, state)
        df = (f_eps - f) / eps
        
        J[:, i] = df[:, 0]
    return J