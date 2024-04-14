"""
compute_ss_model
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/4/2019 - RWB
"""
import sympy as sp
import numpy as np
from scipy.optimize import minimize
from tools.rotations import euler_to_rotation, quaternion_to_euler, euler_to_quaternion
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta


def compute_model(mav, trim_state, trim_input):
    # Note: this function alters the mav private variables
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

    # write transfer function gains to file
    file = open('model_coef.py', 'w')
    file.write('import numpy as np\n')
    file.write('x_trim = np.array([[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]]).T\n' %
               (trim_state.item(0), trim_state.item(1), trim_state.item(2), trim_state.item(3),
                trim_state.item(4), trim_state.item(5), trim_state.item(6), trim_state.item(7),
                trim_state.item(8), trim_state.item(9), trim_state.item(10), trim_state.item(11),
                trim_state.item(12)))
    file.write('u_trim = np.array([[%f, %f, %f, %f]]).T\n' %
               (trim_input.elevator, trim_input.aileron, trim_input.rudder, trim_input.throttle))
    file.write('Va_trim = %f\n' % Va_trim)
    file.write('alpha_trim = %f\n' % alpha_trim)
    file.write('theta_trim = %f\n' % theta_trim)
    file.write('a_phi1 = %f\n' % a_phi1)
    file.write('a_phi2 = %f\n' % a_phi2)
    file.write('a_theta1 = %f\n' % a_theta1)
    file.write('a_theta2 = %f\n' % a_theta2)
    file.write('a_theta3 = %f\n' % a_theta3)
    file.write('a_V1 = %f\n' % a_V1)
    file.write('a_V2 = %f\n' % a_V2)
    file.write('a_V3 = %f\n' % a_V3)
    file.write('A_lon = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lon[0][0], A_lon[0][1], A_lon[0][2], A_lon[0][3], A_lon[0][4],
     A_lon[1][0], A_lon[1][1], A_lon[1][2], A_lon[1][3], A_lon[1][4],
     A_lon[2][0], A_lon[2][1], A_lon[2][2], A_lon[2][3], A_lon[2][4],
     A_lon[3][0], A_lon[3][1], A_lon[3][2], A_lon[3][3], A_lon[3][4],
     A_lon[4][0], A_lon[4][1], A_lon[4][2], A_lon[4][3], A_lon[4][4]))
    file.write('B_lon = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lon[0][0], B_lon[0][1],
     B_lon[1][0], B_lon[1][1],
     B_lon[2][0], B_lon[2][1],
     B_lon[3][0], B_lon[3][1],
     B_lon[4][0], B_lon[4][1],))
    file.write('A_lat = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lat[0][0], A_lat[0][1], A_lat[0][2], A_lat[0][3], A_lat[0][4],
     A_lat[1][0], A_lat[1][1], A_lat[1][2], A_lat[1][3], A_lat[1][4],
     A_lat[2][0], A_lat[2][1], A_lat[2][2], A_lat[2][3], A_lat[2][4],
     A_lat[3][0], A_lat[3][1], A_lat[3][2], A_lat[3][3], A_lat[3][4],
     A_lat[4][0], A_lat[4][1], A_lat[4][2], A_lat[4][3], A_lat[4][4]))
    file.write('B_lat = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lat[0][0], B_lat[0][1],
     B_lat[1][0], B_lat[1][1],
     B_lat[2][0], B_lat[2][1],
     B_lat[3][0], B_lat[3][1],
     B_lat[4][0], B_lat[4][1],))
    file.write('Ts = %f\n' % Ts)
    file.close()


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    mav._state = trim_state
    mav._update_velocity_data()
    Va_trim = mav._Va
    alpha_trim = mav._alpha
    delta_e_trim = trim_state[0]
    throttle = trim_state[3]
    # print("delta_e", delta_e_trim)

    phi, theta_trim, psi = quaternion_to_euler(trim_state[6:10])

    S = MAV.S_wing
    b = MAV.b
    rho = MAV.rho
    c = MAV.c
    Jy = MAV.Jy

    ###### TODO ######
    # define transfer function constants
    a_phi1 = -.5 * rho * Va_trim**2 * S*MAV.C_p_p * b / (2*Va_trim)
    a_phi2 = .5 * rho * Va_trim**2 * S * b * MAV.C_p_delta_a
    a_theta1 = - rho * Va_trim**2 * c * S * MAV.C_m_q * c /(2 * 2 * Va_trim)
    a_theta2 = -rho * Va_trim**2 * c * S  * MAV.C_m_alpha / (2 * Jy)
    a_theta3 = rho * Va_trim**2 * c * S * MAV.C_m_delta_e / (2 * Jy)


    # Compute transfer function coefficients using new propulsion model
    a_V1 = rho * Va_trim * S / MAV.mass * (MAV.C_D_0 + MAV.C_D_alpha * alpha_trim * MAV.C_D_delta_e * delta_e_trim) - dT_dVa(mav,Va_trim,throttle) / MAV.mass
    a_V2 = 1/MAV.mass * dT_ddelta_t(mav,Va_trim, throttle)
    a_V3 = MAV.gravity * np.cos(theta_trim - alpha_trim)

    return Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, a_V1, a_V2, a_V3


def compute_ss_model(mav, trim_state, trim_input):
    x_euler = euler_state(trim_state)
    
    ##### TODO #####
    A = df_dx(mav, x_euler, trim_input)
    B = df_du(mav, x_euler, trim_input) #12x4
    # print(A)

    # extract longitudinal states (u, w, q, theta, pd)
    A_lon = np.zeros((5,5))
    B_lon = np.zeros((5,2))
    E1 = np.array([[0,0,0,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,1,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,1,0],
                      [0,0,0,0,0,0,0,1,0,0,0,0],
                      [0,0,-1,0,0,0,0,0,0,0,0,0]]) 
    A_lon = E1 @ A @ E1.T
    E2 = np.array([[1,0,0,0],
                    [0,0,0,1]])
    B_lon = E1 @ B @ E2.T

    A_lon_eig = np.linalg.eigvals(A_lon)
    # print("A_Lon_eig:",A_lon_eig)


    # change pd to h
    # extract lateral states (v, p, r, phi, psi)
    A_lat = np.zeros((5,5))
    B_lat = np.zeros((5,2))
    E3 = np.array([[0,0,0,0,1,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,1,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,1],
                      [0,0,0,0,0,0,1,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,1,0,0,0]]) 
    E4 = np.array([[0,1,0,0],
                   [0,0,1,0]])
    A_lat = E3 @ A @ E3.T
    B_lat = E3 @ B @ E4.T 
    A_lat_eig = np.linalg.eigvals(A_lat)
    # print("A_Lat_eig:",A_lat_eig)

    return A_lon, B_lon, A_lat, B_lat

def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles
    
    ##### TODO #####
    euler = quaternion_to_euler(x_quat[6:10])
    x_euler = np.zeros((12,1))
    x_euler[0:6] = np.copy(x_quat)[0:6]
    x_euler[6:9,0] = euler
    x_euler[9:12] = np.copy(x_quat[10:13])
    return x_euler

def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions

    ##### TODO #####
    quat = euler_to_quaternion(x_euler[6], x_euler[7], x_euler[8])[:,0,0]
    x_quat = np.zeros((13,1))
    x_quat[0:6] = np.copy(x_euler)[0:6]
    x_quat[6:10,0] = quat
    x_quat[10:13] = np.copy(x_euler[9:12])

    return x_quat

def f_euler(mav, x_euler, delta):
    # return 12x1 dynamics (as if state were Euler state)
    # compute f at euler_state, f_euler will be f, except for the attitude states

    # need to correct attitude states by multiplying f by
    # partial of Quaternion2Euler(quat) with respect to quat
    # compute partial Quaternion2Euler(quat) with respect to quat
    # dEuler/dt = dEuler/dquat * dquat/dt
    x_quat = quaternion_state(x_euler)
    mav._state = x_quat
    mav._update_velocity_data()
    ##### TODO #####

    f_euler_ = np.zeros((12,1))
    forces_moments = mav._forces_moments(delta)
    f = mav._f(x_quat,forces_moments)
    # print(len(x_euler))
    # print(len(x_quat))
    dE_dq = x_euler
    # f_euler_ = df_dx(mav, x_quat * mav._f(x_quat, forces_moments), delta)
    De_Dxq = np.zeros([12,13])
    for i in range(0,6):
        De_Dxq[i, i] = 1
    for i in range(9,12):
        De_Dxq[i,i+1] = 1
    De_Dxq[6:9,6:10] = dtheta_dq(mav, x_quat[6:10])

    f_euler_ =  De_Dxq @ f

    return f_euler_

def dtheta_dq(mav, quat):
    quat = quat.T
    e0, e1, e2, e3 = sp.symbols('e0 e1 e2 e3')

    # Define the functions phi, theta, psi
    phi = sp.atan2(2.0 * (e0 * e1 + e2 * e3), e0**2.0 + e3**2.0 - e1**2.0 - e2**2.0)
    theta = sp.asin(2.0 * (e0 * e2 - e1 * e3))
    psi = sp.atan2(2.0 * (e0 * e3 + e1 * e2), e0**2.0 + e1**2.0 - e2**2.0 - e3**2.0)

    derivatives_phi = [sp.diff(phi, e0), sp.diff(phi, e1), sp.diff(phi, e2), sp.diff(phi, e3)]
    derivatives_theta = [sp.diff(theta, e0), sp.diff(theta, e1), sp.diff(theta, e2), sp.diff(theta, e3)]
    derivatives_psi = [sp.diff(psi, e0), sp.diff(psi, e1), sp.diff(psi, e2), sp.diff(psi, e3)]
    deu_dq = [derivatives_phi, derivatives_theta, derivatives_psi]
    e_values = {e0: quat[0,0], e1: quat[0,1], e2: quat[0,2], e3: quat[0,3]}
    
    deu_dq = np.array([[derivative.subs(e_values).evalf() for derivative in derivative_list] for derivative_list in deu_dq], dtype=float)
    # print(deu_dq)
    
    return deu_dq


def df_dx(mav, x_euler, delta):
    # take partial of f_euler with respect to x_euler
    eps = 0.01  # deviation

    ##### TODO #####
    A = np.zeros((12, 12))  # Jacobian of f wrt x
    x_quat = quaternion_state(x_euler)
    forces_moments = mav._forces_moments(delta)
    f_at_x = f_euler(mav, x_euler, delta)
    for i in range(0,12):
        x_eps = np.copy(x_euler)
        x_eps[i][0] += eps
        f_at_x_eps = f_euler(mav, x_eps, delta)
        df_dxi = (f_at_x_eps - f_at_x) / eps
        A[:,i] = df_dxi[:,0]
    return A


def df_du(mav, x_euler, delta):
    # take partial of f_euler with respect to input
    eps = 0.001  # deviation
    delta_eps = MsgDelta()
    delta_num = np.array([delta.elevator, delta.aileron, delta.rudder, delta.throttle])
    # print(delta_num)
    ##### TODO #####
    B = np.zeros((12, 4))  # Jacobian of f wrt u
    f_at_x = f_euler(mav, x_euler, delta)
    for i in range(0,4):
        x_eps = np.copy(delta_num)
        x_eps[i] += eps
        delta_eps.elevator, delta_eps.aileron, delta_eps.rudder, delta_eps.throttle = x_eps
        f_at_x_eps = f_euler(mav, x_euler, delta_eps)
        df_dxi = (f_at_x_eps - f_at_x) / eps
        B[:,i] = df_dxi[:,0]

    return B


def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    eps = 0.01

    ##### TODO #####
    dT_dVa = 0
    f_at_x = mav._motor_thrust_torque(Va, delta_t)[0]
    x_eps = Va + eps
    f_at_x_eps = mav._motor_thrust_torque(x_eps, delta_t)[0]
    dT_dVa = (f_at_x_eps - f_at_x) / eps
    return dT_dVa


def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    eps = 0.01

    ##### TODO #####
    dT_ddelta_t = 0
    f_at_x = mav._motor_thrust_torque(Va, delta_t)[0]
    x_eps = delta_t + eps
    f_at_x_eps = mav._motor_thrust_torque(Va, x_eps)[0]
    dT_ddelta_t= (f_at_x_eps - f_at_x) / eps
    return dT_ddelta_t
