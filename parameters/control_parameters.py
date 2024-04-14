import numpy as np
# import launch_files.chap05.model_coef as TF
import parameters.aerosonde_parameters as MAV


#### TODO #####
gravity = MAV.gravity  # gravity constant
Va0 = 25 #TF.Va_trim
rho = 1.2682 # density of air
sigma = 0  # low pass filter gain for derivative

#----------roll loop-------------
# get transfer function data for delta_a to phi
wn_roll = 0
zeta_roll = 0
roll_kp = .5
roll_kd = 1.

#----------course loop-------------
wn_course = 0
zeta_course = 0
course_kp = 5.
course_ki = 1.

#----------yaw damper-------------
yaw_damper_p_wo = .5
yaw_damper_kr = .2

#----------pitch loop-------------
wn_pitch = 0
zeta_pitch = 0 
pitch_kp = -0.6
pitch_kd = -.1
K_theta_DC = 0.

#----------altitude loop-------------
wn_altitude = 0
zeta_altitude = 0
altitude_kp = .1
altitude_ki = 0.01
altitude_zone = 75.

#---------airspeed hold using throttle---------------
wn_airspeed_throttle = 0
zeta_airspeed_throttle = 0
airspeed_throttle_kp = .5
airspeed_throttle_ki = .5