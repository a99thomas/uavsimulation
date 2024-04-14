"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
# import models.mav_dynamics as mav
import numpy as np
if __name__ == "__main__":
    import os, sys
    from pathlib import Path
    sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
from tools.transfer_function import TransferFunction
from parameters import aerosonde_parameters as param
# import models.mav_dynamics_control as MAV
# from models.mav_dynamics_control import MavDynamics as MAV

class WindSimulation:
    def __init__(self, Ts, gust_flag = True, steady_state = np.array([[0., 0., 0.]]).T):
        # steady state wind defined in the inertial frame
        self._steady_state = steady_state
        # mav_states = MAV(Ts)  # Provide the appropriate value for Ts
        # Va = 25# mav_states._Va
        Va = param.Va0
        # self.Va = 25
        # Va = np.linalg.norm(mav.MavDynamics.true_state.Va)
        
        ##### TODO #####
        #   Dryden gust model parameters (pg 56 UAV book)
        alt, Lu, Lw, sigma_u, sigma_w = 50, 200, 50, 1.06, 0.7 #Low alt, low turb
        # alt, Lu, Lw, sigma_u, sigma_w = 50, 200, 50, 2.12, 1.4 #Low alt, med turb
        # alt, Lu, Lw, sigma_u, sigma_w = 50, 200, 50, 1.5, 1.5 #Med alt, low turb
        # alt, Lu, Lw, sigma_u, sigma_w = 50, 200, 50, 3.0, 3.0 #Med alt, med turb
        Lv = Lu
        sigma_v = sigma_u

        # Dryden transfer functions (section 4.4 UAV book) - Fill in proper num and den
        # self.u_w = TransferFunction(num=np.array([[0]]), den=np.array([[1,1]]),Ts=Ts)
        self.u_w = TransferFunction(num=np.array([[sigma_u * np.sqrt(2*Va/(np.pi * Lu))]]), den=np.array([[1,Va/Lu]]),Ts=Ts) 
        
        num_coeff = sigma_v * np.sqrt(3 * Va / (np.pi * Lv))
        # print("UW", self.u_w)
        self.v_w = TransferFunction(num=np.array([[num_coeff, num_coeff * Va / (np.sqrt(3) * Lv)]]), den=np.array([[1, 2*Va / Lv, (Va / Lv)**2]]),Ts=Ts)
        
        
        num_coeff = sigma_w * np.sqrt(3 * Va / (np.pi * Lw))
        self.w_w = TransferFunction(num=np.array([[num_coeff, Va / (np.sqrt(3) * Lw)]]), den=np.array([[1, 2 * Va / Lw, (Va / Lw)**2]]),Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        gust = np.array([[self.u_w.update(np.random.randn())],
                         [self.v_w.update(np.random.randn())],
                         [self.w_w.update(np.random.randn())]])
        # print("Gust", gust)
        # gust = np.array([[0],[0],[0]])
        return np.concatenate(( self._steady_state, gust ))


if __name__ == "__main__":
    ws = WindSimulation(.01)
    print(ws.update())
    pass