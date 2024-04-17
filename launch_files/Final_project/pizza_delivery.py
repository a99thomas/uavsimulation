"""
mavsim_python
    - Chapter 11 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        3/26/2019 - RWB
        2/27/2020 - RWB
        1/5/2023 - David L. Christiansen
        7/13/2023 - RWB
"""
import os, sys
# insert parent directory at beginning of python search path
from pathlib import Path
sys.path.insert(0,os.fspath(Path(__file__).parents[2]))
# use QuitListener for Linux or PC <- doesn't work on Mac
#from tools.quit_listener import QuitListener
import numpy as np
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN
from models.mav_dynamics_sensors import MavDynamics
from models.box_dynamics import BoxDynamics
from message_types.msg_state import MsgState
from models.wind_simulation import WindSimulation
from controllers.autopilot import Autopilot
from estimators.observer import Observer
from planners.path_follower import PathFollower
# from chap11.path_manager_cycle import PathManager
from planners.delivery_manager import DeliveryManager
from planners.delivery_location import Delivery_Zone
from viewers.view_box_manager import ViewManager

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
box = BoxDynamics(SIM.ts_simulation)
box_state = MsgState()
autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)
path_follower = PathFollower()
path_manager = DeliveryManager()
viewers = ViewManager(animation=True, data=False, path = False, waypoint=True)
#quitter = QuitListener()

# waypoint definition
from message_types.msg_waypoints import MsgWaypoints
waypoints = MsgWaypoints()
# waypoints.type = 'straight_line'
waypoints.type = 'fillet'
# waypoints.type = 'dubins'
Va = PLAN.Va0
waypoints.add(np.array([[0, 0, -100]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1000, 500, -100]]).T, Va, np.radians(45), np.inf, 0, 0)


waypoints.add(np.array([[1000, -600, -100]]).T, Va, np.radians(45), np.inf, 0, 0)

# waypoints.add(np.array([[1000, -200, -250]]).T, Va, np.radians(45), np.inf, 0, 0)
# waypoints.add(np.array([[1000, 100, -100]]).T, Va, np.radians(45), np.inf, 0, 0)
# waypoints.add(np.array([[1000, 1000, -100]]).T, Va, np.radians(-135), np.inf, 0, 0)
# print(waypoints.ned[0:3,-1])
delivery_destination = np.array([1000, 300, -100])
waypoints.add(np.reshape(delivery_destination,(3,1)), Va, np.radians(45), np.inf, 0, 0)

waypoints.add(np.array([[0, 600, -10]]).T, Va, np.radians(0), np.inf, 0, 0)


delivery_zone = Delivery_Zone(delivery_destination)

# initialize the simulation time
sim_time = SIM.start_time
end_time = 400

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    estimated_state = observer.update(measurements)  # estimate states from measurements
    estimated_state = mav.true_state  # uses true states in the control
    # box_state.north = mav.true_state.north
    # box_state.east = mav.true_state.east
    # box_state.altitude = mav.true_state.altitude - 20


    # -------path manager-------------
    path = path_manager.update(waypoints, PLAN.R_min, estimated_state)

    # -------path follower-------------
    autopilot_commands = path_follower.update(path, estimated_state)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics
    box.update(box_state, estimated_state, waypoints, sim_time)

    # -------update viewer-------------
    viewers.update(
        sim_time,
        true_state=mav.true_state,  # true states
        estimated_state=estimated_state,  # estimated states
        commanded_state=commanded_state,  # commanded states
        delta=delta,  # inputs to aircraft
        path=path, # path
        delivery_zone = delivery_zone,
        waypoints=waypoints, # waypoints
        box_state=box.true_state
    )

    # -------Check to Quit the Loop-------
    # if quitter.check_quit():
    #     break

    # -------increment time-------------
    sim_time += SIM.ts_simulation
    # print(mav.true_state.north)

# close viewers
viewers.close(dataplot_name="ch11_data_plot")



