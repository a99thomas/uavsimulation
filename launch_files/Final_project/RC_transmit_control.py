"""
mavsim_python
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        2/2/2019 - RWB
        1/5/2023 - David L. Christiansen
        7/13/2023 - RWB
"""
import os, sys
# insert parent directory at beginning of python search path
from pathlib import Path
sys.path.insert(0,os.fspath(Path(__file__).parents[2]))
# use QuitListener for Linux or PC <- doesn't work on Mac
#from tools.quit_listener import QuitListener
import pyqtgraph as pg
import numpy as np
import parameters.simulation_parameters as SIM
from viewers.fpv_viewer import MavViewer
from viewers.data_viewer import DataViewer
from message_types.msg_state import MsgState
from models.mav_dynamics_control import MavDynamics
from models.box_dynamics import BoxDynamics
from models.wind_simulation import WindSimulation
from models.trim import compute_trim
from models.compute_models import compute_model
from tools.signals import Signals
from controllers.joystick_controller import Joystick


#quitter = QuitListener()

VIDEO = False
PLOTS = False #True
ANIMATION = True  #True
SAVE_PLOT_IMAGE = False
COMPUTE_MODEL = True

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap5_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_videoslow)

#initialize the visualization
if ANIMATION or PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    mav_view = MavViewer(app=app)  # initialize the mav viewer
if PLOTS:
    # initialize view of data plots
    data_view = DataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simslow)
mav = MavDynamics(SIM.ts_simslow)
box = BoxDynamics(SIM.ts_simulation)
joystick = Joystick()
box_state = MsgState()

# use compute_trim function to compute trim state and trim input
Va = 25.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(mav, Va, gamma)
mav._state = trim_state  # set the initial state of the mav to the trim state
delta = trim_input  # set input to constant constant trim input

# # compute the state space model linearized about trim
if COMPUTE_MODEL:
    compute_model(mav, trim_state, trim_input)

# this signal will be used to excite modes
input_signal = Signals(amplitude=0.3,
                       duration=0.3,
                       start_time=5.0)
delta_e_trim = delta.elevator
delta_a_trim = delta.aileron
delta_r_trim = delta.rudder

# initialize the simulation time
sim_time = SIM.start_time
end_time = 120

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < end_time:

    # -------physical system-------------
    #current_wind = wind.update()  # get the new wind vector
    current_wind = np.zeros((6, 1))
    delta.throttle, delta.rudder, delta.aileron, delta.elevator, trigger = joystick.get_transmitter_values()
    # delta.rudder = delta.rudder #+  trim_input.rudder
    # delta.aileron = delta.aileron #+ trim_input.aileron
    # delta.elevator = delta.elevator #+ trim_input.elevator
    

    mav.update(delta, current_wind)  # propagate the MAV dynamics
    box.update2(box_state, mav.true_state, sim_time, trigger=trigger)

    # -------update viewer-------------
    if ANIMATION:
        mav_view.update(mav.true_state, box.true_state)  # plot body of MAV
    if PLOTS:
        plot_time = sim_time
        data_view.update(mav.true_state,  # true states
                            None,  # estimated states
                            None,  # commanded states
                            delta)  # inputs to aircraft
    if ANIMATION or PLOTS:
        app.processEvents()
    if VIDEO is True:
        video.update(sim_time)
        
    # -------Check to Quit the Loop-------
    # if quitter.check_quit():
    #     break

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if SAVE_PLOT_IMAGE:
    data_view.save_plot_image("ch5_plot")

if VIDEO is True:
    video.close()



