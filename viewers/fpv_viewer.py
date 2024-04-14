"""
Mav Viewer, Aaron Thomas 1/16/24
"""
import pyqtgraph.opengl as gl
import pyqtgraph.Vector as Vector
from viewers.draw_mav import DrawMav
from viewers.draw_box import DrawBox
from time import time

class MavViewer():
    def __init__(self, app, ts_refresh = 1/30.):
        self.scale = 100
        self.fast = 0
        # initialize Qt gui application and window
        self.app = app  # initialize QT, external so that only one QT process is running
        self.window = gl.GLViewWidget()  # initialize the view object
        self.window.setWindowTitle('Spacecraft Viewer')
        grid = gl.GLGridItem() # make a grid to represent the ground
        grid.scale(20, 20, 20) # set the size of the grid (distance between each line)
        self.window.addItem(grid) # add grid to viewer
        # self.window.setCameraPosition(distance=200) # distance from center of plot to camera
        self.window.setBackgroundColor('k')  # set background color to black
        self.window.setGeometry(0, 0, 1500, 1000)  # args: upper_left_x, upper_right_y, width, height

        ##For FPV Mode
        center = self.window.cameraPosition()
        center.setX(250)
        center.setY(500)
        center.setZ(0)
        self.window.setCameraPosition(pos=center, distance=200, elevation=0, azimuth=-90)
        self.window.show()  # display configured window
        self.window.raise_() # bring window to the front
        # self.window.set


        self.plot_initialized = False # has the mav been plotted yet?
        self.sc_plot = []
        self.ts_refresh = ts_refresh
        self.t = time()
        self.t_next = self.t

    def update(self, state, boxstate):
        # initialize the drawing the first time update() is called
        center = self.window.cameraPosition()
        center.setX(250)
        center.setY(500)
        center.setZ(0)
        self.window.setCameraPosition(pos=center, distance=250, elevation=0, azimuth=(-90 - state.chi*180./3.14))
        if not self.plot_initialized:
            self.sc_plot = DrawMav(state, self.window)
            self.box_plot = DrawBox(boxstate, self.window)
            self.plot_initialized = True
        # else update drawing on all other calls to update()
        else:
            if self.fast == 1:
                t = time()
                if t-self.t_next > 0.0:
                    self.sc_plot.update(state)
                    self.box_plot.update(boxstate)
                    self.t = t
                    self.t_next = t + self.ts_refresh
            else:
                self.sc_plot.update(state)
                self.box_plot.update(boxstate)
        # update the center of the camera view to the spacecraft location
        view_location = Vector(state.east, state.north, state.altitude)  # defined in ENU coordinates
        self.window.opts['center'] = view_location
        # 
        # redraw

    def process_app(self):
        self.app.processEvents()

    def clear_viewer(self):
        self.window.clear()