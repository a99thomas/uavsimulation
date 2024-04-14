"""
mavsim_python: drawing tools
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - RWB
        3/30/2022 - RWB
        7/13/2023 - RWB
        3/25/2024 - RWB
"""

import numpy as np
from planners.dubins_parameters import DubinsParameters
from message_types.msg_state import MsgState
from message_types.msg_path import MsgPath
from message_types.msg_waypoints import MsgWaypoints
from tools.rotations import euler_to_rotation
import time


class PathManager:
    '''
        Path manager

        Attributes
        ----------
        path : MsgPath
            path message sent to path follower
        num_waypoints : int
            number of waypoints
        ptr_previous : int
            pointer to previous waypoint
            MAV is traveling from previous to current waypoint
        ptr_current : int
            pointer to current waypoint
        ptr_next : int
            pointer to next waypoint
        halfspace_n : np.nparray (3x1)
            the normal vector that defines the current halfspace plane
        halfspace_r : np.nparray (3x1)
            the inertial vector that defines a point on the current halfspace plane
        manager_state : int
            state of the manager state machine
        manager_requests_waypoints : bool
            a flag for handshaking with the path planner
            True when new waypoints are needed, i.e., at the end of waypoint list.
        dubins_path : DubinsParameters
            A class that defines a dubins path      

        Methods
        -------
        update(waypoints, radius, state)

        _initialize_pointers() :
            initialize the points to 0(previous), 1(current), 2(next)  
        _increment_pointers() :  
            add one to every pointer - currently does it modulo num_waypoints          
        _inHalfSpace(pos):
            checks to see if the position pos is in the halfspace define

        _line_manager(waypoints, state):
            Assumes straight-line paths.  Transition is from one line to the next
            _construct_line(waypoints): 
                used by line manager to construct the next line path

        _fillet_manager(waypoints, radius, state):
            Assumes straight-line waypoints.  Constructs a fillet turn between lines.
            _construct_fillet_line(waypoints, radius):
                used by _fillet_manager to construct the next line path
            _construct_fillet_circle(waypoints, radius):
                used by _fillet_manager to construct the fillet orbit
            
        _dubins_manager(waypoints, radius, state):
            Assumes dubins waypoints.  Constructs Dubin's path between waypoints
            _construct_dubins_circle_start(waypoints, dubins_path):
                used by _dubins_manager to construct the start orbit
            _construct_dubins_line(waypoints, dubins_path):
                used by _dubins_manager to construct the middle line
            _construct_dubins_circle_end(waypoints, dubins_path):
                used by _dubins_manager to construct the end orbit
    '''
    def __init__(self):
        self._path = MsgPath()
        self._num_waypoints = 0
        self._ptr_previous = 0
        self._ptr_current = 0
        self._ptr_next = 1
        self._halfspace_n = np.inf * np.ones((3,1))
        self._halfspace_r = np.inf * np.ones((3,1))
        self._manager_state = 1
        self.manager_requests_waypoints = True
        self.dubins_path = DubinsParameters()
        self.rewind = False


    def update(self, 
               waypoints: MsgWaypoints, 
               radius: float, 
               state: MsgState) -> MsgPath:
        if waypoints.num_waypoints == 0:
            self.manager_requests_waypoints = True
        if self.manager_requests_waypoints is True \
                and waypoints.flag_waypoints_changed is True:
            self.manager_requests_waypoints = False
        if waypoints.type == 'straight_line':
            self.line_manager(waypoints, state)
        elif waypoints.type == 'fillet':
            self._fillet_manager(waypoints, radius, state)
        elif waypoints.type == 'dubins':
            self._dubins_manager(waypoints, radius, state)
        else:
            print('Error in Path Manager: Undefined waypoint type.')
        return self._path


            
    def line_manager(self, waypoints: MsgWaypoints, state):
        

        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # i = self._ptr_current
        if self._inHalfSpace(mav_pos):
            if np.linalg.norm(mav_pos - waypoints.ned[:,self._ptr_current:self._ptr_current+1]) < 100:
                if self._ptr_next < waypoints.num_waypoints:
                    waypoints.flag_waypoints_changed = True
                else:
                    self.flag_manager_requests_waypoints = True
        
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed:
            waypoints.flag_waypoints_changed = False
            self.flag_manager_requests_waypoints = False
            self._num_waypoints = waypoints.num_waypoints
            self._increment_pointers()
            self.construct_line(waypoints)


    def _fillet_manager(self,  
                        waypoints: MsgWaypoints, 
                        radius: float, 
                        state: MsgState):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if self._manager_state == 1:
            if self._inHalfSpace(mav_pos):
                # if np.linalg.norm(mav_pos - waypoints.ned[:,self._ptr_current:self._ptr_current+1]) < radius+50:
                if self._ptr_next < waypoints.num_waypoints:
                    self._manager_state = 2
                    self._construct_fillet_circle(waypoints,radius)

        elif self._manager_state == 2:
            if self._inHalfSpace(mav_pos):
                if self._ptr_next < waypoints.num_waypoints:
                    waypoints.flag_waypoints_changed = True
                else:
                    self.flag_manager_requests_waypoints = True
        
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed:
            waypoints.flag_waypoints_changed = False
            self.flag_manager_requests_waypoints = False
            self._manager_state = 1
            self._num_waypoints = waypoints.num_waypoints
            self._increment_pointers()
            self._construct_fillet_line(waypoints, radius)
      

    def _dubins_manager(self,  
                        waypoints: MsgWaypoints, 
                        radius: float, 
                        state: MsgState):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        
        # if waypoints.flag_waypoints_changed == True:

        #     self._initialize_pointers()
        #     self._manager_state = 1
        # print(self._manager_state)
        
        if waypoints.flag_waypoints_changed:
            self.flag_manager_requests_waypoints = False
            waypoints.flag_waypoints_changed = False
            self._num_waypoints = waypoints.num_waypoints
            self._initialize_pointers()
            ps = waypoints.ned[:, self._ptr_previous:self._ptr_previous+1]
            pe = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
            pnext = waypoints.ned[:, self._ptr_next:self._ptr_next+1]
            # chie = pnext - pe
            # chis = state.chi
            chis = waypoints.course.item(self._ptr_previous)
            chie = waypoints.course.item(self._ptr_current)
            # np.arctan2(pnext.item(1) - pe.item(1), pnext.item(0) - pe.item(0))
            self._num_waypoints = waypoints.num_waypoints
            self.dubins_path.update(ps, chis, pe, chie, radius)
            self._halfspace_n = -self.dubins_path.n1
            self._halfspace_r = self.dubins_path.r1

            self._construct_dubins_circle_start(waypoints, self.dubins_path)

        if self._manager_state == 1: #Circle part 1
            # waypoints.flag_waypoints_changed = False
            if self._inHalfSpace(mav_pos):
                self._halfspace_n = self.dubins_path.n1
                self._manager_state = 2

        elif self._manager_state == 2: #Circle part 2
            if self._inHalfSpace(mav_pos):
                self._construct_dubins_line(waypoints, self.dubins_path)
                self._manager_state = 3

        elif self._manager_state == 3: #Straight
            
            if self._inHalfSpace(mav_pos):
                # self._halfspace_n = self.dubins_path.n2
                self._construct_dubins_circle_end(waypoints, self.dubins_path)
                self._manager_state = 4

        elif self._manager_state == 4: #Circle 2 part 1
            if self._inHalfSpace(mav_pos):
                self._halfspace_n = self.dubins_path.n3
                self._manager_state = 5

        elif self._manager_state == 5: # Circle 2 part 2/Circle 1 part 1
            if self._inHalfSpace(mav_pos):
                self._increment_pointers()
                pprev = waypoints.ned[:, self._ptr_previous:self._ptr_previous+1]
                pe = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
                pnext = waypoints.ned[:, self._ptr_next:self._ptr_next+1]
                # chie = pnext - pe
                # chis = np.arctan2(pe.item(1) - pprev.item(1), pe.item(0) - pprev.item(0))
                # chie = np.arctan2(pnext.item(1) - pe.item(1), pnext.item(0) - pe.item(0))
                chis = waypoints.course.item(self._ptr_previous)
                chie = waypoints.course.item(self._ptr_current)
                self._num_waypoints = waypoints.num_waypoints
                self.dubins_path.update(pprev, chis, pe, chie, radius)
                # self._construct_dubins_circle_start(waypoints, self.dubins_path)
                self._construct_dubins_circle_start(waypoints, self.dubins_path)
                self._manager_state = 1

        
        # if waypoints.flag_waypoints_changed:
        #     waypoints.flag_waypoints_changed = False
        #     self.flag_manager_requests_waypoints = False
        ##### TODO #####
        # Use functions - self._initialize_pointers(), self._dubins_path.update(),
        # self._construct_dubins_circle_start(), self._construct_dubins_line(),
        # self._inHalfSpace(), self._construct_dubins_circle_end(), self._increment_pointers(),

        # Use variables - self._num_waypoints, self._dubins_path, self._ptr_current,
        # self._ptr_previous, self._manager_state, self.manager_requests_waypoints,
        # waypoints.__, radius


    def _initialize_pointers(self):
        
        if self._num_waypoints >= 3:
            self._ptr_previous = 0
            self._ptr_current = 1
            self._ptr_next = 2
        #     self._ptr_previous = self._ptr_current
        #     self._ptr_current = self._ptr_next
        #     if self._ptr_current+1 < self._num_waypoints:
        #         self._ptr_next = self._ptr_next + 1
        #     else:
        #         self._ptr_next = 0
        #         self.rewind = True
        else:
            print('Error Path Manager: need at least three waypoints')
            print(self._num_waypoints)

    def _increment_pointers(self):
        ##### TODO #####
        # self._ptr_previous = 0
        # self._ptr_current = 0
        # self._ptr_next = 0
        if self._num_waypoints >= 3:
            self._ptr_previous = self._ptr_current
            self._ptr_current = self._ptr_next
            if self._ptr_current+1 < self._num_waypoints:
                self._ptr_next = self._ptr_next + 1
            else:
                self._ptr_next = 0
                self.rewind = True
        else:
            print('Error Path Manager: need at least three waypoints')
            print(self._num_waypoints)

    def construct_line(self, waypoints):
        previous = waypoints.ned[:, self._ptr_previous:self._ptr_previous+1]
        current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        next = waypoints.ned[:, self._ptr_next:self._ptr_next+1]
        if self._ptr_current == 9999:
            current = previous + 100*self._path.line_direction
        else:
            current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        if self._ptr_next == 9999:
            next = waypoints.ned[:, self._ptr_current:self._ptr_next]
        self._path.plot_updated = False
        self._path.type = 'line'  
        self._path.airspeed = waypoints.airspeed.item(self._ptr_current), 
        self._path.line_origin = previous
        self._path.line_direction = current - previous
        qi = (next - current) / np.linalg.norm(next-current)
        ni = (previous + current) / np.linalg.norm(previous + current)

        self._halfspace_r = current
        self._halfspace_n = ni
        print(self._halfspace_n)
        print(self._halfspace_r)
        

    def _construct_fillet_line(self, 
                               waypoints: MsgWaypoints, 
                               radius: float):
        previous = waypoints.ned[:, self._ptr_previous:self._ptr_previous+1]
        current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        next = waypoints.ned[:, self._ptr_next:self._ptr_next+1]
        if self._ptr_current == 9999:
            current = previous + 100*self._path.line_direction
        else:
            current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        if self._ptr_next == 9999:
            next = waypoints.ned[:, self._ptr_current:self._ptr_next]
        self._path.plot_updated = False
        self._path.type = 'line'  
        self._path.airspeed = waypoints.airspeed.item(self._ptr_current), 
        self._path.line_origin = previous
        self._path.line_direction = current - previous
        qi = (next - current) / np.linalg.norm(next-current)
        qi_1 = (current - previous) / np.linalg.norm(current - previous)
        ni = (previous + current) / np.linalg.norm(previous + current)
        rho = np.arccos(-qi_1.T @ qi)
        z = current - (radius / np.tan(rho / 2)) * qi_1
        self._halfspace_r = z
        self._halfspace_n = qi_1
        print(self._halfspace_n)
        print(self._halfspace_r)

    def _construct_fillet_circle(self, 
                                 waypoints: MsgWaypoints, 
                                 radius: float):
        previous = waypoints.ned[:, self._ptr_previous:self._ptr_previous+1]
        current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        next = waypoints.ned[:, self._ptr_next:self._ptr_next+1]
        # print("Current", self._ptr_current)
        ##### TODO #####
        if self._ptr_current == 9999:
            current = previous + 100*self._path.line_direction
        else:
            current = waypoints.ned[:, self._ptr_current:self._ptr_current+1]
        if self._ptr_next == 9999:
            next = waypoints.ned[:, self._ptr_current:self._ptr_next]
    
        qi = (next - current) / np.linalg.norm(next-current)
        qi_1 = (current - previous) / np.linalg.norm(current - previous)
        rho = np.arccos(-qi_1.T @ qi).item(0)
        z = current + (radius / np.tan(rho / 2)) * qi

        self._path.plot_updated = False
        self._path.type = 'orbit'  
        self._path.airspeed = waypoints.airspeed.item(self._ptr_current)
        self._path.orbit_center = current - (radius / np.sin(rho/2))*((qi_1-qi)/np.linalg.norm(qi_1-qi))
        direction = np.sign(qi_1[0] * qi[1] - qi_1[1] * qi[0])
        self._path.orbit_radius = radius

        if direction == 1:
            self._path.orbit_direction = 'CW'
        else:
            self._path.orbit_direction = 'CCW'
        self._halfspace_r = z
        self._halfspace_n = qi
        ##### TODO #####
        # current = ?
        # next = ?

        # update halfspace variables
        # self._halfspace_n =
        # self._halfspace_r = 
        
        # Update path variables
        # self._path.__ =

    def _construct_dubins_circle_start(self, 
                                       waypoints: MsgWaypoints, 
                                       dubins_path: DubinsParameters):
        ##### TODO #####
        # update halfspace variables

        self._halfspace_n = -self.dubins_path.n1
        print(self._halfspace_n)
        self._halfspace_r = self.dubins_path.r1

        self._path.plot_updated = False
        self._path.type = 'orbit'  
        self._path.orbit_center = self.dubins_path.center_s
        self._path.orbit_center[2] = waypoints.ned[2,self._ptr_previous]
        if self.dubins_path.dir_s == 1:
            self._path.orbit_direction = 'CW'
            print("CW")
        else:
            self._path.orbit_direction = 'CCW'
            print("CCW")
        self._path.orbit_radius = self.dubins_path.radius
        
        # quit()
        # Update path variables
        # self._path.__ =
        pass

    def _construct_dubins_line(self, 
                               waypoints: MsgWaypoints, 
                               dubins_path: DubinsParameters):
        ##### TODO #####
        # update halfspace variables
        # self._halfspace_n =
        # self._halfspace_r = 

        self._halfspace_n = self.dubins_path.n1
        print(self._halfspace_n)
        self._halfspace_r = self.dubins_path.r2
        self._path.plot_updated = False
        self._path.type = 'line' 
        self._path.line_direction = self.dubins_path.q1
        # print("path", self._path.line_direction)
        self._path.line_origin = self.dubins_path.r1 #waypoints.ned[:, self._ptr_previous]
        self._path.line_origin[2] = waypoints.ned[2,self._ptr_previous]
        # Update path variables
        # self._path.__ =
        pass

    def _construct_dubins_circle_end(self, 
                                     waypoints: MsgWaypoints, 
                                     dubins_path: DubinsParameters):
        ##### TODO #####
        # update halfspace variables
        self._halfspace_n = self.dubins_path.n3
        self._halfspace_r = self.dubins_path.r3
        print(self._halfspace_n)
        self._path.plot_updated = False
        self._path.type = 'orbit' 
        self._path.orbit_center = self.dubins_path.center_e
        if self.dubins_path.dir_e == 1:
            self._path.orbit_direction = 'CW'
        else:
            self._path.orbit_direction = 'CCW'
        self._path.orbit_radius = self.dubins_path.radius
        
        # Update path variables
        # self._path.__ =
        pass

    def _inHalfSpace(self, 
                     pos: np.ndarray)->bool:
        
        '''Is pos in the half space defined by r and n?'''
        # print((pos-self._halfspace_r).T @ self._halfspace_n)
        if (pos-self._halfspace_r).T @ self._halfspace_n >= 0:
            return True
        else:
            return False