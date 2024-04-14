# rrt dubins path planner for mavsim_python
import numpy as np
from message_types.msg_waypoints import MsgWaypoints
from planners.dubins_parameters import DubinsParameters


class RRTDubins:
    def __init__(self):
        self.segment_length = 450  # standard length of path segments
        self.dubins_path = DubinsParameters()

    def update(self, start_pose, end_pose, Va, world_map, radius):
        self.segment_length = 4 * radius
        tree = MsgWaypoints()
        tree.type = 'dubins'
        waypoints_not_smooth = MsgWaypoints()
        waypoints = MsgWaypoints()        

        ##### TODO #####
        # add the start pose to the tree
        print("START", start_pose[0:3,:])
        course = np.arctan2(end_pose.item(1) - start_pose.item(1),
                            end_pose.item(0) - start_pose.item(0)) #Get course
        # waypoints.add(start_pose[0:3,:], Va, start_pose[3,:])
        waypoints.add(start_pose[0:3,:], Va, [0])
        tree.add(start_pose[0:3,:], Va, course)
        
        # check to see if start_pose connects directly to end_pose
        self.extendTree(tree, end_pose, Va, world_map, radius)   

        tree.add(end_pose[0:3,:], Va, end_pose[3,:])
        # find path with minimum cost to end_node
        waypoints_not_smoothed = findMinimumPath(tree, end_pose)
        
        # find path with minimum cost to end_node
        # waypoints_not_smoothed = tree
        waypoints = waypoints_not_smoothed
        # waypoints = self.smoothPath(waypoints_not_smoothed, world_map, radius)
        # waypoints = tree
        print(waypoints.ned)
        self.waypoint_not_smooth = waypoints_not_smoothed
        self.tree = tree
        waypoints.type = 'dubins'
        return waypoints

    def extendTree(self, tree, end_pose, Va, world_map, radius):
        # extend tree by randomly selecting pose and extending tree toward that pose
        
        ##### TODO #####
        flag = False
        # while flag == False: 
        #Generate random configuration, find  closest leaf
        x = 0
        counter = 0
        
        while counter < 2:
            random = randomPose(world_map, end_pose.item(2))
            tmp = tree.ned-np.tile(random[0:3], (1, tree.num_waypoints)) #Distance between random and all leaves
            # print(tmp)
            tmp1 = np.diag(tmp.T @ tmp)
            idx = np.argmin(tmp1) #Get index of closest leaf
            # print(idx)
            dist = np.sqrt(tmp1.item(idx)) #Distance of closest leaf to new point
            # print(dist)
            L = np.max([np.min([dist, self.segment_length])]) #get the shorter length
            # print("L",L)
            
            cost = tree.cost.item(idx) + L #Add total cost
            tmp = random[0:3] - column(tree.ned, idx) #Get new point
            new_ned = column(tree.ned, idx) + L * (tmp / np.linalg.norm(tmp)) #Get NED coordinates
            new_chi = np.arctan2(new_ned.item(1) - tree.ned[1,idx],
                                    new_ned.item(0) - tree.ned[0,idx]) #Get course
            # new_chi = np.arctan2(end_pose.item(1) - tree.ned[1,idx],
            #                         end_pose.item(0) - tree.ned[0,idx]) #Get course
            
            # new_pose = np.concatenate((new_ned, np.array([[new_chi]])), axis = 0)
            tree_pose = np.concatenate((column(tree.ned, idx),np.array([[tree.course.item(idx)]])), axis = 0)
            
            if L > 2*radius:
                if self.collision(column(tree.ned, idx), new_ned, world_map, radius) == False:
                    
                    print("Yo")
                    if np.linalg.norm(new_ned[0:3,:] - end_pose[0:3,:]) <= self.segment_length:   
                        tree.add(new_ned, Va, parent = idx, course = new_chi, connect_to_goal = 1, cost = cost)
                        # print("Waypoints:", tree.num_waypoints)
                        # print(idx)
                        flag = True
                        counter += 1
                        print("# of Solutions:", counter)
                    else:
                        tree.add(new_ned, Va, parent = idx, course = new_chi, cost = cost, connect_to_goal = 0)  

        return flag

    def collision(self, start_pose, end_pose, world_map, radius):
        # check to see of path from start_pose to end_pose colliding with world_map
        
        ##### TODO #####
        collision_flag = False
        points = points_along_path(start_pose, end_pose, 500)
        for i in range(points.shape[1]):
            if heightAboveGround(world_map, column(points, i)) <= 0:
                collision_flag = True
        return collision_flag
    
    
    
    def process_app(self):
        self.planner_viewer.process_app()

    def smoothPath(self, waypoints, world_map, radius):
        
        ##### TODO #####
        # smooth the waypoint path
        smooth = [0]  # add the first waypoint
        # ptr = 1
        # while ptr <= waypoints.num_waypoints - 2:
        #     start_pose = np.concatenate((column(waypoints.ned, smooth[-1]), np.array([[waypoints.course]])))
        
        # construct smooth waypoint path
        smooth_idx=0
        good_idx = 1
        smooth_ned = np.array([[]])
        smooth_Va = []
        smooth_chi = []
        
        # construct smooth waypoint path
        smooth_waypoints = MsgWaypoints()
        smooth_waypoints.type = 'dubins'
        smooth_ned = np.append(smooth_ned, waypoints.ned[:,0]) #First waypoint
        smooth_ned = np.reshape(smooth_ned, (3,1))
        smooth_Va = np.append(smooth_Va, waypoints.airspeed[0])
        smooth_chi = np.append(smooth_chi, waypoints.course[0])
        
        last_smooth = np.zeros((3,1))

        #Get waypoint 
        while good_idx+2 <= waypoints.num_waypoints:
            # print(smooth_ned)
            last_smooth[:,0] = smooth_ned[:,-1]
            prev_point = np.reshape(waypoints.ned[:,good_idx], (3,1))
            next_point = np.reshape(waypoints.ned[:,good_idx+1], (3,1))
            
            #Check collision
            flag = self.collision(last_smooth, next_point, world_map, radius)
            #If no collision, update index get next waypoint
            if flag == False:
                good_idx += 1
                # print("Good",good_idx)
            else:
                smooth_ned = np.append(smooth_ned, prev_point, 1)
                smooth_Va = np.append(smooth_Va, waypoints.airspeed[good_idx])
                smooth_chi = np.append(smooth_chi, waypoints.course[good_idx])
                smooth_idx += 1

        for i in range(0,smooth_idx+1):
            smooth_waypoints.add(column(smooth_ned, i),
                            smooth_Va[i],
                            course = smooth_chi,
                            cost= np.inf,
                            parent = i-1,
                            connect_to_goal=1)
        
        smooth_waypoints.add(column(waypoints.ned, -1),
                            smooth_Va[i],
                            course=smooth_chi,
                            parent = smooth_idx,
                            connect_to_goal=1)
        # smooth_waypoints.add(column(waypoints.ned, -1),
        #                     smooth_Va[i],
        #                     course=smooth_chi,
        #                     parent = smooth_idx,
        #                     connect_to_goal=1)
        return smooth_waypoints

def points_along_path(start_pose, end_pose, N):
    # returns points along path separated by Del
    points = start_pose
    q = (end_pose - start_pose)
    L = np.linalg.norm(q)
    q = q / L
    w = start_pose
    for i in range(1, N):
        w = w + (L / N) * q
        points = np.append(points, w, axis=1)
    return points

def findMinimumPath(tree, end_pose):
    # find the lowest cost path to the end node

    # find nodes that connect to end_node
    connecting_nodes = []
    for i in range(tree.num_waypoints):
        if tree.connect_to_goal.item(i) == 1:
            connecting_nodes.append(i)
    # find minimum cost last node
    idx = np.argmin(tree.cost[connecting_nodes])
    # construct lowest cost path order
    path = [connecting_nodes[idx]]  # last node that connects to end node
    parent_node = tree.parent.item(connecting_nodes[idx])
    while parent_node >= 1:
        path.insert(0, int(parent_node))
        parent_node = tree.parent.item(int(parent_node))
    path.insert(0, 0)
    print("path",path)
    # construct waypoint path
    waypoints = MsgWaypoints()
    for i in path:
        waypoints.add(column(tree.ned, i),
                      tree.airspeed.item(i),
                      tree.course.item(i),
                      np.inf,
                      np.inf,
                      np.inf)
    waypoints.add(end_pose[0:3],
                  tree.airspeed[-1],
                  end_pose.item(3),
                  np.inf,
                  np.inf,
                  np.inf)
    # waypoints.add(column(tree.ned, 0),
    #                   tree.airspeed.item(0),
    #                   [np.pi*1.25],
    #                   np.inf,
    #                   np.inf,
    #                   np.inf)
    
    waypoints.type = tree.type
    return waypoints


def distance(start_pose, end_pose):
    # compute distance between start and end pose
    d = np.linalg.norm(start_pose[0:3] - end_pose[0:3])
    return d


def heightAboveGround(world_map, point):
    # find the altitude of point above ground level
    point_height = -point.item(2)
    tmp = np.abs(point.item(0)-world_map.building_north)
    d_n = np.min(tmp)
    idx_n = np.argmin(tmp)
    tmp = np.abs(point.item(1)-world_map.building_east)
    d_e = np.min(tmp)
    idx_e = np.argmin(tmp)
    if (d_n<world_map.building_width) and (d_e<world_map.building_width):
        map_height = world_map.building_height[idx_n, idx_e]
    else:
        map_height = 0
    h_agl = point_height - map_height
    return h_agl


def randomPose(world_map, pd):
    # generate a random pose
    pn = world_map.city_width * np.random.rand()
    pe = world_map.city_width * np.random.rand()
    chi = 2*np.pi * np.random.rand()
    pose = np.array([[pn], [pe], [pd], [chi]])
    return pose


def mod(x):
    # force x to be between 0 and 2*pi
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


def column(A, i):
    # extracts the ith column of A and return column vector
    tmp = A[:, i]
    col = tmp.reshape(A.shape[0], 1)
    return col