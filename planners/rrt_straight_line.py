# rrt straight line path planner for mavsim_python
import numpy as np
from message_types.msg_waypoints import MsgWaypoints

class RRTStraightLine:
    def __init__(self):
        self.segment_length = 300 # standard length of path segments

    def update(self, start_pose, end_pose, Va, world_map, radius):
        tree = MsgWaypoints()
        waypoints = MsgWaypoints()
        waypoints_not_smoothed = MsgWaypoints()
        #tree.type = 'straight_line'
        tree.type = 'fillet'

        ###### TODO ######
        # add the start pose to the tree
        waypoints.add(start_pose, Va)
        tree.add(start_pose, Va)
        # check to see if start_pose connects directly to end_pose

        self.extend_tree(tree, end_pose, Va, world_map)   
        # find path with minimum cost to end_node
        self.tree = tree
        waypoints_not_smoothed = find_minimum_path(tree, end_pose)
        waypoints = smooth_path(waypoints_not_smoothed, world_map)
        # waypoints = waypoints_not_smoothed
        # self.tree = tree
        self.waypoints_not_smoothed = waypoints_not_smoothed
        # print(waypoints_not_smoothed.ned)
        # print(waypoints.ned)
        return waypoints

    def extend_tree(self, tree, end_pose, Va, world_map):
        # extend tree by randomly selecting pose and extending tree toward that pose
        flag = False
        counter = 0
        # while flag == False: 
        #Generate random configuration, find  closest leaf
        x = 0
        while counter < 10:
            
            random = random_pose(world_map, end_pose.item(2))
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
            # new_chi = np.arctan2(new_ned.item(1) - tree.ned[1,idx],
            #                         new_ned.item(0) - tree.ned[0,idx]) #Get course
            # new_pose = np.concatenate((new_ned, np.array([[new_chi]])), axis = 0)
            # tree_pose = np.concatenate((column(tree.ned, idx),np.array([[tree.course.item(idx)]])), axis = 0)
            # x = x+1
            # print(x)
            if collision(column(tree.ned, idx), new_ned, world_map) == False:
                # print("Yo")
                if np.linalg.norm(new_ned - end_pose) <= self.segment_length:   
                    tree.add(end_pose, Va, parent = idx, connect_to_goal = 1, cost = cost)
                    # print("Waypoints:", tree.num_waypoints)
                    # print(idx)
                    flag = True
                    counter += 1
                    print("# of Solutions:", counter)
                else:
                    tree.add(new_ned, Va, parent = idx, cost = cost, connect_to_goal = 0)  

        # return flag
        
    # def process_app(self):
    #     self.planner_viewer.process_app()

def smooth_path(waypoints, world_map):

    ##### TODO #####
    # smooth the waypoint path
    smooth = [0]  # add the first waypoint
    smooth_idx=0
    good_idx = 1
    smooth_ned = np.array([[]])
    smooth_Va = []
    # construct smooth waypoint path
    smooth_waypoints = MsgWaypoints()
    smooth_waypoints.type = 'fillet'
    smooth_ned = np.append(smooth_ned, waypoints.ned[:,0]) #First waypoint
    smooth_ned = np.reshape(smooth_ned, (3,1))
    smooth_Va = np.append(smooth_Va, waypoints.airspeed[0])
    last_smooth = np.zeros((3,1))

    #Get waypoint 
    while good_idx+2 <= waypoints.num_waypoints:
        # print(smooth_ned)
        last_smooth[:,0] = smooth_ned[:,-1]
        prev_point = np.reshape(waypoints.ned[:,good_idx], (3,1))
        next_point = np.reshape(waypoints.ned[:,good_idx+1], (3,1))
        
        #Check collision
        flag = collision(last_smooth, next_point, world_map)
        #If no collision, update index get next waypoint
        if flag == False:
            good_idx += 1
            # print("Good",good_idx)
        else:
            smooth_ned = np.append(smooth_ned, prev_point, 1)
            smooth_Va = np.append(smooth_Va, waypoints.airspeed[good_idx])
            smooth_idx += 1

    for i in range(0,smooth_idx+1):
        smooth_waypoints.add(column(smooth_ned, i),
                        smooth_Va[i],
                        np.inf,
                        np.inf,
                        parent = i-1,
                        connect_to_goal=1)
    smooth_waypoints.add(column(waypoints.ned, -1),
                        smooth_Va[i],
                        parent = smooth_idx,
                        connect_to_goal=1)
    smooth_waypoints.add(column(waypoints.ned, -1),
                        smooth_Va[i],
                        parent = smooth_idx,
                        connect_to_goal=1)
    
    # smooth_waypoints.add(column(waypoints.ned, 0),
    #                     smooth_Va[i],
    #                     parent = smooth_idx,
    #                     connect_to_goal=1)
    # smooth_waypoints.add(column(waypoints.ned, -1),
    #                 waypoints.airspeed[-1],
    #                 np.inf,
    #                 np.inf,
    #                 np.inf,
    #                 1)
    #If collision, stop
    # smoothed_waypoints = waypoints
    return smooth_waypoints


def find_minimum_path(tree, end_pose):
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
    # construct waypoint path
    waypoints = MsgWaypoints()
    for i in path:
        waypoints.add(column(tree.ned, i),
                      tree.airspeed.item(i),
                      np.inf,
                      np.inf,
                      np.inf,
                      np.inf)
    waypoints.add(end_pose,
                  tree.airspeed[-1],
                  np.inf,
                  np.inf,
                  np.inf,
                  np.inf)
    waypoints.type = tree.type
    return waypoints


def random_pose(world_map, pd):
    # generate a random pose
    pn = world_map.city_width * np.random.rand()
    pe = world_map.city_width * np.random.rand()
    pose = np.array([[pn], [pe], [pd]])
    return pose


def distance(start_pose, end_pose):
    # compute distance between start and end pose
    d = np.linalg.norm(start_pose - end_pose)
    return d


def collision(start_pose, end_pose, world_map):
    # check to see of path from start_pose to end_pose colliding with map
    collision_flag = False
    points = points_along_path(start_pose, end_pose, 1000)
    for i in range(points.shape[1]):
        if height_above_ground(world_map, column(points, i)) <= 0:
            collision_flag = True

    return collision_flag


def height_above_ground(world_map, point):
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


def column(A, i):
    # extracts the ith column of A and return column vector
    tmp = A[:, i]
    col = tmp.reshape(A.shape[0], 1)
    return col