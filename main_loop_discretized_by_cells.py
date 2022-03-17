#!/usr/bin/env python
import scipy
from scipy import spatial
from scipy.spatial import ConvexHull

import numpy as np
import ast
import math
from numpy import linalg as LA


class Dislocation_Segment:
    def __init__(self, start, end, burgers):
        self.start = start
        self.stop = end
        self.start_voroid = None
        self.stop_voroid = None

        self.burgers = burgers
        self.mid = 0.5 * (self.start + self.stop)
        self.length = LA.norm(self.stop - self.start)
        self.t = (self.stop - self.start) / self.length
        b_outer_t = np.outer(self.burgers, self.t)
        self.color = l1(b_outer_t)
        self.gnd = gnd(b_outer_t)
        # identifies if the segment is cut by Voronoi cell boundaries
        self.iscut = 0
        self.idx_cell_belong = None

    def indexing_dislocation(self, index):
        self.index = index


def gnd(matrix):
    m = matrix.shape[0]
    n = matrix.shape[1]
    # temp used to calcualate l1 norm
    temp = []
    for i in range(n):
        temp.append(0)
        for j in range(m):
            temp[i] = temp[i] + matrix[j, i]
        temp[i] = abs(temp[i])

    compare = temp[0]
    for i in range(n):
        if temp[i] > compare:
            compare = temp[i]
    return compare


def l1(matrix):
    m = matrix.shape[0]
    n = matrix.shape[1]
    # temp used to calcualate l1 norm
    # temp_non_abs used to calculate l1 value without absolute value, to help nye tensor cancelling out
    temp = []
    temp_non_abs = []
    for i in range(n):
        temp.append(0)
        for j in range(m):
            temp[i] = temp[i] + matrix[j, i]
        temp_non_abs.append(temp[i])
        temp[i] = abs(temp[i])

    compare = temp[0]
    compare_non_abs = temp_non_abs[0]
    for i in range(n):
        if temp[i] > compare:
            compare = temp[i]
            compare_non_abs = temp_non_abs[i]
    return compare_non_abs


def removing_duplicates(points):
    # removing duplicates among a set of points
    print('removing duplicates among all points')
    x_points = []
    y_points = []
    z_points = []
    for point in points:
        x_points.append(str(round(point[0], 6)))
        y_points.append(str(round(point[1], 6)))
        z_points.append(str(round(point[2], 6)))

    d = {'x': x_points, 'y': y_points, 'z': z_points}
    df = pd.DataFrame(data=d)
    df = df.drop_duplicates()

    x_points_wo_duplicates = df['x'].tolist()
    y_points_wo_duplicates = df['y'].tolist()
    z_points_wo_duplicates = df['z'].tolist()

    points_wo_duplicates = []
    for idx, x_point in enumerate(x_points_wo_duplicates):
        temp = [float(x_point), float(y_points_wo_duplicates[idx]), float(z_points_wo_duplicates[idx])]
        points_wo_duplicates.append(temp)

    points_wo_duplicates = np.array(points_wo_duplicates)
    return points_wo_duplicates


def createSegmentUsingEndPoints(start, stop, loop_discretize_num):
    numstep = loop_discretize_num + 1
    start_x = start[0]
    start_y = start[1]
    start_z = start[2]
    stop_x = stop[0]
    stop_y = stop[1]
    stop_z = stop[2]
    segment_x = list(np.linspace(start_x, stop_x, numstep))
    segment_y = list(np.linspace(start_y, stop_y, numstep))
    segment_z = list(np.linspace(start_z, stop_z, numstep))
    segment = []
    for i in range(len(segment_x)):
        segment.append([segment_x[i], segment_y[i], segment_z[i]])
    # returns a list, with each item being the coordinates of the segment
    return segment


def dislocation_vertices_generation(dislocation_file, burgers_file, bound):
    print('######  start of dislocation_vertices_generation   ######')

    dislocation_data_withinBound = []
    burgers_data_withinBound = []

    # choose the dislocation type of interest to produce GND signal

    dislocation_data_original = dislocation_file.read()
    burgers_data_original = burgers_file.read()
    dislocation_lit_eval = list(ast.literal_eval(dislocation_data_original))
    burgers_lit_eval = list(ast.literal_eval(burgers_data_original))

    # for fully absorbed case, restrict  x: -112 ~ 65, y: -184 ~ 48, z: -198 ~ -187
    for idx in range(len(dislocation_lit_eval)):
        dislocation_data_withinBound.append(dislocation_lit_eval[idx])
        burgers_data_withinBound.append(burgers_lit_eval[idx])

    dislocation_seg_object_list = []
    discretize_num = 2

    #########################################
    # generate dislocation loop vertices
    #########################################

    # dislocation_data: a list of lists, each inner list represent a complete dislocation line, which can be discretized into dislocation segments
    for i, dislocation_seg_vertices in enumerate(dislocation_data_withinBound):
        # dislocation_seg_vertices -> burgers_data_withinBound[i]
        num_of_vertices = len(dislocation_seg_vertices)
        dislocation_seg_burgers = burgers_data_withinBound[i]

        # generate small segments along the loop and form a list
        for idx in range(num_of_vertices - 1):  # number of segment = num_of_vertices - 1
            start = dislocation_seg_vertices[idx]
            stop = dislocation_seg_vertices[idx + 1]

            seg = createSegmentUsingEndPoints(start, stop, discretize_num)
            seg_obj = Dislocation_Segment(np.array(seg[0]),
                                          np.array(seg[1]), dislocation_seg_burgers)

            dislocation_seg_object_list.append(seg_obj)

    # a list of coordinates for the vertices of dislocation lines, used for:
    # (1) the determination extremities for simulation box
    loop_vertices = []
    for seg_object in dislocation_seg_object_list:
        loop_vertices.append(seg_object.start.tolist())
    loop_vertices = np.array(loop_vertices)
    print('len loop_vertices=', len(loop_vertices))
    print('######  end of dislocation_vertices_generation   ######')
    print()
    return loop_vertices, dislocation_seg_object_list


# Face class is used for the combining of multiple coplanar faces of a polyhedron for better visualization and performance
class Faces():
    def __init__(self, tri, sig_dig=3, method="convexhull"):
        self.method = method
        self.tri = np.around(np.array(tri), sig_dig)
        self.grpinx = list(range(len(tri)))  # group index
        norms = np.around([self.norm(s) for s in self.tri], sig_dig)  #
        _, self.inv = np.unique(norms, return_inverse=True, axis=0)  # self.inv:indices of the unique array

    def norm(self, sq):
        cr = np.cross(sq[2] - sq[0], sq[1] - sq[0])
        return np.abs(cr / np.linalg.norm(cr))

    def isneighbor(self, tr1, tr2):
        a = np.concatenate((tr1, tr2), axis=0)
        return len(a) == len(np.unique(a, axis=0)) + 2

    def order(self, v):
        v = np.unique(v, axis=0)
        n = self.norm(v[:3])
        y = np.cross(n, v[1] - v[0])
        y = y / np.linalg.norm(y)
        c = np.dot(v, np.c_[v[1] - v[0], y])
        if self.method == "convexhull":
            h = scipy.spatial.ConvexHull(c)
            return v[h.vertices]
        else:
            mean = np.mean(c, axis=0)
            d = c - mean
            s = np.arctan2(d[:, 0], d[:, 1])
            return v[np.argsort(s)]

    def simplify(self):
        for i, tri1 in enumerate(self.tri):
            for j, tri2 in enumerate(self.tri):
                if j > i:
                    # if two tri are neighbors and their norm is the same, override bigger group index with smaller one
                    if self.isneighbor(tri1, tri2) and \
                            np.isclose(self.inv[i], self.inv[j], atol=1e-03):
                        self.grpinx[j] = self.grpinx[i]
        groups = []
        for i in np.unique(self.grpinx):  # iterate through unique group indecies
            u = self.tri[self.grpinx == i]
            u = np.concatenate([d for d in u])
            u = self.order(u)

            groups.append(u)

        return groups


def vertices_within(vertices, xmax, xmin, ymax, ymin, zmax, zmin):
    vertices_within = []
    for element in vertices:
        if element[0] > xmax + 0.2 * abs(xmax):
            continue
        if element[0] < xmin - 0.2 * abs(xmin):
            continue
        if element[1] > ymax + 0.2 * abs(ymax):
            continue
        if element[1] < ymin - 0.2 * abs(ymin):
            continue
        if element[2] > zmax + 0.2 * abs(zmax):
            continue
        if element[2] < zmin - 0.2 * abs(zmin):
            continue
        vertices_within.append(element.tolist())
    vertices_within = np.array(vertices_within)

    return vertices_within


def check_outside(points, xmax, xmin, ymax, ymin, zmax, zmin):
    outside = 0
    if any(points[:, 0] > xmax + 0.2 * abs(xmax)):
        outside = 1
    if any(points[:, 0] < xmin - 0.2 * abs(xmin)):
        outside = 1
    if any(points[:, 1] > ymax + 0.2 * abs(ymax)):
        outside = 1
    if any(points[:, 1] < ymin - 0.2 * abs(ymin)):
        outside = 1
    if any(points[:, 2] > zmax + 0.2 * abs(zmax)):
        outside = 1
    if any(points[:, 2] < zmin - 0.2 * abs(zmin)):
        outside = 1
    return outside


def extremity(points):
    #########################################
    # find max and min for data
    #########################################
    xmax = points[0][0]
    xmin = points[0][0]
    ymax = points[0][1]
    ymin = points[0][1]
    zmax = points[0][2]
    zmin = points[0][2]

    for point in points:
        if point[0] < xmin:
            xmin = point[0]
        if point[1] < ymin:
            ymin = point[1]
        if point[2] < zmin:
            zmin = point[2]
        if point[0] > xmax:
            xmax = point[0]
        if point[1] > ymax:
            ymax = point[1]
        if point[2] > zmax:
            zmax = point[2]

    return (xmax), (xmin), (ymax), (ymin), (zmax), (zmin)


def meshmask(color, x_mesh_size, y_mesh_size):
    mesh_range = 1  # 20 was used
    num_included_segs_matrix = np.zeros(x_mesh_size, y_mesh_size)

    color_copy = color.copy()
    for index_x in range(x_mesh_size):
        for index_y in range(y_mesh_size):
            if color[index_x, index_y] != 0:
                start_x = index_x - mesh_range
                end_x = index_x + mesh_range
                x_range = []
                y_range = []
                if start_x < end_x:
                    # unpack the result
                    x_range.extend(range(start_x, end_x))
                    # Append the last value
                    x_range.append(end_x)

                start_y = index_y - mesh_range
                end_y = index_y + mesh_range
                if start_y < end_y:
                    # unpack the result
                    y_range.extend(range(start_y, end_y))
                    # Append the last value
                    y_range.append(end_y)

                for i in x_range:
                    for j in y_range:
                        if (i < 0) | (j < 0) | (i >= x_mesh_size) | (j >= y_mesh_size):
                            continue
                        if (color[i, j] != 0) & (i != index_x) & (j != index_y):
                            color[i, j] = color[index_x, index_y] + color_copy[i, j]
                            num_included_segs_matrix[i, j] = num_included_segs_matrix[i, j] + 1

    for i in range(x_mesh_size):
        for j in range(y_mesh_size):
            if num_included_segs_matrix[i, j] != 0:
                color[i, j] = color[i, j] / num_included_segs_matrix[i, j]


def voronoiVerticesReadAndScale(vorvx_file_name, vornb_file_name, dislocation_vertices, expansion):
    print('read and scale Voronoi cell vertices')
    vorvx_file = open(vorvx_file_name, "r")
    vorvx_data_original = vorvx_file.read()
    vorvx_data = list(ast.literal_eval(vorvx_data_original))

    vornb_file = open(vornb_file_name, "r")
    vornb_data_original = vornb_file.read()

    xmax, xmin, ymax, ymin, zmax, zmin = extremity(dislocation_vertices)

    x_edge_verticecube = xmax - xmin
    y_edge_verticecube = ymax - ymin
    z_edge_verticecube = zmax - zmin
    center_vertices = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)]
    center_vorvx = [0.5, 0.5, 0.5]
    x_edge_ratio = x_edge_verticecube / 1
    y_edge_ratio = y_edge_verticecube / 1
    z_edge_ratio = z_edge_verticecube / 1
    vorvx_data_scaled = []
    for vorvx in vorvx_data:
        vorvx_scaled = []
        for vertice in vorvx:
            vertice_modified = [expansion * x_edge_ratio * (vertice[0] - center_vorvx[0]) + center_vertices[0],
                                expansion * y_edge_ratio * (vertice[1] - center_vorvx[1]) + center_vertices[1],
                                expansion * z_edge_ratio * (vertice[2] - center_vorvx[2]) + center_vertices[2]]
            vorvx_scaled.append(vertice_modified)
        vorvx_scaled = np.array(vorvx_scaled)
        vorvx_data_scaled.append(vorvx_scaled)

    # vorvx_data_scaled = []
    # for vorvx in vorvx_data:
    #     vorvx_scaled = []
    #     for vertice in vorvx:
    #         vertice_modified = [expansion * (vertice[0] * x_edge_verticecube + x_shift),
    #                             expansion * (vertice[1] * y_edge_verticecube + y_shift),
    #                             expansion * (vertice[2] * z_edge_verticecube + z_shift)]
    #         vorvx_scaled.append(vertice_modified)
    #     vorvx_scaled = np.array(vorvx_scaled)
    #     vorvx_data_scaled.append(vorvx_scaled)
    return vorvx_data_scaled


def associatePointsWithCells(subsegmentEndpoint, vorvx_data_scaled):
    idx_vorvx_answer = 0
    # nested loop determines which Voronoi cell the start point belongs to
    for idx_vorvx, vorvx in enumerate(vorvx_data_scaled):
        hull_start = scipy.spatial.ConvexHull(vorvx, incremental=True)
        volume_start1 = hull_start.volume

        trialPoint_start = np.array(subsegmentEndpoint)
        hull_start.add_points([trialPoint_start])
        volume_start2 = hull_start.volume

        # find start point coincidence first to save time
        if (volume_start1) == (volume_start2):
            idx_vorvx_answer = idx_vorvx

    return idx_vorvx_answer


def generateSubsegmentEndpoints(seg, cell_length_average):
    # linear search to find the number of cells within the span of the segment
    start = seg.start
    stop = seg.stop
    seg_length = np.linalg.norm(start - stop)
    num_discretize = math.ceil(seg_length / cell_length_average)

    # create subsegments based on the ratio between length of dislocation segment and characteristic cell length
    subsegment = createSegmentUsingEndPoints(start, stop, num_discretize)

    subsegment_endpoints = []
    for item in subsegment:
        subsegment_endpoints.append(item)
    return subsegment_endpoints


def addColorToCells(voro_cell_list, dislocation_seg_object_list):
    # add color of dislocation segment to voronoi cells
    print('add color of dislocation segment to voronoi cells')
    voro_color_list = []
    for voro_cell in voro_cell_list:
        for idx_included in voro_cell.dislocation_start_included:
            voro_cell.color = voro_cell.color + dislocation_seg_object_list[idx_included].color
        voro_cell.color = abs(voro_cell.color)
        voro_color_list.append(voro_cell.color)
    return voro_color_list


class Voronoi_Cell:
    def __init__(self, vorvx):
        self.vorvx = vorvx
        self.dislocation_start_included = []
        self.dislocation_end_included = []
        self.color = 0


def extractPlaneVertices(subsegment, vorvx_data_scaled):
    # for each seg, its partitioned into subsegments, and the vertices have corresponding cell IDs
    cell_id_0 = associatePointsWithCells(subsegment[0], vorvx_data_scaled)
    cell_id_1 = associatePointsWithCells(subsegment[1], vorvx_data_scaled)

    # reduce duplicates for cell0 within vorvx_eachcell
    cell_0_vorvx = []
    cell_0_vertice_sum = []  # use the sum of each vorvx data to check if there are duplicate vertices for each cell
    for vorvx_cell0 in vorvx_data_scaled[cell_id_0]:
        sum_round = round(vorvx_cell0[0] + vorvx_cell0[1] + vorvx_cell0[2], 4)
        if sum_round not in cell_0_vertice_sum:
            cell_0_vorvx.append(vorvx_cell0.tolist())
            cell_0_vertice_sum.append(sum_round)

    # reduce duplicates for cell1 within vorvx_eachcell
    cell_1_vorvx = []
    cell_1_vertice_sum = []  # use the sum of each vorvx data to check if there are duplicate vertices for each cell
    for vorvx_cell1 in vorvx_data_scaled[cell_id_1]:
        sum_round = round(vorvx_cell1[0] + vorvx_cell1[1] + vorvx_cell1[2], 4)
        if sum_round not in cell_1_vertice_sum:
            cell_1_vorvx.append(vorvx_cell1.tolist())
            cell_1_vertice_sum.append(sum_round)

    # commonVerticies: if the sum of vertices appear in both two cells, it is assumed that those vertices composes the common plane
    commonVerticies_sum = [value for value in cell_0_vertice_sum if value in cell_1_vertice_sum]
    idx_cell0 = []
    for commonVertice in commonVerticies_sum:
        idx_cell0.append(cell_0_vertice_sum.index(commonVertice))
    idx_cell1 = []
    for commonVertice in commonVerticies_sum:
        idx_cell1.append(cell_1_vertice_sum.index(commonVertice))

    # idx_cell0 and idx_cell1 are indices for each cell, indicating the common vertices.
    commonPlaneVertices = list(map(cell_0_vorvx.__getitem__, idx_cell0))

    return commonPlaneVertices


def calDistance_point2plane(a, b, c, d, point):
    # a,b,c,d are parameters for equation of plane: ax+by+cz+d=0
    distance = abs(a * point[0] + b * point[1] + c * point[2] + d) / math.sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))

    return distance


def calIntersection(commonPlaneVertices, subsegment):
    subsegment = np.array(subsegment)
    # pick three points in plane
    point1 = commonPlaneVertices[0]
    point2 = commonPlaneVertices[1]
    point3 = commonPlaneVertices[2]
    x1 = point1[0]
    y1 = point1[1]
    z1 = point1[2]
    x2 = point2[0]
    y2 = point2[1]
    z2 = point2[2]
    x3 = point3[0]
    y3 = point3[1]
    z3 = point3[2]

    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1

    # reference: https://www.geeksforgeeks.org/program-to-find-equation-of-a-plane-passing-through-3-points/
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)

    distance0 = calDistance_point2plane(a, b, c, d, subsegment[0])
    distance1 = calDistance_point2plane(a, b, c, d, subsegment[1])
    subsegment_vector = subsegment[1] - subsegment[0]

    intersectionPoint = subsegment[0] + distance0 / (distance0 + distance1) * subsegment_vector

    return intersectionPoint


def determine_iscut(dislocation_seg_object_list, vorvx_data_scaled):
    for seg in dislocation_seg_object_list:
        seg.start_voroid = associatePointsWithCells(seg.start, vorvx_data_scaled)

        # stop point has high probability to coincide with the same cell
        vorvx = vorvx_data_scaled[seg.start_voroid]
        hull_stop = scipy.spatial.ConvexHull(vorvx, incremental=True)
        volume_stop1 = hull_stop.volume

        trialPoint_stop = np.array(seg.stop)
        hull_stop.add_points([trialPoint_stop])
        volume_stop2 = hull_stop.volume

        if volume_stop1 == volume_stop2:
            seg.stop_voroid = seg.start_voroid
        else:
            seg.stop_voroid = associatePointsWithCells(seg.stop, vorvx_data_scaled)

        if seg.start_voroid != seg.stop_voroid:
            seg.iscut = 1


def main():
    # use dislocation information modify range of voronoi cells span
    dislocation_file = open("dislocation111.txt", "r")
    burgers_file = open("dislocation111_burgers.txt", "r")
    # scenarioFlag 1: for loop before absorption, restrict  x: -0 ~ 0, y: -0 ~ 0, z: -0 ~ -0
    # scenarioFlag 2: for loop partially in GB, restrict  x: -0 ~ 0, y: -0 ~ 0, z: -0 ~ -0
    # scenarioFlag 3: for loop in GB, restrict  x: -112 ~ 65, y: -184 ~ 48, z: -198 ~ -187

    scenarioFlag = 2
    if scenarioFlag == 1:  # for file 30240
        # bound = [[-80, 77], [-72, 66], [-140, -30]]
        bound = []

    if scenarioFlag == 2:  # for file 30240
        # bound = [[-106, 3], [-153, -30], [-184, -110]]
        bound = []

    if scenarioFlag == 3:
        bound = [[-112, 65], [-184, 48], [-198, -187]]

    dislocation_vertices, dislocation_seg_object_list = dislocation_vertices_generation(dislocation_file, burgers_file,
                                                                                        bound)

    gnd_sum = 0
    for dislocation in dislocation_seg_object_list:
        gnd_sum = gnd_sum + dislocation.gnd
    print('gnd_sum=', gnd_sum)

    # read and scale Voronoi cell vertices
    vorvx_file = "n100_vorvx1.txt"
    vornb_file = "n100_vornb1.txt"

    vorvx_data_scaled = voronoiVerticesReadAndScale(vorvx_file, vornb_file, dislocation_vertices, expansion=3)

    # a list storing the volume of all cells
    cell_list = []
    cell_volume_list = []
    for idx, vorvx in enumerate(vorvx_data_scaled):
        hull = scipy.spatial.ConvexHull(vorvx, incremental=True)
        cell_list.append(hull)
        cell_volume_list.append(hull.volume)

    cell_volume_average = 0
    for volume in cell_volume_list:
        cell_volume_average = cell_volume_average + volume

    cell_volume_average = cell_volume_average / len(cell_volume_list)
    cell_length_average = pow(cell_volume_average, 1 / 3)
    print('cell_volume_average=', cell_volume_average)
    print('cell_length_average=', cell_length_average)

    # associate start and stop point with Voronoi cells (initialize 'iscut' in the meantime)
    print('associate start and stop point with Voronoi cells')
    # determine which segs are cut
    determine_iscut(dislocation_seg_object_list, vorvx_data_scaled)

    # generate dislocation subsegments (subsegment) if the parent dislocation is cut by cell boundaries
    # use subsegments to determine the intersection point with cell walls
    print('generate sub-dislocation segments')
    # dislocation_seg_object_list_with_intersection is the new list used to calculate GND signal
    dislocation_seg_object_list_with_intersection = []
    for idx_seg, seg in enumerate(dislocation_seg_object_list):
        # complete segs are directly added to the new list
        if seg.iscut != 1:
            dislocation_seg_object_list_with_intersection.append(seg)

        # segs with intersections are partitioned into subsegments
        if seg.iscut == 1:
            subsegmentEndpoint_list = generateSubsegmentEndpoints(seg, cell_length_average)

            # create list of subsegments for intersection calculation
            # subsegment_list contains subsegments created by partitioning the seg of current for-loop
            subsegment_list = []
            for idx_subsegEndpoint, subsegEndpoint in enumerate(subsegmentEndpoint_list):
                # cell_id_list[idx_subsegEndpoint] is the ID of cell that contains subsegEndpoint
                subsegment_temp = []
                if idx_subsegEndpoint == len(subsegmentEndpoint_list) - 1:
                    break
                subsegment_temp.append(subsegEndpoint)
                subsegment_temp.append(subsegmentEndpoint_list[idx_subsegEndpoint + 1])
                subsegment_list.append(subsegment_temp)

            # seg_new_endPoints contains the endpoints of original seg and intersection points that are obtained
            # through loop below
            seg_new_endPoints = []
            seg_new_endPoints.append(seg.start)
            for subsegment in subsegment_list:
                # for each subsegment, find its corresponding cell first
                commonPlaneVertices = extractPlaneVertices(subsegment, vorvx_data_scaled)

                # then calculate intersection with common plane
                # if commonPlaneVertices is empty, it means the two endpoints of subsegments are not in neighboring cells, ignore this scenario
                if commonPlaneVertices:
                    intersectionPoint = calIntersection(commonPlaneVertices, subsegment)
                    seg_new_endPoints.append(intersectionPoint)

            seg_new_endPoints.append(seg.stop)

            # use new endpoints to generate new segments
            for idx in range(len(seg_new_endPoints) - 1):
                seg_new = Dislocation_Segment(seg_new_endPoints[idx], seg_new_endPoints[idx + 1], seg.burgers)

                # problem! some id_cell_new are identical
                dislocation_seg_object_list_with_intersection.append(seg_new)

                # num of segs = num of endpoints-1
                if idx == len(seg_new_endPoints) - 1:
                    break

    # calculate all segments' color contribution to cells
    print('len(dislocation_seg_object_list_with_intersection)=', len(dislocation_seg_object_list_with_intersection))
    print('len(dislocation_seg_object_list)=', len(dislocation_seg_object_list))

    # calculate segment-cell relation
    segment_cell_list = []
    for seg in dislocation_seg_object_list_with_intersection:
        seg.idx_cell_belong = associatePointsWithCells(seg.start, vorvx_data_scaled)
        segment_cell_list.append(seg.idx_cell_belong)


    print('len segment_cell_list=', len(segment_cell_list))
    print('segment_cell_list=', segment_cell_list)

    voro_cell_list = []
    for idx_vorvx, vorvx in enumerate(vorvx_data_scaled):
        voro_cell = Voronoi_Cell(vorvx)
        voro_cell_list.append(voro_cell)

        # determine each cell contains which dislocation segments
        for idx_seg, seg in enumerate(dislocation_seg_object_list_with_intersection):
            if (seg.start_voroid == idx_vorvx) & (seg.stop_voroid == idx_vorvx):
                voro_cell.dislocation_start_included.append(idx_seg)
                continue

    voro_color_list = addColorToCells(voro_cell_list, dislocation_seg_object_list_with_intersection)
    print('voro_color_list=', voro_color_list)
    print('len voro_color_list=', len(voro_color_list))
    print('len(hull_volume_list)=', len(cell_volume_list))

    f_color = open("plot_voro_color.dat", "w+")
    for color in voro_color_list:
        f_color.write('%f' % color)
        f_color.write('\n')

    f_volume = open("plot_voro_volume.dat", "w+")
    for volume in cell_volume_list:
        f_volume.write('%f' % volume)
        f_volume.write('\n')

    f_dislocation_cell_belong = open("plot_dislocation_seg_cell_belong.dat", "w+")
    for seg in dislocation_seg_object_list_with_intersection:
        flag_end = 0
        if seg == dislocation_seg_object_list_with_intersection[-1]:
            flag_end = 1
        line_start = ['%f,' % abs(seg.idx_cell_belong)]
        f_dislocation_cell_belong.writelines(line_start)
        if flag_end == 0:
            f_dislocation_cell_belong.write('\n')

    f_dislocation_color = open("plot_dislocation_seg_color.dat", "w+")
    for seg in dislocation_seg_object_list_with_intersection:
        flag_end = 0
        if seg == dislocation_seg_object_list_with_intersection[-1]:
            flag_end = 1
        line_start = ['%f,' % abs(seg.color)]
        f_dislocation_color.writelines(line_start)
        if flag_end == 0:
            f_dislocation_color.write('\n')

    f_dislocation_start = open("plot_dislocation_seg_start.dat", "w+")
    for seg in dislocation_seg_object_list_with_intersection:
        flag_end = 0
        if seg == dislocation_seg_object_list_with_intersection[-1]:
            flag_end = 1
        line_start = ['%f,' % seg.start[0], '%f,' % seg.start[1], '%f' % seg.start[2]]
        f_dislocation_start.writelines(line_start)
        if flag_end == 0:
            f_dislocation_start.write('\n')

    f_dislocation_end = open("plot_dislocation_seg_end.dat", "w+")
    for seg in dislocation_seg_object_list_with_intersection:
        flag_end = 0
        if seg == dislocation_seg_object_list_with_intersection[-1]:
            flag_end = 1
        line_end = ['%f,' % seg.stop[0], '%f,' % seg.stop[1], '%f' % seg.stop[2]]
        f_dislocation_end.writelines(line_end)
        if flag_end == 0:
            f_dislocation_end.write('\n')

    f_vorvx_scaled = open("plot_vorvx_scaled.dat", "w+")
    for vorvx in vorvx_data_scaled:
        for vertice in vorvx:
            line_vertice = ['%f,' % vertice[0], '%f,' % vertice[1], '%f' % vertice[2]]
            f_vorvx_scaled.writelines(line_vertice)
            f_vorvx_scaled.write('\n')
        f_vorvx_scaled.write('new_seg\n')


if __name__ == "__main__":
    main()
