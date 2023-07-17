#!/usr/bin/env python
import scipy
from scipy import spatial
from scipy.spatial import ConvexHull

import datetime
import gnd_functions
import os
from os import path



#
# # Face class is used for the combining of multiple coplanar faces of a polyhedron for better visualization and performance
# class Faces():
#     def __init__(self, tri, sig_dig=3, method="convexhull"):
#         self.method = method
#         self.tri = np.around(np.array(tri), sig_dig)
#         self.grpinx = list(range(len(tri)))  # group index
#         norms = np.around([self.norm(s) for s in self.tri], sig_dig)  #
#         _, self.inv = np.unique(norms, return_inverse=True, axis=0)  # self.inv:indices of the unique array
#
#     def norm(self, sq):
#         cr = np.cross(sq[2] - sq[0], sq[1] - sq[0])
#         return np.abs(cr / np.linalg.norm(cr))
#
#     def isneighbor(self, tr1, tr2):
#         a = np.concatenate((tr1, tr2), axis=0)
#         return len(a) == len(np.unique(a, axis=0)) + 2
#
#     def order(self, v):
#         v = np.unique(v, axis=0)
#         n = self.norm(v[:3])
#         y = np.cross(n, v[1] - v[0])
#         y = y / np.linalg.norm(y)
#         c = np.dot(v, np.c_[v[1] - v[0], y])
#         if self.method == "convexhull":
#             h = scipy.spatial.ConvexHull(c)
#             return v[h.vertices]
#         else:
#             mean = np.mean(c, axis=0)
#             d = c - mean
#             s = np.arctan2(d[:, 0], d[:, 1])
#             return v[np.argsort(s)]
#
#     def simplify(self):
#         for i, tri1 in enumerate(self.tri):
#             for j, tri2 in enumerate(self.tri):
#                 if j > i:
#                     # if two tri are neighbors and their norm is the same, override bigger group index with smaller one
#                     if self.isneighbor(tri1, tri2) and \
#                             np.isclose(self.inv[i], self.inv[j], atol=1e-03):
#                         self.grpinx[j] = self.grpinx[i]
#         groups = []
#         for i in np.unique(self.grpinx):  # iterate through unique group indecies
#             u = self.tri[self.grpinx == i]
#             u = np.concatenate([d for d in u])
#             u = self.order(u)
#
#             groups.append(u)
#
#         return groups



def main():
    case_string = 'voronoi'

    number_seeds_list = []
    n = int(input("Enter length of seeds' list : "))
    for i in range(0, n):
        ele = int(input())
        number_seeds_list.append(ele)  # adding the element

    for number_seeds in number_seeds_list:
        current_path = os.getcwd()
        target_path = current_path + '/' + case_string
        max_filenum = 1

        my_file = target_path + str(number_seeds) + '/n' + str(number_seeds) + '_vorvx' + str(
            max_filenum) + '.txt'
        print('my_file=', my_file)
        print('path.exists(my_file)=', path.exists(my_file))

        if path.exists(my_file):
            max_filenum = max_filenum + 1
            flag_max_filenum_found = 0
            my_file = target_path + str(number_seeds) + '/n' + str(number_seeds) + '_vorvx' + str(
                max_filenum) + '.txt'
        else:
            print('no file associated with this number of seed is in the directory')
            flag_max_filenum_found = 1

        while flag_max_filenum_found == 0:
            if path.exists(my_file):
                max_filenum = max_filenum + 1
                my_file = target_path + str(number_seeds) + '/n' + str(number_seeds) + '_vorvx' + str(
                    max_filenum) + '.txt'
            else:
                flag_max_filenum_found = 1

        print('max_filenum=', max_filenum)
        for ii in range(1, max_filenum):
            begin_time = datetime.datetime.now()

            # use dislocation information modify range of voronoi cells span
            dislocation_file = open('MDdata_' + case_string + '/dislocation111.txt', "r")
            burgers_file = open('MDdata_' + case_string + '/dislocation111_burgers.txt', "r")

            if 'after' in case_string:
                bound = [[-112, 65], [-184, 48], [-191.8, -187]]
            else:
                bound = None
            dislocation_vertices, dislocation_seg_object_list = gnd_functions.dislocation_vertices_generation(dislocation_file,
                                                                                                burgers_file, bound)
            print('len(dislocation_seg_object_list)=',len(dislocation_seg_object_list))

            gnd_sum = 0
            for dislocation in dislocation_seg_object_list:
                gnd_sum = gnd_sum + dislocation.gnd

            gnd_nonabs_sum = 0
            for dislocation in dislocation_seg_object_list:
                gnd_nonabs_sum = gnd_nonabs_sum + dislocation.gnd_nonabs

            total_length = 0
            for dislocation in dislocation_seg_object_list:
                total_length = total_length + dislocation.length

            # read and scale Voronoi cell vertices
            vorvx_file = case_string + str(number_seeds) + '/' + 'n{}_vorvx{}.txt'.format(number_seeds, ii)
            print('vorvx_file=', vorvx_file)
            vorvx_data_scaled = gnd_functions.voronoiVerticesRead(vorvx_file)

            # a list storing the volume of all cells
            cell_list = []
            cell_volume_list = []
            cell_length_list = []
            for idx, vorvx in enumerate(vorvx_data_scaled):
                hull = scipy.spatial.ConvexHull(vorvx, incremental=True)
                cell_list.append(hull)
                cell_volume_list.append(hull.volume)
                cell_length_list.append(pow(hull.volume, 1 / 3))

            cell_volume_average = 0
            for volume in cell_volume_list:
                cell_volume_average = cell_volume_average + volume

            mean_cell_length_list = sum(cell_length_list) / len(cell_length_list)
            cell_volume_average = cell_volume_average / len(cell_volume_list)
            cell_length_average = pow(cell_volume_average, 1 / 3)

            # associate start and stop point with Voronoi cells (initialize 'iscut' in the meantime)
            print('associate start and stop point with Voronoi cells')
            # determine which segs are cut
            gnd_functions.determine_iscut(dislocation_seg_object_list, vorvx_data_scaled)

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
                    subsegmentEndpoint_list = gnd_functions.generateSubsegmentEndpoints(seg, cell_length_average)

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
                        commonPlaneVertices = gnd_functions.extractPlaneVertices(subsegment, vorvx_data_scaled)

                        # then calculate intersection with common plane
                        # if commonPlaneVertices is empty, it means the two endpoints of subsegments are not in neighboring cells, ignore this scenario
                        if len(commonPlaneVertices) > 3:
                            intersectionPoint = gnd_functions.calIntersection(commonPlaneVertices, subsegment)
                            seg_new_endPoints.append(intersectionPoint)

                    seg_new_endPoints.append(seg.stop)

                    # use new endpoints to generate new segments
                    for idx in range(len(seg_new_endPoints) - 1):
                        seg_new = gnd_functions.Dislocation_Segment(seg_new_endPoints[idx], seg_new_endPoints[idx + 1], seg.burgers)

                        # problem! some id_cell_new are identical
                        dislocation_seg_object_list_with_intersection.append(seg_new)

                        # num of segs = num of endpoints-1
                        if idx == len(seg_new_endPoints) - 1:
                            break

            # associate start and stop points for each new segments with cells
            for seg in dislocation_seg_object_list_with_intersection:
                seg.start_voroid = gnd_functions.associatePointsWithCells(seg.start, vorvx_data_scaled)
                seg.stop_voroid = gnd_functions.associatePointsWithCells(seg.start, vorvx_data_scaled)

            # calculate segment-cell relation
            segment_cell_list = []
            for seg in dislocation_seg_object_list_with_intersection:
                seg.idx_cell_belong = gnd_functions.associatePointsWithCells(seg.start, vorvx_data_scaled)
                segment_cell_list.append(seg.idx_cell_belong)

            voro_cell_list = []
            for idx_vorvx, vorvx in enumerate(vorvx_data_scaled):
                voro_cell = gnd_functions.Voronoi_Cell(vorvx)
                voro_cell_list.append(voro_cell)

                # determine each cell contains which dislocation segments
                for idx_seg, seg in enumerate(dislocation_seg_object_list_with_intersection):
                    if seg.start_voroid != seg.stop_voroid:
                        print('found trouble')
                        start = seg.start
                        stop = seg.stop
                        middle_point = [0.5 * (start[0] + stop[0]), 0.5 * (start[1] + stop[1]),
                                        0.5 * (start[2] + stop[2])]
                        middle_point_vorvxid = gnd_functions.associatePointsWithCells(middle_point, vorvx_data_scaled)

                        if middle_point_vorvxid == idx_vorvx:
                            voro_cell.dislocation_start_included.append(idx_seg)

                for idx_seg, seg in enumerate(dislocation_seg_object_list_with_intersection):
                    if (seg.start_voroid == idx_vorvx) & (seg.stop_voroid == idx_vorvx):
                        voro_cell.dislocation_start_included.append(idx_seg)

            voro_gnd_list = gnd_functions.addGNDToCells(voro_cell_list, dislocation_seg_object_list_with_intersection)

            count_voro_gnd_nonzero = 0
            for voro_gnd in voro_gnd_list:
                if voro_gnd != 0:
                    count_voro_gnd_nonzero = count_voro_gnd_nonzero + 1

            f_dislocation_cell_belong = open(
                (case_string + str(number_seeds) + '/' + "{}_plot_dislocation_seg_cell_belong{}.dat".format(
                    number_seeds, ii)),
                "w+")
            for seg in dislocation_seg_object_list_with_intersection:
                flag_end = 0
                if seg == dislocation_seg_object_list_with_intersection[-1]:
                    flag_end = 1
                line_start = ['%f,' % abs(seg.idx_cell_belong)]
                f_dislocation_cell_belong.writelines(line_start)
                if flag_end == 0:
                    f_dislocation_cell_belong.write('\n')

            f_gnd = open((case_string + str(number_seeds) + '/' + "{}_plot_voro_gnd{}.dat".format(number_seeds, ii)),
                         "w+")
            for gnd in voro_gnd_list:
                f_gnd.write('%f' % gnd)
                f_gnd.write('\n')

            f_volume = open(
                (case_string + str(number_seeds) + '/' + "{}_plot_voro_volume{}.dat".format(number_seeds, ii)),
                "w+")
            for volume in cell_volume_list:
                f_volume.write('%f' % volume)
                f_volume.write('\n')

            f_dislocation_start = open(
                (case_string + str(number_seeds) + '/' + "{}_plot_dislocation_seg_start{}.dat".format(number_seeds,
                                                                                                      ii)), "w+")
            for seg in dislocation_seg_object_list:
                flag_end = 0
                if seg == dislocation_seg_object_list[-1]:
                    flag_end = 1
                line_start = ['%f,' % seg.start[0], '%f,' % seg.start[1], '%f' % seg.start[2]]
                f_dislocation_start.writelines(line_start)
                if flag_end == 0:
                    f_dislocation_start.write('\n')

            f_dislocation_end = open(
                (case_string + str(number_seeds) + '/' + "{}_plot_dislocation_seg_end{}.dat".format(number_seeds, ii)),
                "w+")
            for seg in dislocation_seg_object_list:
                flag_end = 0
                if seg == dislocation_seg_object_list[-1]:
                    flag_end = 1
                line_end = ['%f,' % seg.stop[0], '%f,' % seg.stop[1], '%f' % seg.stop[2]]
                f_dislocation_end.writelines(line_end)
                if flag_end == 0:
                    f_dislocation_end.write('\n')

            f_vorvx_scaled = open(
                (case_string + str(number_seeds) + '/' + "{}_plot_vorvx_scaled{}.dat".format(number_seeds, ii)), "w+")
            for vorvx in vorvx_data_scaled:
                for vertice in vorvx:
                    line_vertice = ['%f,' % vertice[0], '%f,' % vertice[1], '%f' % vertice[2]]
                    f_vorvx_scaled.writelines(line_vertice)
                    f_vorvx_scaled.write('\n')
                f_vorvx_scaled.write('new_seg\n')

            f_vorvx_scaled = open(
                (case_string + str(number_seeds) + '/' + "{}_cpu_time{}.dat".format(number_seeds, ii)), "w+")
            f_vorvx_scaled.writelines(str(datetime.datetime.now() - begin_time))
            f_vorvx_scaled.write('\n')

            f_config = open((case_string + str(number_seeds) + '/' + "{}_config{}.dat".format(number_seeds, ii)), "w+")
            f_config.writelines('cell_length_average = ' + str(cell_length_average) + '\n')
            f_config.writelines('cell_volume_average = ' + str(cell_volume_average) + '\n')
            f_config.writelines('gnd_sum = ' + str(gnd_sum) + '\n')
            f_config.writelines('gnd_nonabs_sum = ' + str(gnd_nonabs_sum) + '\n')
            f_config.writelines('total_length = ' + str(total_length) + '\n')
            f_config.writelines('number_seeds = ' + str(number_seeds) + '\n')

            print('cell_length_average = ', cell_length_average)
            print('cell_volume_average = ', cell_volume_average)
            print('gnd_sum = ', gnd_sum)
            print('gnd_nonabs_sum = ', gnd_nonabs_sum)
            print('total_length = ', total_length)
            print('num_seeds = ', number_seeds)


if __name__ == "__main__":
    main()
