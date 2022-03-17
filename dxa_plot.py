#!/usr/bin/env python

import numpy as np
import math
from numpy import linalg as LA
# import matplotlib.pyplot as plt
import ast


class dislocation_segment:
    def __init__(self, start, end, burgers):
        self.start = start
        self.end = end
        self.burgers = burgers
        self.mid = 0.5 * (self.start + self.end)
        self.length = LA.norm(self.end - self.start)
        # self.t = (self.end - self.start) / self.length
        print('projection z!!!!')
        self.t = [self.end[0]-self.start[0],0,self.end[2]-self.start[2]]

        b_outer_t = np.outer(self.burgers, self.t)

        self.color = l1(b_outer_t)


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


dislocation_data = []
burgers_data = []

# dislocation_file = open("dislocation100.txt", "r")
# burgers_file = open("dislocation100_burgers.txt", "r")
# dislocation_data_original = dislocation_file.read()
# burgers_data_original = burgers_file.read()
# dislocation_temp = list(ast.literal_eval(dislocation_data_original))
# burgers_temp = list(ast.literal_eval(burgers_data_original))
# for idx in range(len(dislocation_temp)):
#     dislocation_data.append(dislocation_temp[idx])
#     burgers_data.append(burgers_temp[idx])
# print('after 100 dislocation:', len(dislocation_data))
# print('after 100 burgers:', len(burgers_data))
#
# dislocation_file = open("dislocation110.txt", "r")
# burgers_file = open("dislocation110_burgers.txt", "r")
# dislocation_data_original = dislocation_file.read()
# burgers_data_original = burgers_file.read()
# dislocation_temp = list(ast.literal_eval(dislocation_data_original))
# burgers_temp = list(ast.literal_eval(burgers_data_original))
# for idx in range(len(dislocation_temp)):
#     dislocation_data.append(dislocation_temp[idx])
#     burgers_data.append(burgers_temp[idx])
# print('after 110 dislocation:', len(dislocation_data))
# print('after 110 burgers:', len(burgers_data))
#
dislocation_file = open("dislocation111.txt", "r")
burgers_file = open("dislocation111_burgers.txt", "r")
dislocation_data_original = dislocation_file.read()
burgers_data_original = burgers_file.read()
dislocation_temp = list(ast.literal_eval(dislocation_data_original))
burgers_temp = list(ast.literal_eval(burgers_data_original))
for idx in range(len(dislocation_temp)):
    dislocation_data.append(dislocation_temp[idx])
    burgers_data.append(burgers_temp[idx])
print('after 111 dislocation:', len(dislocation_data))
print('after 111 burgers:', len(burgers_data))
#
# dislocation_file = open("dislocationother.txt", "r")
# burgers_file = open("dislocationother_burgers.txt", "r")
# dislocation_data_original = dislocation_file.read()
# burgers_data_original = burgers_file.read()
# dislocation_temp = list(ast.literal_eval(dislocation_data_original))
# burgers_temp = list(ast.literal_eval(burgers_data_original))
# for idx in range(len(dislocation_temp)):
#     dislocation_data.append(dislocation_temp[idx])
#     burgers_data.append(burgers_temp[idx])
# print('after other dislocation:', len(dislocation_data))
# print('after other burgers:', len(burgers_data))

# print('dislocation_data=', dislocation_data)
# print('burgers_data=', burgers_data)
print()
loop_discrete = []
loop_discretize_num = 10

# for each dislocation line, use its vertices' information and burgers vector to create dislocation_segment object
for i in range(len(dislocation_data)):
    vertices_temp = dislocation_data[i]
    dislocation_vertices_temp = []
    num_of_dislocations = int(len(vertices_temp))

    dislocation_vertices = vertices_temp
    dislocation_burgers = burgers_data[i]

    ######### loop #########
    loop_coord = []
    for item in dislocation_vertices:
        loop_coord.append(item)

    # generate small segments along the loop and form a list
    for idx in range(num_of_dislocations - 1):
        start = loop_coord[idx]
        stop = loop_coord[idx + 1]
        numstep = loop_discretize_num
        start_x = start[0]
        start_y = start[1]
        start_z = start[2]
        stop_x = stop[0]
        stop_y = stop[1]
        stop_z = stop[2]
        segment_x = list(np.linspace(start_x, stop_x, numstep))
        segment_y = list(np.linspace(start_y, stop_y, numstep))
        segment_z = list(np.linspace(start_z, stop_z, numstep))
        # segment = np.linspace(start, stop, numstep)
        segment = []
        for i in range(len(segment_x)):
            segment.append([segment_x[i], segment_y[i], segment_z[i]])
        for idx_temp_seg in range(len(segment) - 1):
            temp_seg_discrete = dislocation_segment(np.array(segment[idx_temp_seg]),
                                                    np.array(segment[idx_temp_seg + 1]), dislocation_burgers)
            loop_discrete.append(temp_seg_discrete)

print('loop_discrete=', loop_discrete)
# print('len(loop_discrete)=',len(loop_discrete))
# for item in loop_discrete:
#     print(item.start,item.end)
#
# fig_temp = plt.figure()
# ax = fig_temp.add_subplot(111,projection='3d')
# for item in loop_discrete:
#     ax.scatter(item.start[0],item.start[1],item.start[2])
#
#
# ax.set_xlim(-95, 95)
# ax.set_ylim(-85, 85)
# ax.set_zlim(-125, 125)
#
# plt.show()


# projection
project_x = 0
project_y = 2

# definition for color map to generate color index for each pixel
# cell = np.array([[192.198, 0., 0., - 96.099],
#                  [0., 171.826, 0., - 85.913],
#                  [0., 0., 255.891, - 127.9455]])

# modify the size of investigated area to increase efficiency

x_min = -250
x_max = 250

y_min = -250
y_max = 250

z_min = -200
z_max = -180

loop_project_upperlimit_x = -100
loop_project_lowerlimit_x = 100
loop_project_upperlimit_y = -200
loop_project_lowerlimit_y = 50
loop_project_upperlimit_z = -200
loop_project_lowerlimit_z = -180

# x_min = -500
# x_max = 500
#
# y_min = -500
# y_max = 500
#
# z_min = -200
# z_max = -180

# x_min = cell[0, 3]
# x_max = cell[0, 3] + cell[0, 0]
#
# y_min = cell[1, 3]
# y_max = cell[1, 3] + cell[1, 1]
#
# z_min = cell[2, 3]
# z_max = cell[2, 3] + cell[2, 2]

print('x_min=', x_min)
print('x_max=', x_max)
print('y_min=', y_min)
print('y_max=', y_max)
print('z_min=', z_min)
print('z_max=', z_max)
mesh_number = 600  # 200 was used

# for seg in loop_discrete:
#     if seg.start[project_x] > x_max:
#         x_max = seg.start[project_x]
#     if seg.start[project_y] > y_max:
#         y_max = seg.start[project_y]
#     if seg.start[project_x] < x_min:
#         x_min = seg.start[project_x]
#     if seg.start[project_y] < y_min:
#         y_min = seg.start[project_y]
#
# print('x_max,y_max=', x_max, y_max)
# print('x_min,y_min=', x_min, y_min)
#
# enlarge_boundary = 1
# x_min = x_min * enlarge_boundary - 50
# x_max = x_max * enlarge_boundary + 50
# y_min = y_min * enlarge_boundary - 50
# y_max = y_max * enlarge_boundary + 50


x_mesh = np.linspace(x_min, x_max, mesh_number)
x_step = x_mesh[1] - x_mesh[0]
x_mesh_size = len(x_mesh)

y_mesh = np.linspace(y_min, y_max, mesh_number)
y_step = y_mesh[1] - y_mesh[0]
y_mesh_size = len(y_mesh)

# generate mormalization color intensity
color_min_normolize = 0
color_max_normolize = 0
for seg in loop_discrete:
    if abs(seg.color) > color_max_normolize:
        color_max_normolize = abs(seg.color)

color_min_normolize = color_max_normolize
for seg in loop_discrete:
    if abs(seg.color) < color_min_normolize:
        color_min_normolize = abs(seg.color)

first_mesh_range = 1
first_mesh_range_calculate = 1  # 20 was used

num_included_segs_matrix = np.zeros((x_mesh_size, y_mesh_size))
color = np.zeros((x_mesh_size, y_mesh_size))
# for each segment, determine it's color contribution
for seg1 in loop_discrete:
    index_x = int((seg1.start[project_x] - x_min) / x_step)
    index_y = int((seg1.start[project_y] - y_min) / y_step)

    # check if the segment is out of xmin xmax range
    if (index_x >= x_mesh_size) | (index_y >= y_mesh_size) | (index_x < 0) | (index_y < 0):
        continue

    color[index_x, index_y] = seg1.color
    num_included_segs_matrix[int(index_x), int(index_y)] = num_included_segs_matrix[
                                                               int(index_x), int(index_y)] + 1

color_copy = color.copy()
for index_x in range(x_mesh_size):
    for index_y in range(y_mesh_size):
        if color[index_x, index_y] != 0:
            start_x = index_x - first_mesh_range_calculate
            end_x = index_x + first_mesh_range_calculate
            x_range = []
            y_range = []
            if start_x < end_x:
                # unpack the result
                x_range.extend(range(start_x, end_x))
                # Append the last value
                x_range.append(end_x)

            start_y = index_y - first_mesh_range_calculate
            end_y = index_y + first_mesh_range_calculate
            if start_y < end_y:
                # unpack the result
                y_range.extend(range(start_y, end_y))
                # Append the last value
                y_range.append(end_y)

            # x_range = [*range(index_x - first_mesh_range_calculate, index_x + first_mesh_range_calculate, 1)]
            # y_range = [*range(index_y - first_mesh_range_calculate, index_y + first_mesh_range_calculate, 1)]
            for i in x_range:
                for j in y_range:
                    # print(index_x,index_y,i,j)
                    # if (i > 0) | (j > 0) | (i < x_mesh_size) | (j < y_mesh_size):
                    if (i < 0) | (j < 0) | (i >= x_mesh_size) | (j >= y_mesh_size):
                        continue
                    if (color[i, j] != 0) & (i != index_x) & (j != index_y):
                        # color[index_x, index_y] = color[index_x, index_y] + abs(color[i, j])
                        color[i, j] = color[index_x, index_y] + color_copy[i, j]
                        num_included_segs_matrix[i, j] = num_included_segs_matrix[i, j] + 1

# for index_x in range(x_mesh_size):
#     for index_y in range(y_mesh_size):
#         if color[index_x, index_y] != 0:
#             x_range = [*range(index_x - first_mesh_range_calculate, index_x + first_mesh_range_calculate, 1)]
#             y_range = [*range(index_y - first_mesh_range_calculate, index_y + first_mesh_range_calculate, 1)]
#             for i in x_range:
#                 for j in y_range:
#                     if (i > 0) | (j > 0) | (i < x_mesh_size) | (j < y_mesh_size):
#                         continue
#                     if color[i, j] != 0:
#                         color[index_x, index_y] = color[index_x, index_y] + color[i, j]
#                         num_included_segs_matrix[int(index_x), int(index_y)] = num_included_segs_matrix[
#                                                                                    int(index_x), int(index_y)] + 1

print('color=', color)
# for seg2 in loop_discrete:
#     deltax = seg1.start[project_x] - seg2.start[project_x]
#     deltay = seg1.start[project_y] - seg2.start[project_y]
#     if pow(deltax, 2) + pow(deltay, 2) < pow(first_mesh_range_calculate * x_step, 2) + pow(
#             first_mesh_range_calculate * y_step, 2):
#         color[index_x, index_y] = color[index_x, index_y] + seg2.color
#         num_included_segs_matrix[int(index_x), int(index_y)] = num_included_segs_matrix[
#                                                                    int(index_x), int(index_y)] + 1


for i in range(x_mesh_size):
    for j in range(y_mesh_size):
        if num_included_segs_matrix[i, j] != 0:
            color[i, j] = color[i, j] / num_included_segs_matrix[i, j]
# for i in range(x_mesh_size):
#     for j in range(y_mesh_size):
#         num_included_segs = 0
#         temp = 0
#         flag_found_seg = 0
#
#         # determine whether the pixel needs to be calculated for its color
#         for seg in loop_discrete:
#             if flag_found_seg == 1:
#                 break
#             if (pow((seg.start[project_x] - i * x_step - x_min), 2) + pow((seg.start[project_y] - j * y_step - y_min),
#                                                                           2)) < \
#                     (pow((first_mesh_range * x_step), 2) + pow((first_mesh_range * y_step), 2)):
#                 flag_found_seg = 1
#
#         # calculate color for such pixel
#         if flag_found_seg == 1:
#             for seg_calculate in loop_discrete:
#                 if (pow((seg_calculate.start[project_x] - i * x_step - x_min), 2) + pow(
#                         (seg_calculate.start[project_y] - j * y_step - y_min), 2)) < (
#                         pow((first_mesh_range_calculate * x_step), 2) + pow((first_mesh_range_calculate * y_step), 2)):
#                     temp = temp + (seg_calculate.color)
#                     num_included_segs = num_included_segs + 1
#
#         if num_included_segs != 0:
#             temp = abs(temp / num_included_segs)
#         else:
#             temp = 0
#         color.append(temp)


color = np.array(color)
color.resize(x_mesh_size, y_mesh_size)

# diffusing
span = 2
decay = 1
color_temp = color.copy()
for i in range(span, x_mesh_size - 1 - span):
    for j in range(span, y_mesh_size - 1 - span):
        for k in range(1, span):
            for l in range(1, span):
                color[i + k, j + l] = color_temp[i + k, j + l] + color_temp[i, j] \
                                      * math.exp(-decay * k ** 2) * math.exp(-decay * l ** 2)  # [100]
                color[i - k, j - l] = color_temp[i - k, j - l] + color_temp[i, j] \
                                      * math.exp(-decay * k ** 2) * math.exp(-decay * l ** 2)  # [100]
                color[i + k, j - l] = color_temp[i + k, j - l] + color_temp[i, j] \
                                      * math.exp(-decay * k ** 2) * math.exp(-decay * l ** 2)  # [100]
                color[i - k, j + l] = color_temp[i - k, j + l] + color_temp[i, j] \
                                      * math.exp(-decay * k ** 2) * math.exp(-decay * l ** 2)  # [100]

for i in range(span, x_mesh_size - 1 - span):
    for j in range(span, y_mesh_size - 1 - span):
        color[i, j] = abs(color[i, j])

# second mesh
second_mesh_range = 3
color_copy = color.copy()

print('x_mesh_size==', x_mesh_size)
print('y_mesh_size==', y_mesh_size)
print(color_copy.shape[0])
print(color_copy.shape[1])
for i in range(x_mesh_size):
    for j in range(y_mesh_size):
        temp = 0
        for k in range(1, second_mesh_range):
            for l in range(1, second_mesh_range):
                if (i + k > x_mesh_size - 1) | (j + l > y_mesh_size - 1) | (i - k < 1) | (j - l < 1):
                    continue
                temp = temp + abs(color_copy[i + k, j + l]) + abs(color_copy[i + k, j - l]) + abs(
                    color_copy[i - k, j + l]) + abs(color_copy[i - k, j - l])

        temp = temp + abs(color_copy[i, j])
        color[i, j] = temp

color_min = np.amin(color)
color_max = np.amax(color)
# color_max = 18.08073901035275  # xy
# color_max = 17.765063780448166   # xz
for i in range(x_mesh_size):
    for j in range(y_mesh_size):
        color[i, j] = color[i, j] - color_min
        color[i, j] = color[i, j] / color_max

print('final color_max 1 mesh size =', color_max)
print('color=', color)
print(type(color))

np.savetxt("color_mesh_20.txt", color)

y = np.loadtxt("color_mesh_20.txt")
#
# fig2 = plt.figure(figsize=(10, 6))
# ax5 = fig2.add_subplot(111)
#
# cmap = plt.get_cmap('jet')
# rgba_img = cmap(color)
# ax5.imshow(rgba_img)
#
# plt.show()
