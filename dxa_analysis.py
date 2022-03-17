#!/usr/bin/env python

import math
import ovito
from ovito import modifiers
from ovito.io import import_file, export_file
from ovito.modifiers import DislocationAnalysisModifier
from ovito.data import DislocationNetwork
from ovito.vis import Viewport
from ovito.vis import ParticlesVis
from ovito.modifiers import SelectTypeModifier, DeleteSelectedModifier
import os.path

from ovito.io.ase import ovito_to_ase


def write_cell_to_text(cell_data):
    cell_file_name = 'cell.txt'
    cell_file_exists = os.path.isfile(cell_file_name)
    if (cell_file_exists):
        os.remove(cell_file_exists)
    cell_file = open(cell_file_name, 'a')
    print('cell_data=',cell_data)
    print('cell_data=',type(cell_data))
    cell_file.write(cell_data)
    cell_file.close()

def write_dislocations_to_text(segment_temp, burgerstext, true_burgers_vector, idx_among):
    dislocation_file_name = 'dislocation{}.txt'.format(burgerstext)
    dislocationfile_exists = os.path.isfile(dislocation_file_name)

    # writing dislocation pivots into file
    if (dislocationfile_exists) & (idx_among == 0):
        os.remove(dislocation_file_name)

    disFile = open(dislocation_file_name, 'a')
    disFile.write('[')
    for idx_among_all_segments, list in enumerate(segment_temp):
        flag_end = 0

        # check if the end is reached
        if idx_among_all_segments == len(segment_temp) - 1:
            flag_end = 1
        string = '[' + str(list[0]) + ',' + str(list[1]) + ',' + str(list[2]) + ']'
        disFile.write(string)

        # if it's not the end, use comma to separate
        if flag_end == 0:
            disFile.write(',')

    disFile.write(']')
    disFile.write(',')
    disFile.close()

    # writing burgers vector into file
    burgers_file_name = 'dislocation{}_burgers.txt'.format(burgerstext)
    burgerfile_exists = os.path.isfile(burgers_file_name)

    if (burgerfile_exists) & (idx_among == 0):
        os.remove(burgers_file_name)

    burgerFile = open(burgers_file_name, 'a')
    string = '[' + str(true_burgers_vector[0]) + ',' + str(true_burgers_vector[1]) + ',' + str(
        true_burgers_vector[2]) + ']'
    burgerFile.write(string)
    burgerFile.write(',')
    burgerFile.close()


def main():
    # precision for burgers vector determination
    epsilon = 1e-3

    pipeline = import_file('before_absorption.xyz')
    # pipeline = import_file('bigerloop_10frame.xyz')
    # pipeline = import_file('frame19.xyz')
    print('import complete')
    ####### modifiers #######
    # Extract dislocation lines from a crystal with BCC structure:
    modifier_dxa = DislocationAnalysisModifier()
    modifier_dxa.input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC
    pipeline.modifiers.append(modifier_dxa)

    # selecting atoms
    modifier_select = SelectTypeModifier(property='Particle Type', types={'Fe'})
    pipeline.modifiers.append(modifier_select)
    print('atom selection complete')

    # deleting atoms so that dislocation can be viewed
    modifier_delete = DeleteSelectedModifier()
    pipeline.modifiers.append(modifier_delete)
    print('atom delete complete')

    data = pipeline.compute()
    print('data.cell[0]', '\n', data.cell[...])
    # write_cell_to_text(data.cell[...])


    total_line_length = data.attributes['DislocationAnalysis.total_line_length']
    cell_volume = data.attributes['DislocationAnalysis.cell_volume']
    print("Dislocation density: %f" % (total_line_length / cell_volume))

    dislocation_other = []
    dislocation_111 = []
    dislocation_100 = []
    dislocation_110 = []

    idx_among_100 = 0
    idx_among_110 = 0
    idx_among_111 = 0
    idx_among_other = 0

    # Print list of dislocation lines:
    print("Found %i dislocation segments" % len(data.dislocations.segments))
    for idx_among_all_segments, segment in enumerate(data.dislocations.segments):
        modulus_of_B = math.sqrt(pow(segment.true_burgers_vector[0], 2) + pow(segment.true_burgers_vector[1], 2) + pow(
            segment.true_burgers_vector[2], 2))

        segment_temp = segment.points.tolist()
        # [100]
        if abs(modulus_of_B - 1) < epsilon:
            idx_among_100 = idx_among_100 + 1
            print('writing 100')
            write_dislocations_to_text(segment_temp, '100', segment.true_burgers_vector, idx_among_100)
            dislocation_100.append(segment_temp)
            print(segment.true_burgers_vector)

        # [110]
        elif abs(modulus_of_B - math.sqrt(2)) < epsilon:
            print('writing 110')

            idx_among_110 = idx_among_110 + 1
            write_dislocations_to_text(segment_temp, '110', segment.true_burgers_vector, idx_among_110)
            dislocation_110.append(segment_temp)
            print(segment.true_burgers_vector)

        elif abs(modulus_of_B - (math.sqrt(3) / 2)) < epsilon:
            idx_among_111 = idx_among_111 + 1
            write_dislocations_to_text(segment_temp, '111', segment.true_burgers_vector, idx_among_111)
            dislocation_111.append(segment_temp)
            print(segment.true_burgers_vector)

        # other
        else:
            idx_among_other = idx_among_other + 1
            write_dislocations_to_text(segment_temp, 'other', segment.true_burgers_vector, idx_among_other)
            dislocation_other.append(segment_temp)
            print(segment.true_burgers_vector)

    print('number of 100: ', len(dislocation_100))
    print('number of 110: ', len(dislocation_110))
    print('number of 111: ', len(dislocation_111))
    print('number of other: ', len(dislocation_other))

    #
    # print('data.objects=', data.objects)
    # print('total num of dislocation detected = ',
    #       len(dislocation_100) + len(dislocation_110) + len(dislocation_111) + len(dislocation_other))

    pipeline.add_to_scene()
    data.cell.vis.rendering_color = (1, 0, 0)
    data.cell.vis.enabled = True

    vp_cell = Viewport()
    vp_cell.type = Viewport.Type.Perspective
    vp_cell.camera_pos = (-100, -150, 150)
    vp_cell.camera_dir = (2, 3, -3)
    vp_cell.fov = math.radians(60.0)
    vp_cell.zoom_all()
    vp_cell.render_image(size=(800, 600), filename="cell.png", background=(0, 0, 0), frame=8)

    print(data.attributes)
    # data.dislocations.vis.enabled = True

    vis_element = pipeline.source.data.particles.vis
    vis_element.shape = ParticlesVis.Shape.Square

    vp_dislocation = Viewport()
    vp_dislocation.type = Viewport.Type.Perspective
    vp_dislocation.camera_pos = (-100, -150, 150)
    vp_dislocation.camera_dir = (2, 3, -3)
    vp_dislocation.fov = math.radians(60.0)
    vp_dislocation.zoom_all()
    vp_dislocation.render_image(size=(800, 600), filename="dislocation.png", background=(0, 0, 0), frame=8)

    export_file(pipeline, "output.xyz", "xyz", columns=
    ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"]
                )

    # printing 100 dislocations
    # for segment in dislocation_100:
    #     print('[')
    #     for i in range(3):
    #         print(segment.points[i][0], ',', segment.points[i][1], ',', segment.points[i][2],',')
    #     print(']')

    # printing 100 dislocations
    # for segment in dislocation_110:
    #     print('[')
    #     for i in range(3):
    #         print(segment.points[i][0], ',', segment.points[i][1], ',', segment.points[i][2],',')
    #     print(']')
    #

if __name__ == "__main__":
    main()