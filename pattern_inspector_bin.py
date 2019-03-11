# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 18:03:49 2018

Description: finds the number of counts and their pattern in a single 2d-frame.
             Neighbours are all the remaining events in a frame relative to
             the event in question. The event in question is either the original
             event, i.e. the pixel in which the interaction took place, or a
             pixel that is adjacent to the original event and then, later on for
             bigger multi-events, adjacent events that are adjacent to previous
             adjacent pixels. 

@author: Kjell koch-mehrin
"""
import numpy as np


def pattern_inspector(frame, singles, doubles, triples, quads, gt_quads):

    # location of events in frame
    idx = np.where(frame > 0)
    eventcors = np.array([idx[0],idx[1]]).T # [row,col] format i.e. [y, x]
    
    # loop over each event
    while len(eventcors) > 0:
        # get coordinates of other events relative to respective event
        relneighbours = eventcors - eventcors[0]
        # get total absolute values for relative coordinates
        relneighbours_abs = np.sqrt(relneighbours**2).astype('int').sum(axis=1)
        
        # check if event is isolated
        if 1 not in relneighbours_abs:
            singles.append(frame[eventcors[0][0],eventcors[0][1]])
            # remove original event
            eventcors = eventcors[1:]
            
        # find pattern size of multi-event and sum all energies within multi
        else:
            # create holder for multi-event
            adjacentcors = []
            total_energy = []
             # add original pixel to total energy of event
            total_energy.append(frame[eventcors[0][0],eventcors[0][1]])
            # find index of neighbours that are adjacent
            adjacent_idx = relneighbours_abs == 1
            # append all neighbour event coordinates that are adjacent
            adjacentcors.append(eventcors[adjacent_idx])
            # remove adjacent events
            eventcors = eventcors[~adjacent_idx]
            # remove original event
            eventcors = eventcors[1:]
            
            # keep appending coordinates that are neighbours of adjacents, etc.
            while len(adjacentcors) > 0:
                # loop over first set of adjacent events
                for i in range(len(adjacentcors[0])):
                    # check if there are any other events left
                    if len(eventcors) > 0:
                        # get the relative coordinate of all neighbours to adjacent
                        relneighbours = eventcors - adjacentcors[0][i]
                        # get absolute values for relative position of neighbours
                        relneighbours_abs = np.sqrt(relneighbours**2).astype('int').sum(axis=1)
                        # get index of neighbours that are adjacent
                        adjacent_idx = relneighbours_abs == 1
                        # add the new adjacent neighbours set to the adjacent list
                        if len(adjacent_idx) > 0:
                            adjacentcors.append(eventcors[adjacent_idx])
                        # remove adjacent events of original adjacent
                        eventcors = eventcors[~adjacent_idx]
                    # add subject adjacent pixel to total energy of event
                    total_energy.append(frame[adjacentcors[0][i][0], adjacentcors[0][i][1]])
                # remove original adjacent event
                adjacentcors = adjacentcors[1:]
            
            # get summed energy of all charges and pattern type
            pattern = len(total_energy)
            if pattern == 2:
                doubles.append(sum(total_energy))
            elif pattern == 3:
                triples.append(sum(total_energy))
            elif pattern == 4:
                quads.append(sum(total_energy))
            elif pattern > 4:
                gt_quads.append(sum(total_energy))
        
    return singles, doubles, triples, quads, gt_quads
        