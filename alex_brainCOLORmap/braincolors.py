#! /usr/bin/env python
"""                                                      
http://www.braincolor.org

braincolors.py takes in an Excel file with an adjacency matrix, 
where each value signifies adjacency between regions, and outputs the 
optimal assignment of colors to each group of regions on the command line, 
where optimal means maximally distinguishable colors within a neighborhood 
of similar colors in a color space:
 
1. Read in an Excel file with a binary (or weighted) adjacency matrix,
   where each row or column represents a region, and each value signifies 
   whether (or the degree to which) a given pair of regions are adjacent.
   Example: (a) column 0 = region abbreviation
            (b) column 1 & row 0 = full region name
            (c) column 2 = group number (each region is assigned to a group)
2. Create a colormap for the number of regions, with hues that are sampled
   from the (approx. perceptually uniform) CIELch cylindrical color space.
3. Convert the matrix from #1 to a graph, where each node represents a region
   and each edge represents the adjacency value between its connected nodes.
4. Break up the graph in #3 into subgraphs, where each subgraph represents
   a group of adjacent regions (assigned the same group number in #1c).
5. Compute every permutation of colors for the nodes of each subgraph in #4,
   with adjacent colors in the color space.
6. Assign each edge in each subgraph the value of the color difference 
   between the colors assigned to its pair of connected nodes in #5.
   (Multiply the connection matrix for each subgraph by
    the color difference matrix for each permutation.)
7. Find the optimal colors for the subgraph nodes that maximizes the sum 
   of the edge values from #6.
8. Plot the colormap, the whole graph, or individual colored subgraphs.
9. Optional: Replace RGB colors in an XML file.  
   The program recolor_eps.[csh,py] takes the output XML to recolor EPS files.

(c) Copyright 2010 . arno klein . arno@binarybottle.com . MIT license
"""

import sys
import xlrd
import re
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools
from colormath.color_objects import LCHuvColor,LabColor, sRGBColor
#from colormath.color_objects import 
from colormath.color_diff import delta_e_cie2000
from colormath.color_conversions import convert_color
from elementtree import ElementTree as et
import pdb


##################
# SET PARAMETERS #
##################

# Choose plotting procedure(s)
plot_colormap = 1  # Plot colormap
plot_graph = 0  # Plot whole graph
plot_subgraphs = 0  # Plot each individual colored subgraph

# Output color text file, save plots, and modify XML file
save_colors = 1
save_plots = 1
make_xml = 1

# Files
in_dir = 'input/'
out_dir = 'output/'
out_images = out_dir
out_colors = out_dir + 'region_colors.txt'
in_xml = in_dir + 'labels.xml'
out_xml = out_dir + 'labels.xml'
in_table = in_dir + 'region_adjacency_matrix.xls'
col_abbr = 0  # column with region abbreviations
col_group = 2  # column with region group numbers
col_start_data = 3  # first column with data
row_start_data = 1  # first row with data
everyother = 2  # use <everyother> alternate row(s);
                # set to 2 for redundant labels across, e.g. brain hemispheres
                    
# Color parameters
init_angle = 0 #22.5
chroma = 70  # color "saturation" level
Lumas_init = np.array([60,75,90])  # vary luminance values for adjacent colors
repeat_hues = 1  # repeat each hue for the different Lumas

# Graph layout
graph_node_size = 1000
graph_edge_width = 2
graph_font_size = 10
subgraph_node_size = 3000
subgraph_edge_width = 5
subgraph_font_size = 18
axis_buffer = 10

# Use weights in input adjacency matrix
use_input_weights = 0  

# Plot whole graph with colored subgraphs
plot_graph_color = 0  

#########
# BEGIN #
#########

# Convert weighted connection matrix to weighted graph
book = xlrd.open_workbook(in_table)
sheet = book.sheets()[0]
roi_abbrs = sheet.col_values(col_abbr)[1:sheet.ncols:everyother]
roi_abbrs = [str(s).strip() for s in roi_abbrs]
roi_numbers = sheet.col_values(2)[1:sheet.ncols:everyother] 
code_min = np.int(min(roi_numbers))
code_max = np.int(max(roi_numbers))
code_step = 1
    
iA = 0
A = np.zeros(((sheet.nrows-row_start_data)/everyother,(sheet.ncols-col_start_data)/everyother))
for irow in range(row_start_data,sheet.nrows,everyother):
    Arow = [s.value for s in sheet.row(irow)[col_start_data:]]
    A[iA] = Arow[0:len(Arow):everyother]
    iA += 1
A = A/np.max(A)  # normalize weights
G = nx.from_numpy_matrix(A)
Ntotal = G.number_of_nodes()
for inode in range(Ntotal):
    G.node[inode]['abbr'] = roi_abbrs[inode] 
    G.node[inode]['code'] = roi_numbers[inode]

# Secondary parameters
color_angle = 360.0 / Ntotal
if repeat_hues:
    nLumas = len(Lumas_init)
    nangles = np.ceil(Ntotal / np.float(nLumas))
    color_angle = nLumas * color_angle
else: 
    nangles = Ntotal
        
# Plot the colormap for the whole graph    
if plot_colormap:
    plt.figure(Ntotal+1,figsize=(5,10))
    # Define colormap as uniformly distributed colors in CIELch color space
    Lumas = Lumas_init.copy()
    while len(Lumas) < Ntotal: 
        Lumas = np.hstack((Lumas,Lumas_init))
    hues = np.arange(init_angle, init_angle + nangles*color_angle, color_angle)
    if repeat_hues:
        hues = np.hstack((hues * np.ones((nLumas,1))).transpose())
    for iN in range(Ntotal):
        ax = plt.subplot(Ntotal, 1, iN+1)
        plt.axis("off")
        lch = LCHuvColor(Lumas[iN], chroma, hues[iN]) #print(lch)
        #rgb = lch.convert_to('rgb', debug=False)
        rgb = convert_color(lch, sRGBColor)
        plt.barh(0,50,1,0, color=[rgb.rgb_r/255.,rgb.rgb_g/255.,rgb.rgb_b/255.])
    if save_plots:
        plt.savefig(out_images + "colormap.png")
     
# Plot graph
if plot_graph:
    plt.figure(Ntotal+2)
    labels={}
    for i in range(Ntotal):
        labels[i] = G.node[i]['abbr']
    pos = nx.graphviz_layout(G,prog="neato")
    nx.draw(G,pos,node_color='cyan',node_size=graph_node_size,width=graph_edge_width,with_labels=False)
    nx.draw_networkx_labels(G, pos, labels, font_size=graph_font_size, font_color='black')
    plt.axis('off')
    
if make_xml:
    tree = et.ElementTree(file=in_xml)
    
# Loop through subgraphs
if plot_graph_color + plot_subgraphs + make_xml + save_colors > 0:
    f = open(out_colors,'w')
    run_permutations = 1
    for code_start in range(code_min,code_max+code_step,code_step):   
        glist = [n for n,d in G.nodes_iter(data=True) \
                   if (np.int(d['code'])>=code_start) and \
                      (np.int(d['code'])<code_start+code_step)] 
        N = len(glist)
        if N > 0:
            g = G.subgraph(glist)
            
            # Define colormap as uniformly distributed colors in CIELch color space
            Lumas = Lumas_init.copy()
            while len(Lumas) < N: 
                Lumas = np.hstack((Lumas,Lumas_init))
        
            if repeat_hues:
                nangles_g = np.ceil(N / np.float(nLumas))
            else: 
                nangles_g = N
            hues = np.arange(init_angle, init_angle + nangles_g*color_angle, color_angle)
            if repeat_hues:
                hues = np.hstack((hues * np.ones((nLumas,1))).transpose())

            init_angle += nangles_g*color_angle
        
            # Compute the differences between every pair of colors in the colormap
            if run_permutations:
                # Convert subgraph into an adjacency matrix (1 for adjacent pair of regions)
                neighbor_matrix = np.array(nx.to_numpy_matrix(g,nodelist=glist))
                if use_input_weights:
                    pass
                else:
                    neighbor_matrix = (neighbor_matrix > 0).astype(np.uint8)
                
                # Compute permutations of colors and color pair differences
                DEmax = 0
                permutations = [np.array(s) for s in itertools.permutations(range(0,N),N)]
                permutation_max = np.zeros(N)
                for ipermutations in range(len(permutations)):
                    permutation = permutations[ipermutations]
                    color_delta_matrix = np.zeros(np.shape(neighbor_matrix))   
                    for i1, icolor1 in enumerate(permutation):
                        lch1 = LCHuvColor(Lumas[icolor1],chroma,hues[icolor1]) 
                        for i2, icolor2 in enumerate(permutation):
                            if (i2 > i1) and (neighbor_matrix[i1,i2] > 0):
                                lch2 = LCHuvColor(Lumas[icolor2],chroma,hues[icolor2])
                                #pdb.set_trace()
                                #DE = lch1.delta_e(lch2, mode='cie2000')
                                labCol1 = convert_color(lch1, LabColor)
                                labCol2 = convert_color(lch2, LabColor)
                                DE = delta_e_cie2000(labCol1, labCol2);
                                color_delta_matrix[i1,i2] = DE   
                    DE = np.sum((color_delta_matrix * neighbor_matrix))
                    # Store the color permutation with the maximum adjacency cost
                    if DE > DEmax:
                        DEmax = DE
                        permutation_max = permutation
                        #color_delta_matrix_max = color_delta_matrix
     
            # Color subgraphs
            if plot_graph_color:
                plt.figure(Ntotal+2)
                for iN in range(N):
                    ic = np.int(permutation_max[iN])
                    lch = LCHuvColor(Lumas[ic],chroma,hues[ic]) #print(lch)
                    #rgb = lch.convert_to('rgb', debug=False)
                    rgb = convert_color(lch, sRGBColor)
                    color = [rgb.rgb_r/255.,rgb.rgb_g/255.,rgb.rgb_b/255.]
                    nx.draw_networkx_nodes(g,pos,node_size=graph_node_size,nodelist=[g.node.keys()[iN]],node_color=color)

            # Draw a figure of the colored subgraph
            if plot_subgraphs or save_colors:
                if plot_subgraphs:
                    plt.figure(code_start)
                labels={}
                for iN in range(N):
                    labels[g.nodes()[iN]] = g.node[g.nodes()[iN]]['abbr']
                if plot_subgraphs:
                    pos = nx.graphviz_layout(g,prog="neato")
                    nx.draw(g,pos,node_size=subgraph_node_size,width=subgraph_edge_width,alpha=0.5,with_labels=False)
                    nx.draw_networkx_labels(g,pos,labels,font_size=subgraph_font_size,font_color='black')
                    plt.axis('off')
                print("")
                f.write("\n")
                for iN in range(N):
                    ic = np.int(permutation_max[iN])
                    lch = LCHuvColor(Lumas[ic],chroma,hues[ic]) #print(lch)
                    #rgb = lch.convert_to('rgb', debug=False)
                    rgb = convert_color(lch, sRGBColor)
                    color = [rgb.rgb_r/255.,rgb.rgb_g/255.,rgb.rgb_b/255.]
                    # Print optimal colors to the command line
                    print(g.node[g.nodes()[iN]]['abbr'] + ': ' + str([rgb.rgb_r,rgb.rgb_g,rgb.rgb_b]))
                    f.write(g.node[g.nodes()[iN]]['abbr'] + ': ' + str([rgb.rgb_r,rgb.rgb_g,rgb.rgb_b]) + '\n')
                    if plot_subgraphs:
                        nx.draw_networkx_nodes(g,pos,node_size=subgraph_node_size,nodelist=[g.node.keys()[iN]],node_color=color)
                if plot_subgraphs:
                    ax = plt.gca().axis()
                    plt.gca().axis([ax[0]-axis_buffer,ax[1]+axis_buffer,ax[2]-axis_buffer,ax[3]+axis_buffer])
                    if save_plots:
                        plt.savefig(out_images + "subgraph" + str(int(g.node[g.nodes()[0]]['code'])) + ".png")
                    plt.show()
                
            # Replace RGB colors in an XML file       
            """
            <LabelList>
            <Label>
              <Name>3rd Ventricle</Name>
              <Number>4</Number>
              <RGBColor>204 182 142</RGBColor>
            </Label>
            """
            if make_xml:
                for iN in range(N):
                    ic = np.int(permutation_max[iN])
                    lch = LCHuvColor(Lumas[ic],chroma,hues[ic]) #print(lch)
                    #rgb = lch.convert_to('rgb', debug=False)
                    rgb = convert_color(lch, sRGBColor)
                    color = [rgb.rgb_r, rgb.rgb_g, rgb.rgb_b]
                    color = ' '.join([str(s) for s in color])
                    for elem in tree.getiterator()[0]:
                        if g.node[g.nodes()[iN]]['abbr'] in elem.getchildren()[0].text:
                            elem.getchildren()[2].text = color
                            #print(g.node[g.nodes()[iN]]['abbr'],color)
                    
if plot_graph:
    if save_plots:
        plt.savefig(out_images + "graph.png")

if make_xml:
    tree.write(out_xml)
    
if save_colors:
    f.close()
