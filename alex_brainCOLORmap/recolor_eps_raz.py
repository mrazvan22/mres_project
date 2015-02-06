#! /usr/bin/env python
#                                                      22 October 2010
#		Andrew J. Worth
#		andy@neuromorphometrics.com
#
#		Neuromorphometrics, Inc.
#		22 Westminster Street
#		Somerville, MA  02144-1630  USA
#
#		http://neuromorphometrics.com
#
#*********************************************************************
#
#  (c) Copyright 2010 Neuromorphometrics, Inc.   All rights reserved
#
#*********************************************************************

import sys
import string
import os
import re
import array
import math
import glob
import datetime
import numpy as np
from optparse import OptionParser
from xml.dom import minidom

usage = "%prog [options] inFile.eps outFile.eps\n\
\n\
This program identifies and recolors items in an eps file.\n\
\n\
Examples:\n\
  %prog -p --color 'RRR GGG BBB' inFile.eps preparedFile.eps \n\
  %prog --mapFile labelFileName.xml preparedFile.eps finalFile.eps"

parser = OptionParser(usage)

parser.add_option("-p", "--prepare",
		  action="store_true", dest="doPrepare", default=False,
                  help="prepare the file for coloring "
                  "[default: %default]")
parser.add_option("-c", "--color",
                  action="store", type="string", dest="findMeColor",
                  default='255 255 0',
                  help="color used to identify colored items "
                       "[default: %default]")
parser.add_option("-m", "--mapFile",
                  action="store", type="string", dest="mapFileName",
                  default='parcLabels.xml',
                  help="Name to color mapping file "
                       "[default: %default]")
parser.add_option("-d", "--csvFile",
                  action="store", type="string", dest="csvFileName",
                  default='GMM_matrix.csv',
                  help="Name to CSV matrix file of the variance for EBM, GMM or the 3 cluster matrices for Dirichet processes"
                       "[default: %default]")
parser.add_option("-o", "--overrideCount",
                  action="store", type="int", dest="overrideCount",
                  default='0',
                  help="Override number of colors found (for testing) "
                       "[default: %default]")

parser.add_option("-v", action="store_true", dest="verbose", default=False,
                  help="say what is going on "
                  "[default: %default]")
parser.add_option("-q", action="store_false", dest="verbose",
                  help="don't say what is going on")
(options, args) = parser.parse_args()

if len(args) != 2:
	parser.error("2 arguments are required")

# Get required arguments
inFileName  = args[0]
outFileName = args[1]

if options.verbose:
	print "inFileName  = %s"     % inFileName
	print "outFileName = %s"     % outFileName
	print "color       = \"%s\"" % options.findMeColor
	print "mapFile     = %s"     % options.mapFileName

# This program takes an eps file and can do one of two things:
# 1) Prepare an eps file to be recolored by identifying the items in
#    the eps file that are colored, or
# 2) Recolor a prepared eps file given a color mapping file.
#
# To prepare an eps file, first "eps2eps" is run on it to remove the
# cruft.  Then the file is written out for each colored item found with
# that colored item set to a given color (yellow by default), opened
# for viewing, and then asks for an identifier for that colored item.
# Once all desired colored items are identified, the eps file is written
# out with a comment indicating the identifier for each color.
#
# To recolor the prepared eps file, a mapping file is read that gives
# the color for each of the identifiers that might be found in the
# eps file comments.  The mapping file is XML like this:
#
# 	<?xml version="1.0"?>
# 	<LabelList>
#	<Label>
#	  <Name>Unlabeled</Name>
#	  <Number>1</Number>
#	  <RGBColor>0 0 0</RGBColor>
#	</Label>
#	...
#	<Label>
#	  <Name>Vitamin E Tablet</Name>
#	  <Number>74</Number>
#	  <RGBColor>236 217 151</RGBColor>
#	</Label>
#	</LabelList>
#
# The <Number> tag is not used here, but it is used by NVM!
#
# After reading the mapping file, the eps file is read and the colors
# in the mapping file are substituted for the original colors for each
# of the identifiers, and then the final recolored eps file is written
# out.

labelList = [ ]
labelNums = [ ]

#RED = (255, 0, 0)
#YELLOW = (255 255 0)
GRAY = '128 128 128'


# stages 6,12,18,24,36
ABNORMAL_ORDER = [170, 171, 132, 133, 154, 155, 106, 107, 200, 201, 144, 145, 122, 123, 180, 181, 202, 203, 152, 153, 102, 103, 118, 119, 172, 173, 166, 167, 190, 191, 160, 161, 128, 129, 168, 169, 142, 143, 198, 199, 195, 196, 184, 185]
stageIndexAbnormal = [0] + [ABNORMAL_ORDER.index(x) + 1 for x in [133, 133, 181, 169 ]]
STAGE_NR_LABELS = [6,12,18,24,36]

# stage values start from 0: ABETA142 = 0, PTAU181P = 1, TAU = 2, ..
NUMS_TO_EVENT = {170:6, 171:6, 132:11, 133:11, 154:18, 155:18, 106:19, 107:19, 200:20, 201:20, 144:21, 145:21, 122:22, 123:22, 180:23, 181:23, 202:24, 203:24, 152:26, 153:26, 102:27, 103:27, 118:28, 119:28, 172:30, 173:30, 166:31, 167:31, 190:32, 191:32, 160:33, 161:33, 128:34, 129:34, 168:35, 169:35, 142:36, 143:36, 198:37, 199:37, 195:39, 196:39, 184:40, 185:40}


CSV_MATRICES = ['EBM_matrix.csv', 'GMM_matrix.csv', 'cluster1_matrix.csv', 'cluster2_matrix.csv', 'cluster3_matrix.csv']
OUT_FOLDERS = [ "ordering_figures_final/images/%s" % x.split("_")[0] for x in CSV_MATRICES]


# 1st_ball=ABETA(event at index 0) 2nd_ball=PTAU(index 0) ... 7th_ball=FDG(index 12)
BALL_ABN_IND = [0,1,2,3,4,5,12]

NR_STAGES = len(stageIndexAbnormal)
NR_MATRICES = len(CSV_MATRICES)
NR_EVENTS = 42
NR_BALLS = 7


def getInterpolatedColor(abn_level):
  # abn_level = 0 -> yellow
  # abn_level = 1 -> red
  print "abn_level: %f" % abn_level
  x = int((1-abn_level)*255)
  print "x:%f" % x
  return "255 %d 0" % x

def create_images():
  if os.access(options.mapFileName,os.F_OK): # then file exists
    # Read color map file
    xmldoc = minidom.parse(options.mapFileName)
    xmlLabelList = xmldoc.getElementsByTagName("Label")
    for xmlLabel in xmlLabelList:
      name = xmlLabel.getElementsByTagName("Name")
      labelList.append(name[0].firstChild.data)
      uNumber = xmlLabel.getElementsByTagName("Number")
      number = int(uNumber[0].firstChild.data)
      labelNums.append(number)

  eventsAbnormalityAll = np.zeros([NR_MATRICES, NR_STAGES, NR_EVENTS], float)
  #for matrixName,outFolder in zip([CSV_MATRICES[0]], [OUT_FOLDERS[0]]):
  for matrixIndex in range(len(CSV_MATRICES)):
    matrixName = CSV_MATRICES[matrixIndex]
    outFolder = OUT_FOLDERS[matrixIndex]
    matrix = np.loadtxt(open(matrixName,"rb"),delimiter=",")
    for stageIndex in range(NR_STAGES):

      #print matrix
      # for each event get the sum of all the probabilities until the current stage
      eventsAbnormality = np.sum(matrix[:,:STAGE_NR_LABELS[stageIndex]],1)
      print eventsAbnormality

      assert(len(eventsAbnormality) == NR_EVENTS)
      eventsAbnormalityAll[matrixIndex, stageIndex, :] = eventsAbnormality


      labelColors = []
      redNums = ABNORMAL_ORDER[:stageIndexAbnormal[stageIndex]]
      yellowNums = set(ABNORMAL_ORDER) - set(redNums)
      for labelNum in labelNums:
        if labelNum in ABNORMAL_ORDER:
          color = getInterpolatedColor(eventsAbnormality[NUMS_TO_EVENT[labelNum]])
          print color
          labelColors.append(color)
        else:
          labelColors.append(GRAY)

      print "\nColoring file %s using matrix %s, stage %d ..." % (inFileName, matrixName, stageIndex)

      # Read input file and write ouput, changing the colors after
      #  finding '% recoloreps LABELSTRING' comments
      ff = open(inFileName,'r')
      contents = ff.read()
      ff.close()
      contentLines = contents.split('\n')

      # open output file for writing
      newOutFileName = "%s/stage_%d.eps" % (outFolder, STAGE_NR_LABELS[stageIndex])
      print newOutFileName
      of = open(newOutFileName,'w')
      #of = open(os.devnull,'w')

      skipNextLine = False
      for line in contentLines:
        if skipNextLine == False:
          h = re.compile('% recoloreps .*')
          hS = h.search(line)
          if hS: # if this is a color comment line
            if options.verbose:
              print 'Found color comment: ', line
            of.write(line+'\n')
            toFind = line[13:]
            if toFind in labelList:
              if options.verbose:
                print 'looking for color for ', toFind
              index = labelList.index(toFind)
              if options.verbose:
                print 'index is ', index
              color = labelColors[index]
              if options.verbose:
                print 'writing color: ', color
              of.write(color+' rG\n')
              skipNextLine = True
            else:
              if options.verbose:
                print toFind + ' is not in labelList\n'
          else: # something other than color, print it out
            of.write(line+'\n')
        else:
          if options.verbose:
            print 'skipped actual color line\n'
          skipNextLine = False

      of.close()

  # --------------------------------------------------------------------
  #
  if options.verbose:
    print "All done, bye."

  print eventsAbnormalityAll
  return eventsAbnormalityAll

  #sys.exit()

def create_latex(eventsAbnormalityAll):
  text = r'''
\documentclass[11pt,a4paper,oneside]{report}

\usepackage{float}
\usepackage{tikz}
\usetikzlibrary{plotmarks}
\usepackage{amsmath,graphicx}
\usepackage{epstopdf}
\usepackage[font=normal,labelfont=bf]{caption}
\usepackage{subcaption}
\usepackage{color}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{scalefnt}


% margin size
\usepackage[margin=1in]{geometry}

% tikz settings
\tikzstyle{state}=[circle,thick,draw=black, align=center, minimum size=2.1cm,
inner sep=0]
\tikzstyle{vertex}=[circle,thick,draw=black]
\tikzstyle{terminal}=[rectangle,thick,draw=black]
\tikzstyle{edge} = [draw,thick]
\tikzstyle{lo} = [edge,dotted]
\tikzstyle{hi} = [edge]
\tikzstyle{trans} = [edge,->]


\begin{document}
\belowdisplayskip=12pt plus 3pt minus 9pt
\belowdisplayshortskip=7pt plus 3pt minus 4pt

% scale parameter for the circles and the gradient
\tikzset{every picture/.append style={scale=0.6}}
% scale parameter for the upper and lower small brain images
\newcommand*{\scaleBrainImg}{0.3}

%col{x}{y}{z} respresents the color for ball z from matrix x at stage y (matrix x, stage y, ball z)
'''

  for matrixIndex in range(NR_MATRICES):
    for stageIndex in range(NR_STAGES):
      for ballIndex in range(NR_BALLS):
        # abnormality=0 -> yellow (greenComponent=1); abnormality=1 -> red (greenComponent=0)
        greenVal =  1 - eventsAbnormalityAll[matrixIndex, stageIndex, BALL_ABN_IND[ballIndex]]
        text += r'''
\definecolor{col''' + "%d%d%d" % (matrixIndex, stageIndex, ballIndex) + '''}{rgb}{1,''' + "%.3f" % greenVal + ''',0}'''

  text += '''\n\n '''
  for matrixIndex in range(NR_MATRICES):
    text += r'''
\begin{figure}[H]
  \centering'''
    for stageIndex in range(NR_STAGES):
      text += r'''
  %\begin{subfigure}[b]{0.15\textwidth}
    \begin{tikzpicture}[scale=1.0,auto,swap]

    % the two brain figures on top
    \node (upper_brain) at (0,1.5) { \includegraphics*[scale=\scaleBrainImg,trim=0 0 240 0]{'''

      text +=  "images/%s/stage_%d.eps" % (CSV_MATRICES[matrixIndex].split("_")[0], STAGE_NR_LABELS[stageIndex]) + r'''}};
    \node (lower_brain) at (0,-1.5) { \includegraphics*[scale=\scaleBrainImg,trim=240 0 0 0]{'''
      text +=  "images/%s/stage_%d.eps" % (CSV_MATRICES[matrixIndex].split("_")[0], STAGE_NR_LABELS[stageIndex]) + r'''}};

    % the 6 circles
    \draw[fill=''' + "col%d%d0" % (matrixIndex, stageIndex) + '''] (-1.6,-3.4) circle [radius=0.33cm] node {\scriptsize A};
    \draw[fill=''' + "col%d%d1" % (matrixIndex, stageIndex) + '''] (-0.7,-3.4) circle [radius=0.33cm] node {\scriptsize P};
    \draw[fill=''' + "col%d%d2" % (matrixIndex, stageIndex) + '''] (0.2,-3.4) circle [radius=0.33cm] node {\scriptsize T};
    \draw[fill=''' + "col%d%d3" % (matrixIndex, stageIndex) + '''] (-1.6,-4.2) circle [radius=0.33cm] node {\scriptsize C1};
    \draw[fill=''' + "col%d%d4" % (matrixIndex, stageIndex) + '''] (-0.7,-4.2) circle [radius=0.33cm] node {\scriptsize C2};
    \draw[fill=''' + "col%d%d5" % (matrixIndex, stageIndex) + '''] (0.2,-4.2) circle [radius=0.33cm] node {\scriptsize C3};

    % the big circle on the right
    \draw[fill=''' + "col%d%d6" % (matrixIndex, stageIndex) + '''] (1.3,-3.8) circle [radius=0.6cm] node {\scriptsize FDG};

    \end{tikzpicture}
  %\end{subfigure}
  % next subfigure
  \hspace{-1.5em}
  ~'''

    text += r'''
  \hspace{1em}
  % the red-to-yellow gradient on the right 
  \begin{tikzpicture}[scale=1.0,auto,swap]
    \shade[top color=red,bottom color=yellow] (0,0) rectangle (0.5,6);
    \node[inner sep=0] (corr_text) at (0.2,6.5) {1};
    \node[inner sep=0] (corr_text) at (0.2,-0.5) {0};
  \end{tikzpicture}
  \caption{''' + CSV_MATRICES[matrixIndex].split("_")[0] +'''}
\end{figure}

'''
  text += r'''

% scale parameter for the font size in the circles
\newcommand*{\scaleLabelImg}{0.7}
\begin{figure}[H]
  \centering
  {\scalefont{\scaleLabelImg}
  \begin{tikzpicture}[scale=1.0,auto,swap]

  % the two brain figures on top
  \node (brain) at (0,1.5) { \includegraphics*[scale=\scaleLabelImg]{images/Mid-Lateral_surface3.eps}};

  % the 6 circles
  \draw[scale=\scaleLabelImg] (-9,-6) circle [radius=2cm] node {ABETA142};
  \draw[scale=\scaleLabelImg] (-4,-6) circle [radius=2cm] node {PTAU181P};
  \draw[scale=\scaleLabelImg] (1,-6) circle [radius=2cm] node {TAU};
  \draw[scale=\scaleLabelImg] (-9,-10.5) circle [radius=2cm] node {RAVLT};
  \draw[scale=\scaleLabelImg] (-4,-10.5) circle [radius=2cm] node {MMSE};
  \draw[scale=\scaleLabelImg] (1,-10.5) circle [radius=2cm] node {ADAS13};

  % the big circle on the right
  \draw[scale=\scaleLabelImg] (8,-8.25) circle [radius=4cm] node {FDG-PET};

  \end{tikzpicture}
  }

  \caption{Labels of the different areas analysed}
\end{figure}

\end{document}

'''

  return text




eventsAbnormalityAll = create_images()
text = create_latex(eventsAbnormalityAll)

out = open(outFileName, 'w')
out.write(text)
out.close()
os.system("cd ordering_figures_final && pdflatex report1_gen.tex")
os.system("cd ordering_figures_final && xelatex report1_gen.tex")
print text
print eventsAbnormalityAll[1, :,2]
