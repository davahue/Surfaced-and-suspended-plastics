# -*- coding: utf-8 -*-
"""
Created on 25-05-2021 in Delft
@author: Antonio Moreno-Rodenas

This script showcases the in-camera color and shape based tracking of arbitrarily shaped blue plastic particles in an example video.

"""

import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import cv2.aruco as aruco
import imutils
import pandas as pd
import json 
from tqdm import tqdm
import math
import time
import collections

#%%
Camera = 1

pathExp = r'.\VideoExampletrackingV3'
Sync = pd.read_excel(r'.\ExtrinsicCalibration\DataFrameSync_Phase2_01_Transport_test _Cups(+ve).xls')

#%GET FLASH SYNC
Synked = Sync[Sync['Experiment'] == r'E:\01_Transport_test _Cups(+ve)\2_PP_undef_50\V3']
SynkedDelay = 0

#%% CAMERA 
pathCamera = os.path.join(pathExp, f'Cam {Camera}')
### Read Video
caps = []
framesnum = []
for videofile in np.sort(os.listdir(pathCamera)):
    capi = cv2.VideoCapture(os.path.join(pathCamera, videofile))
    caps.append(capi)
    framesnum.append(capi.get(cv2.CAP_PROP_FRAME_COUNT))

print('FrameCount:', framesnum)

#%%
# CREATE AOI MASK
ret, frame = capi.read()
h, w = frame.shape[:2]

AOI = np.zeros(frame.shape[:2])

#%% CALIBRATION DATA
#GET intrinsic Calibration
PathIntCal = r'IntrinsicCalibration\Camera {0}\CalibrationCAM{0}.json'.format(Camera)
with open(PathIntCal) as json_file:
    Cal = json.load(json_file)
    
IntMatrix = np.array(Cal['mtx'])
DistPar = np.array(Cal['dist'])

newcameramtx, roi = cv2.getOptimalNewCameraMatrix(IntMatrix,DistPar,(w,h),1,(w,h))

# GET extrinsic Calibration
### Aruco parameterization
aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_1000) # loading ARUCO fiduciary marker dictionary (4x4 code)
parameters_aruco =  aruco.DetectorParameters_create()

## Calibrate extrinsics
ArucoWorld_exterior = pd.read_excel(r"ExtrinsicCalibration\ArUcos_PhaseII.xlsx",
              sheet_name = 'outsideArucos',
              skiprows = [0,1,2],
              usecols=range(16),
              names = ['id','x0', 'y0', 'z0','x4', 'y4', 'z4', 'x3', 'y3', 'z3', 'x2', 'y2', 'z2', 'x1', 'y1', 'z1'], # Aruco coordinates of Biruk are oposite to the aruco convention, so we rename for convenience
              index_col = 'id')

ArucoWorld_interior = pd.read_excel(r"ExtrinsicCalibration\ArUcos_PhaseII.xlsx",
              sheet_name = 'insideArucos',
              skiprows = [0,1,2],
              usecols=range(16),
              names = ['id','x0', 'y0', 'z0', 'x4', 'y4', 'z4', 'x3', 'y3', 'z3', 'x2', 'y2', 'z2' ,'x1', 'y1', 'z1'],
              index_col = 'id')

ArucoWorld_exterior = ArucoWorld_exterior
ArucoWorld_interior = ArucoWorld_interior


## AUX FUNCTIONS
def Get_Extrinsic_frame_to_world_pointsOnlyCentre(ids, corners, WorldCoords):
    """
    Gets the corresponding points in the world and frame CRS for the markers coordinates:[centre, top_left, top_right, bottom_right, bottom_left]
    """
    idsInDic = [int(s.split('_')[1]) for s in WorldCoords.index]
    indexes = [int(np.where(ids == i)[0]) for i in idsInDic if i in ids]
    WorldCoordinates = []
    FrameCoordinates = []
    for ii in indexes:
        FrameCoordinates.append([corners[ii].mean(axis = 1)])
        WorldCoordinates.append(ArucoWorld_exterior.loc[f'id_{ids[ii][0]}', ['x0','y0','z0']].values.reshape(1,3))
    return np.array(WorldCoordinates).reshape(len(indexes),3,1).astype(float), np.array(FrameCoordinates).reshape(len(indexes),2).astype(float)

def ExtrinsicCalibration(frame_gray, frame_display):
    global aruco_dic, parameters_aruco
    # Extrinsic Calibration
    corners, ids, rejectedImgPoints = aruco.detectMarkers(frame_gray, aruco_dict, parameters=parameters_aruco)
    ids_f = []
    corners_f = []
    for i, id_i in enumerate(ids):
        if id_i in [13,14,15,16,17,11]: # Center Arucos
            pts = corners[i].astype(int).reshape((-1,1,2))
            cv2.circle(frame_display, tuple(pts.mean(axis = 0)[0].astype(int)), 1, (0,255,0), 2)
            cv2.putText(frame_display, f'Id_{id_i}', (pts[0,0][0], pts[0,0][1]), cv2.FONT_HERSHEY_SIMPLEX ,0.5, (0, 255, 0) , 2, cv2.LINE_AA) 
            ids_f.append(id_i)
            corners_f.append(corners[i])

    Corresponding_W2pixCoords = Get_Extrinsic_frame_to_world_pointsOnlyCentre(ids, corners, ArucoWorld_exterior)
    
    retval, rvec_ran, tvec_ran = cv2.solvePnP(Corresponding_W2pixCoords[0], Corresponding_W2pixCoords[1], IntMatrix, np.array([0.,0.,0.,0.,0.]))
    return rvec_ran, tvec_ran, ids, corners

def CheckSurfacedCountour(WaterLine, contour):
    '''
    Performs a fast check to define if the object detected touches the pre-defined water plane
    '''
    global surfaced
    ## Get highest point of contour
    highestPoint = contour[np.argmin(contour.reshape(len(contour), 2)[:,1])][0]

    ## Perform check up for position w/r to the water plane
    P1wl = np.array([WaterLine[0],WaterLine[1]])
    P2wl = np.array([WaterLine[2],WaterLine[3]])
    
    vP1_HP = highestPoint - P1wl
    vP1_P2 = P2wl - P1wl
    
    # CrossProduct
    xp = vP1_HP[0]*vP1_P2[1] - vP1_HP[1]*vP1_P2[0]
    
    if xp >= 0: # highest point is above the waterline
        surfaced = 1
    else:
        surfaced = 0
    return surfaced

def UPdateAOICAM1(ids, corners):
    '''
    Fefines an area of interest where we should expect particles
    '''
    id11 = corners[np.where(ids.T[0] == 11)[0][0]].mean(axis = 1)[0] + [70,-170]
    id13 = corners[np.where(ids.T[0] == 13)[0][0]].mean(axis = 1)[0] + [-70,250]

    AOI_1 = [id11[1], id13[1], id11[0], id13[0]]
    
    id11 = corners[np.where(ids.T[0] == 11)[0][0]].mean(axis = 1)[0] + [0,-134]
    id13 = corners[np.where(ids.T[0] == 13)[0][0]].mean(axis = 1)[0] + [0,-145]
    
    WaterLine = [int(id11[0]), int(id11[1]), int(id13[0]), int(id13[1])]
    
    return AOI_1, WaterLine

#%%
Results = []
StartingFrame = 0
SyncDelay = SynkedDelay
StartingFrame = StartingFrame + SyncDelay

cap = caps[0]
framesnumCum = [0] + list(np.cumsum(framesnum)[:-1])

#InitialCalibration
cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
ret, frame = cap.read()
frame_display = frame.copy()
frame = cv2.undistort(frame, IntMatrix, DistPar, None, newcameramtx)
frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
rot, trans, ids, corners = ExtrinsicCalibration(frame_gray, frame_display)

AOI_1, WaterLine = UPdateAOICAM1(ids, corners)

## Initialize 
VisualizationWidth = 1000 # select to adjust visualization window
Display = True # display video feedback during processing

frame_display = imutils.resize(frame, width = VisualizationWidth)
outVideo = cv2.VideoWriter(r'Example_output_detection_processed.avi',cv2.VideoWriter_fourcc('M','J','P','G'), 24, (frame_display.shape[1],frame_display.shape[0])) # records output

# Initialize params
FrameID = StartingFrame
FrameIDSync = 0

DTskip = 5 ## If the particle is not detected, skip n frames to accelerate processing
MemoryDetection = collections.deque(maxlen=30) # memory for particle detection
ObjectDetected = False

accumulatedFrames = 0
for cami, cam in enumerate(caps):
    
    cam.set(cv2.CAP_PROP_POS_FRAMES, FrameID)
    
    while True:

        ObjectDetected = False
        surfaced = 0
        ret, frame = cam.read()
        if ret:
            
            #Undistort
            frame = cv2.undistort(frame, IntMatrix, DistPar, None, newcameramtx)
            
            frame_display = frame.copy()
            frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
            frame_HSV = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV) 
        
            #Extrinsic calibration
            #Every 50 frames
            if (FrameID%50 == 0) | (FrameID < 50):
                rot, trans, ids, corners = ExtrinsicCalibration(frame_gray, frame_display)
                AOI_1, WaterLine = UPdateAOICAM1(ids, corners)

            else: #display previous location of markers
                for i, id_i in enumerate(ids):
                    if id_i in [13,14,15,16,17,11]: # Center Arucos
                        pts = corners[i].astype(int).reshape((-1,1,2))
                        cv2.circle(frame_display, tuple(pts.mean(axis = 0)[0].astype(int)), 1, (0,255,0), 2)
                        cv2.putText(frame_display, f'Id_{id_i}', (pts[0,0][0], pts[0,0][1]), cv2.FONT_HERSHEY_SIMPLEX ,0.5, (0, 255, 0) , 2, cv2.LINE_AA) 

            ##### FILTER FOR THE PARTICLE: Tune for the expected particle characteristics
            # color substraction:
            fgMask = cv2.inRange(frame_HSV, np.array([90,140,100]), np.array([105,255,255]))
            fgMask = cv2.blur(fgMask,(5,5),0)
        
            ### Detect Main Blob
            contours_Background, _ = cv2.findContours(fgMask, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
            
            if len(contours_Background) > 0:
                for i, contour in enumerate(contours_Background):            
                    ((x, y), radius) = cv2.minEnclosingCircle(contour)
                    area = cv2.contourArea(contour)
                    perimeter = cv2.arcLength(contour, True)
                    
                    ###*** AREA FILTER NOTE: SET UPPER LOWER BOUNDARY
                    if (area>500) & (area<100000):
                        circularity = 4*math.pi*(area/(perimeter*perimeter))
                             
                        ###*** formfactor FILTER:
                        if 0.2 < circularity < 1.2:
                            
                            M = cv2.moments(contour)
                            center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))
                            
                            # Check if in AOI:
                            if (AOI_1[0] < center[1]) and (center[1] < AOI_1[1]) and (AOI_1[2] < center[0]) and (center[0] < AOI_1[3]):
                                
                                ## FIT ELLIPSE
                                ellipse = cv2.fitEllipse(contour)
                                poly = cv2.ellipse2Poly((int(ellipse[0][0]), int(ellipse[0][1])), (int(ellipse[1][0] / 2), int(ellipse[1][1] / 2)), int(ellipse[2]), 0, 360, 5)

                                ## Draw axis of ellipse
                                pt1 = (int(ellipse[0][0]),int(ellipse[0][1]))
                                pt2 = (pt1[0] + int(np.sin(np.deg2rad(ellipse[2])) * ellipse[1][1]), pt1[1] - int(np.cos(np.deg2rad(ellipse[2])) * ellipse[1][1]))
                                pt3 = (pt1[0] + int(np.sin(np.deg2rad(ellipse[2] - 90)) * ellipse[1][0]), pt1[1] - int(np.cos(np.deg2rad(ellipse[2] - 90)) * ellipse[1][0]))
                                
                                cv2.line(frame_display, pt1, pt2, (0,255,255), 3)
                                cv2.line(frame_display, pt1, pt3, (0,255,0), 3)

                                ObjectDetected = True
                                
                                centerD = center
                                areaD = area
                                
                                cv2.drawContours(frame_display, contour, -1, (0, 255, 0), 3)
                                cv2.circle(frame_display, center, 5, (0,255,0), 2)

                                # Check if surfaced
                                surfaced = CheckSurfacedCountour(WaterLine, contour)
                                
            # display Waterline
            cv2.line(frame_display,(WaterLine[0],WaterLine[1]),(WaterLine[2],WaterLine[3]),(0,255,0),1)
            
            if surfaced:
                cv2.putText(frame_display, r'Surfaced', (centerD), cv2.FONT_HERSHEY_SIMPLEX, 1 , (0,0,255),2)  
                cv2.circle(frame_display, centerD, 5, (0,0,255), 2)
        
            if not ObjectDetected:
                centerD = [np.nan, np.nan]
                pt2 = [np.nan, np.nan]
                areaD = np.nan
                
            Results.append([FrameIDSync, centerD[0], centerD[1], pt2[0], pt2[1], areaD, *rot.T[0], *trans.T[0], surfaced])

            # Display
            cv2.putText(frame_display, r'Vid {} Frame {:.0f}, {:.0f}'.format(cami, cam.get(cv2.CAP_PROP_POS_FRAMES), FrameIDSync), (30, 30), cv2.FONT_HERSHEY_SIMPLEX, 1 , (0,255,0),2)         
            
            frame_display = imutils.resize(frame_display, width = VisualizationWidth)
            outVideo.write(frame_display)
            
            if Display:
                cv2.imshow('Frame', frame_display)
                
                if cv2.waitKey(10) & 0xFF == ord('q'):
                    FrameID = 0
                    break
                
            ### SKIP n frames if no particle is detected
            MemoryDetection.append(ObjectDetected)
            
            if (not any(MemoryDetection)) & (len(MemoryDetection) == 30):
                FrameID += DTskip
                cam.set(cv2.CAP_PROP_POS_FRAMES, FrameID)
            else:            
                FrameID = cam.get(cv2.CAP_PROP_POS_FRAMES)
                
            if FrameID >= framesnum[cami]:            
                FrameID = 0
                break
                
            FrameIDSync = FrameID - StartingFrame + framesnumCum[cami]
            
            print(FrameIDSync)
            
cv2.destroyAllWindows()
outVideo.release()

#%%
### SAVE RESULTS
Results = pd.DataFrame(Results, columns = ['FrameIDSync', 'centre_x', 'centre_y', 'mainAx_x', 'mainAx_y', 'Area','rot1', 'rot2', 'rot3', 'trans1', 'trans2', 'trans3', 'Surfaced'])
Results.to_csv(r'Example_output_detection.csv')















