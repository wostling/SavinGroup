# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:37:46 2024

@author: willi
"""
import pandas as pd

df = pd.DataFrame([[1,2,3,4,5],
      [6,7,8,9,10],
      [11,12,13,14,15]])

def quartersplit(array, center):
    column_center, row_center= center[0], center[1]
    
    rows, columns = array.shape[0], array.shape[1]
    
    upperleft = 0
    upperright = 0
    lowerleft = 0
    lowerright = 0
    
    for row in range(rows):
        for column in range(columns):
            if (column<=column_center) & (row>row_center):
                lowerleft = lowerleft + array[column][row]
            elif (column>column_center) & (row>row_center):
                lowerright =lowerright + array[column][row]
            elif (column<=column_center) & (row<=row_center):
                upperleft = upperleft + array[column][row]
            else:
                upperright = upperright+ array[column][row]
    
    print("Upper left: " + str(upperleft))
    print("Upper right: " + str(upperright))
    print("Lower left: " + str(lowerleft))
    print ("Lower right: "+ str(lowerright))
    return upperleft, upperright, lowerleft, lowerright