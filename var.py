"""
Define variables used in calculation

Author: Kristian Alfsvag
"""

#Define window to work in
minLat=0
maxLat=65
minLong=-95
maxLong=15
latStep=2
longStep=2

#Coordinates of Bergen
bergenLat = 60
bergenLong = 5



import ha
#Functions to get indices from coordinates
i = lambda x:ha.indFromLong(x,minLong=minLong,longStep=longStep)
j = lambda y:ha.indFromLat(y,maxLat=maxLat,latStep=latStep)
#Functions to get coordinates from indices
x = lambda i:ha.longFromInd(i,minLong=minLong,longStep=longStep)
y = lambda j:ha.latFromInd(j,maxLat=maxLat,latStep=latStep)

bergenJ = j(bergenLat)
bergenI = i(bergenLong)

#The pressure levels the data is given on:
pressure_levels = [1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 400, 300, 200, 100]
