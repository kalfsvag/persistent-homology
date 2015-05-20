"""
Module for handling arrays.

Author: Kristian Alfsvag"""

from dynlib.metopen import metopen

	
def loadU(p = 600, year = 2012):
	filename = 'ei.ans.'+str(year)+'.' + str(p) + '.u'
	f, dat, grid = metopen(filename, 'u')
	return dat
	
def loadV(p = 600, year = 2012):
	filename = 'ei.ans.'+str(year)+'.' + str(p) + '.v'
	f, dat, grid = metopen(filename, 'v')
	return dat

def loadQ(p = 600, year = 2012):
	filename = 'ei.ans.'+str(year)+'.' + str(p) + '.q'
	f, dat, grid = metopen(filename, 'q')
	return dat
	
def loadTcw(year=2000):
	f, dat, grid = metopen('ei.ans.'+str(year)+'.sfc.tw', 'tcw')
	return dat
	
def sliceArray(array, minTime=0, maxTime=1463, timeStep=1, minLat=-90, maxLat=90, latStep=1, minLong=-180, maxLong=179.5, longStep=1):
	"""Slices out portion of array given by time, latitude and longitude."""
	t1 = minTime
	t2 = maxTime+1
	j1 = indFromLat(minLat)+1 #j1>j2 since indices here go in opposite direction from latitude
	j2 = indFromLat(maxLat)
	i1 = indFromLong(minLong)
	i2 = indFromLong(maxLong)+1
	
	return array[t1:t2:timeStep, j2:j1:latStep, i1:i2:longStep]	

def latFromInd(j, maxLat=90, latStep=1):
	"""Return latitude corresponding to index j in array loaded with metopen and then sliced with sliceArray."""
	return maxLat - j*latStep/2

def longFromInd(i, minLong=-180, longStep=1):
	"""Return longitude corresponding to index i in array loaded with metopen."""
	return minLong + i*longStep/2

def indFromLat(y, maxLat=90, latStep=1):
	"""Return index in array loaded with metopen corresponding to latitude y."""
	return int( (maxLat-y)*2 )/latStep

def indFromLong(x, minLong=-180, longStep=1):
	"""Return index in array loaded with metopen corresponding to longitude x."""
	return int( (x-minLong)*2 )/longStep

def loadPrecip(year):
	"""Return list of daily precipitation in Bergen in given year.
	Remember: Entry -1.0 means no precip."""
	import scipy.io
	
	loc = '/Data/gfi/spengler/kal045/data/'
	name = 'rr_'+str(year)+'0101-'+str(year)+'1231.mat'
	
	mat =  scipy.io.loadmat(loc+name)
	mat = mat['st50540'] #The part corresponding to Bergen
	mat=mat[0,0] #Don't know why I have to do this
	return mat[3][0] #This should be the list we're after
	
	



