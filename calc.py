"""
This module contains the main calculations in the thesis.

Author: Kristian Alfsvag

Source code made during Master's thesis "Detecting Atmospheric rivers using persistent homology", Institute of Mathematics, University of Bergen, spring 2015.
"""


def findConnections(data, source, target, moveToEnd = lambda p: False):
	"""
	Uses persistent homology to find persistence diagram that hopefully gives a notion of how connected [:,source,:] is to [targetStart:,target[0],target[1]] in the filtered complex given by data.
	
	data:			3D array
	source:			integer in data.shape[1]
	target:			(integer, integer) in (data.shape[1], data.shape[2])
	moveToEnd:		If moveToEnd(p), then move point p to end of filtration. If removing these points gives a deformation retract of the full complex, then the result will be equivalent to just removing the points alltogether.
	"""
	import numpy as np
	import itertools
	import cfc
	from cells import CubicalCell, Simplex
	
	#Edit points so that points are moved to end of filtration.
	minVal = np.amin(data)
	newMinVal = np.nextafter(minVal, -float('inf') ) #Decrementing the value slightly.
	for p in itertools.product(*[range(d) for d in data.shape]):
		 if moveToEnd(p):
		 	data[p] = newMinVal
		
	
	#Edit points in source so that these are added first. (This will overwrite if some points in source were moved to the end.)
	maxVal = np.amax(data)
	
	for t in xrange(data.shape[0]):
		for i in xrange(data.shape[2]):
			data[t, source, i] = maxVal
				
	f = cfc.make2SkeletonFrom3Darray(data, filtDir='t')
	
	#Add cells to filtration making a connection between target and source
	#Do this by adding one extra vertex, a triangular surface between target and this point, and a single line from this point to Equator.
	
	ev = 'extra_vertex' #My implementation of cells need some implicit ordering of the vertices. Say that this vertex is smaller than all others, e.g. will always be first element in tuple of vertices.
	
	f.addCell( Simplex(ev, maxVal) ) #Important that this is not prioritized. Want this to appear in later filtration step than the source.
		
	#For all timesteps, add line from extra vertex to target with same filtration value as target. (Which is min of ev's and target's)
	for t in xrange(data.shape[0]): 
	
		line = Simplex( ( (t,target[0],target[1]), ev ), data[t,target[0],target[1]] ) #A line from target at time t to ev. Filtration value is value at target at time t.
		f.addCell(line, prioritize=True) #Prioritizing these will make the representatives chosen for the cycles include these if it is possible.
		
	#Add surfaces between these lines to eliminate unwanted cycles
	for t in xrange(data.shape[0]-1):
		filtVal = min(data[t,target[0],target[1]], data[(t+1,target[0],target[1])]) #Is to be added when last vertice is added
		triangle = Simplex( ( (t,target[0],target[1]), (t+1,target[0],target[1]), ev ), filtVal)
		f.addCell(triangle)
	
	#Add single line from source to extra vertex.
	source_vertex = (0,source,0) #One point on the source. Does not matter which.
	extra_line = Simplex( (source_vertex, ev), maxVal)
	extra_line.mark=True
	f.addCell( extra_line, prioritize=True ) #Line from source to ev. This will be first line with ev in boundary.
		
	f.computePersistentHomology(1, countMarks=True, trackProgression=True)
	d = f.getDiagram(1)
	#Only include cycles that include the artificial line from source to ev an odd number of times, so that we only get the rivers and not the noise. (Going back and forth means not going there at all, at least in Z/2.
	
	d.cycles = [c for c in d.cycles if f.markCount[f.index[c.birthCell]]%2 == 1 ]
		
	return d, f
		
	


def findAllRivers(year, minTime, maxTime, moveToEnd = lambda p: False):
	"""
	Try to find all atmospheric rivers starting at equator in [minTime, maxTime] ending in Bergen in a smaller interval determined by bergenFraction
	
	year
	minTime:	see rivers starting at equator after this timestep
	maxTime:	last timestep
	bergenFraction:	Fraction given as a pair of integers.
			Only see rivers hitting bergen after minTime + bergenFraction*(maxTime/minTime).
	
	return (d, f)	d is list of all 1-cycles going through the back wall, f is the filtration
	
	Method:
	Construct filtered complex representing space-time, and add artificial connection between Bergen in [minTimeBergen, maxTime] and equator, so that an actual connection will appear as a 1-cycle. Filtration is given by tcw, direction top-down (add highest values first).
	
	"""
	import time
	import numpy as np
	import ha,cfc,var
	from cells import CubicalCell, Simplex
	
	start = time.time()
	
	tcw = ha.loadTcw(year)
	
	data = ha.sliceArray(tcw,minLat=var.minLat, maxLat=var.maxLat, latStep=var.latStep, minLong=var.minLong, maxLong=var.maxLong, longStep=var.longStep, minTime=minTime, maxTime=maxTime)#, timeStep=(maxTime-minTime)/3)
	
	#Index corresponding to equator
	eqJ = var.j(0)
	#Indices in array corresponding to Bergen
	bergenJ=var.j(var.bergenLat)
	bergenI=var.i(var.bergenLong)
	
	d,f = findConnections(data, eqJ, (bergenJ, bergenI), moveToEnd=moveToEnd)
	end=time.time()
	print('It took '+str((end-start)/60)+' minutes')
	return d, f




def findAllRiversInMonth(year, month, extendFract = 0.0):
	"""
	Try to find all atmospheric rivers from equator to Bergen in month in year, using findAllRivers
	
	year
	month:		a number from 1-12. (Will only approximate month by dividing the year in twelve.)
			note: 1 (january) may or may not be a valid argument. depends on how it's implemented
	extendFract:	A float number between 0 and 1.
	
	Method:
	Construct filtered complex representing space-time, and add artificial connection between Bergen and equator, so that
	actual connection will appear as a 1-cycle. Filtration is given by tcw, direction top-down (add highest values first).
	
	"""
	#Timesteps: range(1464)
	#Month is a number from 2 to 12 (now this does not work with january)
	
	monthInterval = 1464/12
	numStepsToExtend = int(monthInterval*extendFract)
	minTimeBergen = monthInterval*(month-1)
	minTime = minTimeBergen - numStepsToExtend
	if minTime<0:
		raise Exception('Tried to extend out of year. Month: '+str(month))
	maxTime = monthInterval*(month) - 1 #-1 because my sliceArray-function is weird
	
	#Which points to move to end in filtration. 
	import var
	moveToEnd = lambda (t,j,i): t<numStepsToExtend and var.y(j)==var.bergenLat and var.x(i)==var.bergenLong  #var.y(j)*(maxTime-minTime)*bergenFract > (var.maxLat)*t

	
	return findAllRivers(year, minTime, maxTime, moveToEnd)

def Tbec(year, month):
	"""How i computed the Tbec (as far as I remember). Added 08.09.2016
	"""
	return findAllRiversInMonth(year, month, 0.3)
	

def skewComp(year, month):
	"""
	Compute persistent homology as in findAllRivers, with and without skewing.
	"""
	import var
	
	monthInterval = 1464/12
	
	fract = 0.3 
	ds, fs = findAllRiversInMonth(year, month, fract)
	d, f = findAllRiversInMonth(year, month, 0)
	return (ds, d), (fs,f)
	
	
	
def controlComp(year,month,extendFract=0.3):
	"""
	Compute persistent homology as in findAllRiversInMonth, but without altering the complex
	
	"""
	import time
	import ha,cfc,var, inout
	start = time.time()
	
	
	monthInterval = 1464/12
	numStepsToExtend = int(monthInterval*extendFract)
	minTimeBergen = monthInterval*(month-1)
	minTime = minTimeBergen - numStepsToExtend
	if minTime<0:
		raise Exception('Tried to extend out of year. Month: '+str(month))
	maxTime = monthInterval*(month) - 1 #-1 because my sliceArray-function is weird
	
	
	
	tcw = ha.loadTcw(year)
	data = ha.sliceArray(tcw,minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong,latStep=var.latStep, longStep=var.longStep,minTime=minTime, maxTime=maxTime)
	
	#test
	#import numpy as np
	#data = np.random.rand(20,20,20)
	
	f = cfc.make2SkeletonFrom3Darray(data, filtDir='t')
	f.computePersistentHomology(1)
	
	d = f.getDiagram(1)
	inout.save(d, 'No_cells_added/diagram_{:02d}_{:4d}'.format(month,year) )
	end = time.time()
	print('It took '+str((end-start)/60)+' minutes')
	
	return f.getDiagram(1)
	
	
	
def bec(year=2012):
	"""
	At every time-step, find out when Bergen becomes connected to equator in filtration starting with highest humidities.
	Does this over the atlantic ocean.
	Bergen is given by 60 degrees north, 5 degrees east
	Computed using H0
	"""
	import time
	import numpy as np
	import ha, var, cfc
	from cells import Simplex
	data = ha.loadTcw(year)
	#maxTime=1
	data = ha.sliceArray(data, minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong, latStep=var.latStep, longStep=var.longStep) #, maxTime=maxTime) #maxTime only for testing
	maxTime=data.shape[0]
	values = []
	start = time.time()
	for t in range(maxTime):
		print "Iteration number", t+1, "of", maxTime, "in year", year
		a=data[t,:,:]
		#Highest value in array, point that gets added first
		maxValue = np.max(a)
		maxValue = np.nextafter(maxValue, float('inf')) #Make maxvalue slightly larger
		#Edit equator so that this is added first.
		for ind in range(a.shape[1]):
			a[var.j(0),ind] = maxValue
		f = cfc.make1SkeletonFrom2Darray(a,'t')
		#Add extra vertex above Bergen. This will be very first cell in filtration.
		ev = 'extra_vertex'
		f.addCell(Simplex(ev, maxValue), prioritize=True)
		#Add line from Bergen to ev, with filtration value same as Bergen.
		bergen = (var.bergenJ, var.bergenI)
		f.addCell(Simplex( (bergen, ev), filtrationVale=a[bergen] ) )
		
		f.computePersistentHomology(0)
		d0 = f.getDiagram(0)
		
		import heapq
		#Get the two cycles that are born first.
		d = heapq.nlargest(2, d0, key=lambda c: c.birth)
		
		values = values + [min([c.death for c in d])] #bec is the filtration level were the two oldest 0-cycles merge.
	end = time.time()
	print "It took", (end-start)/60, "minutes."
	return values

	
	
	
def compTimeBec(month, year, save=True):
	"""Compare the time calculation to the bec, by calculating persistent homology induced by the bec.
	"""
	import inout, cfc, plot, heapq
	import matplotlib.pyplot as plt
	import numpy as np
	from persistencediagram import PersistenceDiagram
	#Load data
	tbec = inout.load('Time_diagrams/{:4d}/time_diagram_{:02d}_{:4d}'.format(year,month, year))
	bec = inout.load('bec/bec{:4d}'.format(year))
	timeInterval = len(bec)/12
	bec = bec[(month-1)*timeInterval:month*timeInterval]
	#Compute persistent homology given by bec
	f = cfc.makeCubicalComplex(np.array(bec), 't')
	f.computePersistentHomology(0)
	d = f.getDiagram(0)

	persistence=lambda c:c.persistence()
	for cycle in tbec.cycles:
		cycle.color = 'b'
	for cycle in d.cycles:
		cycle.color = 'r'
	bothDiagrams = PersistenceDiagram(d.cycles+tbec.cycles, 't')	
	bothDiagrams.cycles.sort(key = lambda cycle: cycle.persistence())
	#Cut off somewhere to get better looking diagrams
	numCycles = 2*len(d.cycles) #This is the number of cycles that are included in the picture. 
	bothDiagrams.cycles = bothDiagrams.cycles[-numCycles:]
	#Sort by birth
	bothDiagrams.cycles.sort(key=lambda cycle: cycle.getBirth())
	#Plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plot.plotDiagram(ax, bothDiagrams)
	for line in ax.get_lines():
		if line.get_color() == (0,0,1):
			tbec_line = line
			break
	for line in ax.get_lines():
		if line.get_color() == (1,0,0):
			bec_line = line
			break
	plt.legend([tbec_line,bec_line],['The tbec','Induced by bec'],loc='upper center', ncol=2)
	if save:
		plt.savefig(inout.path+'time_vs_bec_diagrams_{:4d}_{:02d}'.format(year,month), dpi=200, bbox_inches = 'tight')
		plt.clf()
	else:
		plt.show()


