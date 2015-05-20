"""
Module with explicit calculations.

Author: Kristian Alfsvag
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
	d = f.getDiagrams(1)
	#Only include cycles that include the artificial line from source to ev an odd number of times, so that we only get the rivers and not the noise. (Going back and forth means not going there at all, at least in Z/2.
	
	d.cycles = [c for c in d.cycles if f.markCount[f.index[c.birthCell]]%2 == 1 ]
		
	return d, f
		
	


def findAllRivers(year, minTime, maxTime, moveToEnd = lambda p: False):
	"""
	Try to find all atmospheric rivers starting at equator in [minTime, maxTime] ending in Bergen in a smaller interval determined by bergenFraction
	
	year
	minTime:		see rivers starting at equator after this timestep
	maxTime:		last timestep
	bergenFraction	Fraction given as a pair of integers.
					Only see rivers hitting bergen after minTime + bergenFraction*(maxTime/minTime).
	
	return (dNoNoise, f)	dNoNoise is list of all 1-cycles going through the artificial connection
							f is the filtration
	
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
	month:			a number from 1-12. (Will only approximate month by dividing the year in twelve.)
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
	
	
	
def controlComp(year):
	"""
	Compute persistent homology as in findAllRivers, but without altering the complex
	
	"""
	import time
	import ha,cfc,var
	
	start = time.time()
	
	tcw = ha.loadTcw(year)
	data = ha.sliceArray(tcw,minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong,latStep=var.latStep, longStep=var.longStep)
	
	#test
	#import numpy as np
	#data = np.random.rand(20,20,20)
	
	f = cfc.make2SkeletonFrom3Darray(data, filtDir='t')
	f.computePersistentHomology(1)

	end = time.time()
	print('It took '+str((end-start)/60)+' minutes')
	
	return f.getDiagrams(1)
	

def check():
	"""
	For doing control examples
	"""
	
	import numpy as np
	import cfc
	import math
	
	
	tMax=12
	jMax=12
	iMax=1
	
	S = [(tMax-1,j) for j in range(jMax)]
	S.extend( [(tMax/2, j) for j in range(jMax/2,jMax)] )
	S.extend( [(0, j) for j in range(jMax/2)] ) 
	S.extend( [(t, jMax/2) for t in range(tMax-1)] )
	

	
	T,J,I = np.mgrid[0:tMax,0:jMax,0:iMax]
	
	#a is given by a(p) = d(p,S) (the distance to the set S)
	a = np.array([(T-x[0])**2 + (J-x[1])**2 for x in S])
	a = np.amin(a, axis=0)
	a[tMax-1,jMax/2-1,0] = 1
	a = np.amax(a)-a

	results = []	
	for stepsToExtend in (0,tMax/3):
		moveToEnd = lambda (t,j,i): t<stepsToExtend and j==a.shape[1]-1 and i==a.shape[2]-1
		results.append(findConnections(a[tMax/3-stepsToExtend:,:,:],0, (a.shape[1]-1,a.shape[2]-1), moveToEnd))
	return results	
	
		#~ 
	#~ S = [(tMax/3, j) for j in range(jMax)]
	#~ S.extend( [(tMax*2/3,j) for j in range(jMax/2,jMax)] )
	#~ S.extend( [(t,jMax/2) for t in range(tMax/3,tMax*2/3)] )
	
	#~ results = []
	#~ for stepsToExtend in range(tMax/2):
		#~ moveToEnd = lambda (t,j,i): t<stepsToExtend and j==a.shape[1]-1 and i==a.shape[2]-1
		#~ results.append(findConnections(a[tMax/2-stepsToExtend:,:,:],0, (a.shape[1]-1,a.shape[2]-1), moveToEnd))
		#~ 
	#~ return results
	
	#moveToEnd = lambda (t,j,i): t<a.shape[0]/3 and j==a.shape[1]-1 and i==a.shape[2]-1
	#return findConnections(a,0,(a.shape[1]-1,a.shape[2]-1) ) #,moveToEnd)
	
	#J,I = np.meshgrid(np.linspace(0,4*math.pi,10), np.linspace(0,4*math.pi,10))
	#A = np.array([J]*10)
	#S = np.sin(A)
	##This seems to go in the wrong direction, try to swap axes
	#S = S.swapaxes(0,2)
	##S = S.swapaxes(1,2)
	#fract=1.1
	#moveToEnd = lambda (t,j,i): j*S.shape[0]*fract > S.shape[1]*t
	#return findConnections(S, 0, (S.shape[1]-1,S.shape[2]-1), moveToEnd) #This gives different results with different fracts. As it should be.	
	
	#a = np.zeros([2,3,3])
	#for t in range(2):
	# 	for j in range(3):
	#		a[t,j,0]=1
	#		a[t,j,2]=2
	#	a[t,2,1]=1
	#
	#return findConnections(a,0,(a.shape[1]-1,a.shape[2]-1),moveToEnd)
	
	
	
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
		d0 = f.getDiagrams(0)
		
		import heapq
		#Get the two cycles that are born first.
		d = heapq.nlargest(2, d0, key=lambda c: c.birth)
		
		values = values + [min([c.death for c in d])] #bec is the filtration level were the two oldest 0-cycles merge.
	end = time.time()
	print "It took", (end-start)/60, "minutes."
	return values

	
def getIntegratedFlux(year):
	"""Return arrays of same dimensions as got from metopen, where value is sum of wind speeds weighted by humidity.
	-return UQ: integral from 1000 to 200 of u*q dp
	-return VQ: integral from 1000 to 200 of v*q dp
	"""
	import inout, ha
	import numpy as np
	pl = var.pressure_levels
	#pl = [1000, 950] test
	
	U = {p : ha.loadU(p,year) for p in pl}
	V = {p : ha.loadV(p,year) for p in pl}
	Q = {p : ha.loadQ(p,year) for p in pl}
	
	shape = U[pl[0]].shape #all these should be equal
	#Calculating integral by trapezoidal rule
	UQ = np.zeros(shape)
	for i in range(len(pl)-1):
		UQ = UQ + 0.5*(U[pl[i]]*Q[pl[i]] + U[pl[i+1]]*Q[pl[i+1]]) * (pl[i]-pl[i+1])
	
	VQ = np.zeros(shape)
	for i in range(len(pl)-1):
		VQ = VQ + 0.5*(V[pl[i]]*Q[pl[i]] + V[pl[i+1]]*Q[pl[i+1]]) * (pl[i]-pl[i+1])
		
	return UQ, VQ


def compBecTime(month):
	"""Compare the bec to the time calculation as follows:
	For all cycles in time calculation, find at what time it hits Bergen (not completely well defined, but...)
	Compare the birth value of the cycle to the bec at this time.
	Only do for 2012 as for now
	"""
	import inout, var
	extendFract = 0.3
	d,f = findAllRiversInMonth(2012, month, extendFract=extendFract)
	bec = inout.load('bec2012')
	monthInterval = len(bec)/12
	#Find correspondence between time steps in vertices of filtered complex, and time index in bec
	#It is given by t_bec = t_time + minTime
	numStepsExtended = int(monthInterval*extendFract)
	minTime = monthInterval*(month-1) - numStepsExtended
	#Dictionary to store the time steps a cycle hits bergen.
	bergenTimes = {}
	#Go through all cycles, find their basis and see if it goes through Bergen.
	for cycle in d:
		bergenTimes[cycle] = set()
		for cell in f.getAllCellsInCycle(cycle):
			for k in range(1):
				vertex = cell.vertices[k]
				if type(vertex) == tuple: #The vertex might be the string 'extra_vertex'
					j = vertex[1]
					i = vertex[2]
					if var.x(i) == var.bergenLong and var.y(j) == var.bergenLat:
						bergenTimes[cycle].add(vertex[0]+minTime)
	#TODO: Make this do what I want.
	return bergenTimes
				
			
	
	
	
	
def compTimeBec(month):
	"""Compare the time calculation to the bec, by calculating persistent homology induced by the bec.
	"""
	import inout, cfc, plot, heapq
	import matplotlib.pyplot as plt
	import numpy as np
	#Load data
	ds = inout.load('Time_diagrams/2012/time_diagram_{:02d}_2012'.format(month))
	bec = inout.load('bec2012')
	timeInterval = len(bec)/12
	bec = bec[(month-1)*timeInterval:month*timeInterval]
	#Compute persistent homology given by bec
	f = cfc.makeCubicalComplex(np.array(bec), 't')
	f.computePersistentHomology(0)
	d = f.getDiagrams(0)
	#A bit more noise in the ds, so choose only the longest living ones, making the two the same length.
	n = len(d.cycles)
	persistence=lambda c:c.lifelength()
	ds.cycles = heapq.nlargest(n, ds.cycles, key=persistence)
	ds.cycles.sort(key=persistence)
	d.cycles.sort(key=persistence)
	#Plot
	fig = plt.figure()
	timeAx = fig.add_subplot(211)
	becAx = fig.add_subplot(212)
	plot.plotDiagrams(timeAx, ds, 'The tbec')
	plot.plotDiagrams(becAx, d, 'Induced by bec')
	#Get the two plots on the same scale
	timeBound = timeAx.get_xbound()
	becBound = becAx.get_xbound()
	lower = min( timeBound[0], becBound[0] )
	upper = max( timeBound[1], becBound[1] )
	timeAx.set_xbound(lower, upper)
	becAx.set_xbound(lower,upper)
	#Save the figure to file
	plt.savefig(inout.path+'time_vs_bec_diagrams_2012_{:02d}'.format(month), dpi=200, bbox_inches = 'tight')
	plt.clf()


