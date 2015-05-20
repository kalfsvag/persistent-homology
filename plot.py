"""
Module for plotting values.

Author: Kristian Alfsvag"""


def plotVelocity(u,v,q=None):
	"""Plot velocity vectors.
	u --- 2D array of u components.
	v --- 2D array of v components.
	q --- 2D array of data (optional)
	"""
	from pylab import quiver
	assert u.shape == v.shape
	x = range(u.shape[1])
	y = range(u.shape[0])[-1::-1]
	if q == None:
		quiver(x,y,u,v)
	else:
		assert q.shape == u.shape
		quiver(x,y,u,v,q)
		
def plotBettiNumbers(ax,d):
	"""Plot betti numbers of the persistence diagram d.
	"""
	sing, rank = d.getBettiNumbers()
	#This should not be necessary with latest version, but for diagrams computed earlier...
	if sing[-1] in (float('inf'),-float('inf')):
		sing.pop()
		rank.pop()
	for i in range(len(sing)-1):
		ax.plot((sing[i],sing[i+1]),(rank[i],rank[i]),'b') #Straight horizontal line.
		ax.plot((sing[i+1],sing[i+1]),(rank[i],rank[i+1]),'b') #Straight vertical line.
	#The last line to infinity is left.
	inf = sing[-1] + (sing[-1] - sing[0])*0.1
	ax.plot((sing[-1],inf), (rank[-1],rank[-1]),'b')
	
	ax.set_xbound(lower=sing[0],upper=inf)



def plotContourAndWind(ax, x, y, tcw, u = None, v = None):
	#import matplotlib.pyplot as plt
	#fig = plt.figure()
	#ax = fig.add_subplot('111')
	ax.contour(x,y,tcw)
	if (u, v)!=(None,None):
		ax.quiver(x,y,u,v)
	import var
	ax.plot(var.bergenLong, var.bergenLat, 'rx')
	ax.annotate(xy=(var.bergenLong,var.bergenLat), s='Bergen')
	ax.set_xlabel('Longitude')
	ax.set_ylabel('Latitude')
	return ax

def plotTimeSeriesAroundInd(ax, data, ind=-1, indRadius=15, label=''):
	"""Plot values. Time interval is 6 hours.
	"""
	if ind>=0:
		start=max([0,ind-indRadius])
		stop=min([len(data), ind+indRadius])	
	else:
		start=stop=None
	cut = slice(start,stop)
	days = [float(i)/4 for i in range(len(data))]
	ax.plot(days[cut],data[cut],label=label)
	if ind>=0:
		ax.plot(days[ind], data[ind], 'rx')
	ax.set_xlabel('Day in year' )
	
	
def plotContourTcwBec(tcw,ind, what, UQ=None, VQ=None):
	"""Plot figures for the thesis.
	"""
	
	import ha, var, inout
	import matplotlib.pyplot as plt
	data = ha.sliceArray(tcw, minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong)
	u = None
	v = None
	uvStr=''
	if (UQ, VQ)!=(None,None):
		u = ha.sliceArray(UQ, minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong)
		u = u[ind]
		v = ha.sliceArray(VQ, minLat=var.minLat, maxLat=var.maxLat, minLong=var.minLong, maxLong=var.maxLong)
		v = v[ind]
		uvStr='UV'
	x = [ha.longFromInd(i, minLong=var.minLong) for i in range(data.shape[2]) ]
	y = [ha.latFromInd(j, maxLat=var.maxLat) for j in range(data.shape[1]) ]
	tcwBergen=[data[t,2*var.j(var.bergenLat),2*var.i(var.bergenLong)] for t in range(1464)]
	bec2012=inout.load('bec2012')
	fig = plt.figure()
	cont =  plt.subplot2grid((3,1), (0,0), rowspan=2)
	grph =  plt.subplot2grid((3,1), (2,0) )
	plotContourAndWind(cont,x,y,data[ind],u,v)
	plotTimeSeriesAroundInd(grph, tcwBergen, ind, label='tcw Bergen' )
	plotTimeSeriesAroundInd(grph, bec2012, ind, label='bec')
	plt.legend()
	saveName=inout.path+'Contour_plots/High_'+what+uvStr+'/2012_'+str(ind)
	plt.savefig(saveName, dpi=200, bbox_inches='tight')
	
	plt.close()


		

def plotBecMonthly(year,monthInd):
	"""Plot the bec of a given year and month and save it.
	monthInd:	index of month (i.e. month in range(12))"""
	
	import matplotlib.pyplot as plt
	import inout
	
	months=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dec']
	
	b = inout.load(str(year)+'/bec'+str(year))
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	ax.plot(b[len(b)/12*monthInd: len(b)/12*(monthInd+1)])
	plt.title('bec '+months[monthInd]+' '+str(year)) 
	
	saveName = '/home/kal045/Masteroppgave/data/'+str(year)+'/bec'+str(year)+'_'+str(monthInd+1)
	
	plt.savefig(saveName, dpi=200,bbox_inches='tight')
	plt.clf()
	
#def compTwoDiagrams(ax, d1, d2, filtDir='b', title='')
	


def plotDiagrams(ax, d, title=''):
	"""
	Plot a persistence diagram.
	
	d:			diagram to be plotted
	"""

	import matplotlib.pyplot as plt
	if d.filtDir == 'b':
		inf = float('inf')
	elif d.filtDir == 't':
		inf = -float('inf')
	
	#Need to find max and min of non-infinite births and deaths in order to plot
	values = []
	for cycle in d.cycles:
		if cycle.birth != -inf:
			values.append(cycle.birth)
		if cycle.death != inf:
			values.append(cycle.death)
	maxVal = max(values)
	minVal = min(values)
	if d.filtDir=='b':
		start = minVal
		end = maxVal + (maxVal-minVal)*0.1 #move end a bit to the right to discern cycles lasting for ever and cycles lasting to maxVal
		infAdd = (maxVal-minVal)*10 #Should plot cycles lasting for ever for more than just in the window, so window can be expanded a little without problem
	if d.filtDir=='t':
		start = maxVal
		end = minVal - (maxVal-minVal)*0.1 #move end a bit to the left to discern cycles lasting for ever and cycles lasting to maxVal
		infAdd = -(maxVal-minVal)*10 #Should plot cycles lasting for ever for more than just in the window, so window can be expanded a little without problem
	
	#~ fig = plt.figure()
	#~ ax = fig.add_subplot(111)
	cycleNumber = 0
	for cycle in d.cycles:
		birth = cycle.birth
		if birth == -inf:
			birth = start - infAdd
		death = cycle.death
		if death == inf:
			death = end + infAdd
		ax.plot([cycle.birth, death], [cycleNumber, cycleNumber], 'b')
		cycleNumber += 1
	ax.axis([start, end, -len(d.cycles)/10, len(d.cycles)+len(d.cycles)/10])
	ax.set_title(title)
	
	#~ if save:
		#~ plt.savefig(saveName, dpi=200,bbox_inches='tight')
		#~ plt.clf()
	#~ else:
		#~ fig.show()
	return
	
def plotRiverDiagram(year, monthInd, (endVal, startVal)=(None,None) ):
	"""Plot a persistence diagram and save it.
	year
	monthInd:	index of month, is in range(12)
	endVal: 	use this if I want same scales on several plots
	startVal:	ditto
	"""
	import matplotlib.pyplot as plt
	import inout, ha
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Aug', 'Sep', 'Okt', 'Nov', 'Dec']
	month = months[monthInd]
	d = inout.load(month+'/riverdiagram'+str(year)+'%02d'%(monthInd+1) )
	#Find total precip in Bergen this month
	bp = ha.loadPrecip(year)
	l = len(bp)
	precip = [x for x in bp[len(bp)/12*monthInd:len(bp)/12*(monthInd+1)]] #precip in month
	#Currently only looking at rivers hitting in the last two thirds of month, so...
	precip = precip[len(precip)/3:-1]
	precip = sum([x for x in precip if x>0]) #Data may contain -1 (no rain) and -99999 (no data)
	
	#Copypasta from plotDiagrams:
	inf = float('inf')	
	
	if (endVal,startVal) == (None,None):
		#Need to find max and min of non-infinite births and deaths in order to plot
		values = [cycle.birth for cycle in d if cycle.birth != inf]+[cycle.death for cycle in d if cycle.death != inf]
		maxVal = max(values)
		minVal = min(values)
		#filtDir='t'
		start = maxVal
		end = minVal - (maxVal-minVal)*0.1 #move end a bit to the left to discern cycles
	else:
		start = startVal
		end = endVal
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	cycleNumber = 0
	for cycle in d:
		birth = cycle.birth
		if birth == inf:
			birth = start
		death = cycle.death
		if death == inf:
			death = end
		ax.plot([cycle.birth, death], [cycleNumber, cycleNumber], 'b')
		cycleNumber += 1
	ax.axis([start, end, -len(d)/10, len(d)+len(d)/10])
	plt.title('Atmospheric rivers '+month+' '+str(year)+'\nTotal precip this month: '+str(precip))
	
	saveName = inout.path+month+'/'+str(year)+'_'+'%02d'%(monthInd+1)+'_rivers'
	plt.savefig(saveName, dpi=200,bbox_inches='tight')
	plt.clf()
	
def plotDiagramsAdv(d, a, n=0, x=lambda i:i, y=lambda j:-j, filtDir='b'): #, plotPair = False):
	"""Plots a persistence diagram together with some ties to the function it is derived from.
	parameters:
		d:	Persistence diagram
		a:	The grid function that induces d (as a numpy array)
		n:	The n longest living cycles will get some focus. 
		"""
	import matplotlib.pyplot as plt
	import heapq
	import numpy as np
	import var

	
	#Variant med matrise. Antar at d er faatt fra matrisen vha getPersistenceFromArray og getDiagrams. Plotter to figurer: en som plotter persistensdiagrammet, og en som tar de ti stoerste og proever aa vise hvor i verden de er.
	
	if filtDir=='b':
		M = np.max
		m = np.min
	elif filtDir=='t':
		M = np.min
		m = np.max
	
	start = m(a)
	end = M(a)*11/10

	fig = plt.figure()
	dia = fig.add_subplot(211)
	dia.set_xlabel('Filtration value')
	mp = fig.add_subplot(212)
	mp.set_xlabel('Longitude')
	mp.set_ylabel('Latitude')
	
	mp.contour([x(i) for i in range(a.shape[1])],[y(j) for j in range(a.shape[0])], a)

	dFocus = heapq.nlargest(n, d, key=lambda c:c.lifelength())	

	cycleNumber = 0
	label = 0
	for cycle in d:
		death = cycle.death
		if death == float('inf'):
			death = end
		dia.plot([cycle.birth, death], [cycleNumber, cycleNumber], 'b')
		
		#Dersom vi har en av de ti stoerste:
		if cycle in dFocus:
			#Finner et hjoerne i simplekset, og koordinatene til dette:
			(j,i) = cycle.birthCell.vertices[0]

			#Setter merkelapp paa linjen.
			dia.annotate(label, xy = (death, cycleNumber))
			mp.plot(x(i), y(j), 'rx')
			mp.annotate(label, (x(i),y(j)))
			
			label += 1

		cycleNumber += 1

	dia.axis([start, end, -len(d)/10, len(d)+len(d)/10])
	mp.axis([x(0), x(a.shape[1]-1), y(a.shape[0]-1), y(0)])


	plt.show()
	

	
