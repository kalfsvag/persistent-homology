"""
Module for constructing filtered complexes

Author: Kristian Alfsvag
"""

def make1SkeletonFrom2Darray(a, filtDir='b'):
	"""Make filtered cubical complex from array
	Complex har vertices given by coordinates of a, filtration value at vertice v is a[v]
	Cells in complex are all vertices, and lines between these.
	
	a:			a 2D numpy array
	filtDir:	string giving direction of filtration.
				't' gives filtration from top
				'b' gives filtration from bottom
	"""
	import numpy as np
	import itertools
	from cells import Simplex
	from filteredcomplex import FilteredComplex
	
	#Make filtration go in right direction
	if filtDir == 't':
		m = min #Filtration value of cell should be minimum of filtration value on boundary
	elif filtDir == 'b':
		m = max #Filtration value of cell should be maximum of filtration value on boundary
	else:
		raise Exception('Invalid filtration direction: '+str(filtDir))
	
	#Start by adding all vertices to list of cells
	cells = [Simplex((v,), a[v]) for v in itertools.product(*[range(d) for d in a.shape]) ]
	
	#Add all 1-cells:
	coords = itertools.product(*[range(d) for d in a.shape])
	for c in coords:
		#Try to go one step in each direction, and if possible add line in that direction
		if c[0] < a.shape[0]-1:
			c1 = (c[0]+1,c[1])
			vertices = (c,c1)
			cell = Simplex(vertices, m(a[c],a[c1])) 
			cells.append(cell)
		if c[1] < a.shape[1]-1:
			c1 = (c[0],c[1]+1)
			vertices = (c,c1)
			cell = Simplex(vertices, m(a[c],a[c1]) )
			cells.append(cell)

	return FilteredComplex(cells, filtDir)
		

def make2SkeletonFrom3Darray(a, filtDir='t'):
	"""Make filtered cubical complex from array
	Complex har vertices given by coordinates of a, filtration value at vertice v is a[v]
	Cells in complex are all vertices, lines, and squares between these.
	
	a:			a 3D numpy array
	filtDir:	string giving direction of filtration.
				't' gives filtration from top
				'b' gives filtration from bottom
	"""
	import numpy as np
	import itertools
	from cells import CubicalCell
	from filteredcomplex import FilteredComplex
		
	#Make filtration go in right direction
	if filtDir == 't':
		m = min #Filtration value of cell should be minimum of filtration value on boundary
	elif filtDir == 'b':
		m = max #Filtration value of cell should be maximum of filtration value on boundary
	else:
		raise Exception('Invalid filtration direction: '+str(filtDir))
	
	#Start by adding all vertices to list of cells
	cells = [CubicalCell(v, a.item(v)) for v in itertools.product(*[range(d) for d in a.shape]) ]
	
	#Add all 1-cells:
	coords = itertools.product(*[range(d) for d in a.shape])
	for c in coords:
		#Try to go one step in each direction, and if possible add line in that direction
		if c[0] < a.shape[0]-1:
			c1 = (c[0]+1,c[1],c[2])
			vertices = np.ndarray(2,object)
			vertices[0] = c
			vertices[1] = c1
			cell = CubicalCell(vertices, m(a.item(c),a.item(c1)) )
			cells.append(cell)
		if c[1] < a.shape[1]-1:
			c1 = (c[0],c[1]+1,c[2])
			vertices = np.ndarray(2,object)
			vertices[0] = c
			vertices[1] = c1
			cell = CubicalCell(vertices, m(a.item(c),a.item(c1)) )
			cells.append(cell)
		if c[2] < a.shape[2]-1:
			c1 = (c[0],c[1],c[2]+1)
			vertices = np.ndarray(2,object)
			vertices[0] = c
			vertices[1] = c1
			cell = CubicalCell(vertices, m(a.item(c),a.item(c1)) )
			cells.append(cell)
	
			
	#Add all 2-cells:
	coords = itertools.product(*[range(d) for d in a.shape])
	for c in coords:
		#Other vertices in squares with lower corner c
		c100 = (c[0]+1,c[1],c[2])
		c010 = (c[0],c[1]+1,c[2])
		c001 = (c[0],c[1],c[2]+1)
		c110 = (c[0]+1,c[1]+1,c[2])
		c101 = (c[0]+1,c[1],c[2]+1)
		c011 = (c[0],c[1]+1,c[2]+1)
		
		#Square of vertices with 3rd component 0
		if c[0] < a.shape[0]-1 and c[1] < a.shape[1]-1:
			#Make array of right shape
			vertices = np.ndarray([2,2],object)
			vertices[0,0]=c
			vertices[1,0]=c100
			vertices[0,1]=c010
			vertices[1,1]=c110
			filtrationValue = m(a.item(c), a.item(c100), a.item(c010), a.item(c110))
			cells.append(CubicalCell(vertices,filtrationValue))
		
		#Square of vertices with 2nd component 0
		if c[0] < a.shape[0]-1 and c[2] < a.shape[2]-1:
			#Make array of right shape
			vertices = np.ndarray([2,2],object)
			vertices[0,0]=c
			vertices[1,0]=c100
			vertices[0,1]=c001
			vertices[1,1]=c101
			filtrationValue = m(a.item(c), a.item(c100), a.item(c001), a.item(c101))
			cells.append(CubicalCell(vertices,filtrationValue))
			
		#Square of vertices with 1st component 0
		if c[1] < a.shape[1]-1 and c[2] < a.shape[2]-1:
			#Make array of right shape
			vertices = np.ndarray([2,2],object)
			vertices[0,0]=c
			vertices[1,0]=c010
			vertices[0,1]=c001
			vertices[1,1]=c011
			filtrationValue = m(a.item(c), a.item(c010), a.item(c001), a.item(c011))
			cells.append(CubicalCell(vertices,filtrationValue))
			
	return FilteredComplex(cells, filtDir)
			


def makeCubicalComplex(a, filtDir='t'):
	"""Make a filtered cubical complex from array.
	Complex has vertices given by coordinates of a, filtration value at vertice v is a[v].
	Filtration value at higher dimensional cell is max value at vertices.
	
	a:			a numpy array
	filtDir:	string giving direction of filtration.
				't' gives filtration from top
				'b' gives filtration from bottom
	
	assumes dimension > 0"""
	
	import numpy as np
	import itertools
	from cells import CubicalCell
	from filteredcomplex import FilteredComplex
	
	#dimension of array (and also of complex)
	dim = a.ndim
	#iterator over all coordinates of a that are lower corners of a cube (those not at the upper boundaries)
	coords = itertools.product(*[range(d-1) for d in a.shape])
	#Set of cells in filtration. (Is a set to avoid same cell added twice)
	cells = set()
	
	#Make filtration go in right direction
	if filtDir == 't':
		m = min #Filtration value of cell should be minimum of filtration value on boundary
	elif filtDir == 'b':
		m = max #Filtration value of cell should be maximum of filtration value on boundary
	else:
		raise Exception('Invalid filtration direction: '+str(filtDir))
	
	for c in coords:
		vertices = np.ndarray([2]*dim, dtype=object)
		#iterator of coordinates of cube with lowest corner in c (relative to c). Is also coordinates of vertices
		locCoords = itertools.product(*[range(2)]*dim)
		#Want to add tuples as vector, so make them into np arrays
		ca = np.array(c)
		for lc in locCoords:
			x = tuple( ca + np.array(lc) )
			#Set the value of vertices array to corresponding vertex
			vertices.itemset(lc, x)
		
		#Highest dimensional cube with lowest corner in c:
		topCell = CubicalCell(vertices)
		#Recursively add topCell and its boundary, will also set filtration value
		addCell(topCell, cells, a, m)
	#Now all cells are added, and have correct filtration value
	return FilteredComplex(list(cells), filtDir)
			
			
def addCell(cell, cell_set, a,  m):
	"""Recursively add cell and its boundary, computing filtration value from boundary.
	cell:		cell to be added
	cell_set: 	set cell will be added to
	a:			array giving filtration value at 0-cells
	m:			min if filtration is from top
				max if filtration is from bottom
	"""
	if cell.dim == 0:
		#Filtration value for 0-dim cell is given by a[vertices]
		filtrationValue = a[cell.vertices]
	else:
		#If higher dimensional cell, add boundary to set and compute filtration value from boundary.
		filtrationValue = m([addCell(f, cell_set, a, m) for f in cell.boundary()])
	cell.filtrationValue = filtrationValue
	cell_set.add(cell)
	return filtrationValue

