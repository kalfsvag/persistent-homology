"""
Contains the class filteredComplex

Author: Kristian Alfsvag
"""



class FilteredComplex:
	"""Class representing a filtered complex.
	 """
	 
	def __init__(self, cells, filtDir = 'b'):
		"""
		cells:		list of cells in filtration. Must be a cellular complex (i.e. closed under boundary operation)
		filtDir:	'b' for bottom up filtration (Usual order on filtration values)
					't' for top down filtration (Reversed order on filtration values)
					None for disregarding filtration values
		"""
		self.cells = cells
		self.maxDimComputed = -1
		self.diagrams = None
		self.filtDir = filtDir
		if filtDir == 'b':
			self.filt_cmp = filtFromBottom_cmp
		elif filtDir == 't':
			self.filt_cmp = filtFromTop_cmp
		else:
			self.filt_cmp = lambda c1,c0: c1.dim-c0.dim #Compare only dimension.
		
		
	def addCell(self, c, prioritize=False):
		"""Add cell to filtration. Marks the homology as not computed.
		
		c:		Cell to be added.
		prioritize:	If True, add cell to beginning of list of cells, instead of appending to the end. This will make this cell appear before another cell in filtration, if they otherwise appear at same filtration step. This should have no effect on the homology, but will have effect on the basis elements found.
		"""
		#Set index where c should be inserted:
		if prioritize:
			i = 0
		else:
			i = -1
		self.cells.insert(i,c)
		self.maxDimComputed = -1
		
		 
	def getCells(self):
		"""Get list of cells in filtration."""
		return self.cells
		
	
	def getDiagrams(self, dim):
		"""Return the persistence diagram of dimension dim (if already computed)
		dim:	an integer >= 0"""
		from persistencediagram import PersistenceDiagram
		if dim<=self.maxDimComputed:
			return PersistenceDiagram(self.diagrams[dim], filtDir=self.filtDir)
		else:
			raise Exception("Tried to get persistent homology of dim "+str(dim)+". This is not yet computed.")
			
			
	def getAllCellsInCycle(self, cycle):
		"""
		Return all cells in cycle.
		cycle:	must be a cycle in diagram gotten from self.computePersistentHomology
		"""
		if self.maxDimComputed<0:
			raise Exception("Homology not yet computed.")
		try:
			i = self.index[cycle.birthCell]
		except KeyError:
			raise Exception("Cycle not in this filtered complex.")
		cellsInCycle = []
		for j in self.basis[i]:
			#cells[j] should be added modulo 2.
			if self.cells[j] in cellsInCycle:
				cellsInCycle.remove(self.cells[j])
			else:
				cellsInCycle.append(self.cells[j])
		return cellsInCycle
		
		
	def computePersistentHomology(self, maxDim, trackBasis=False, countMarks=False, trackProgression=False):
		"""Compute the persistent homology (with coefficients in the field Z/2) of the filtered complex.
		Algorithm given in "Computing Persistent Homology" (Zomorodian, Carlsson, 2005)
		
		maxDim:				maximal dimension to compute homology in
		trackBasis:			if True, keep track of basis for homology. Basis will be stored in self.basis.
		countMarks:			if True, count the number of times a marked cell is added to a cycle. Stored in list self.markCount.
		trackProgression:	if True, visualize how far the computation has progressed
		"""
			
		from persistencecycle import PersistenceCycle
		print 'Computing persistent homology'
		#Sort the cells
		self.cells.sort(self.filt_cmp)
		
		#Dictionary giving index of each cell
		self.index = {self.cells[i]: i for i in xrange(len(self.cells))}
		
		#List telling if a cell with given index does not correspond to pivot element. (I.e. is a generator of the kernel of the boundary. (I.e. corresponds to a cycle)
		nonPivot = [False]*len(self.cells)
		
		#Initialize list to hold all the cycles
		cycles = []
		
		T = [None]*len(self.cells)
		#If i is the maximal index where cells[i] is in the boundary of cells[j] (after removing pivot rows), then T[i][0] = [indices of non-pivotal part of boundary of cells[j] ], and T[i][1]=j). Corresponds to column reduced column in boundary matrix.
		if trackBasis:
			self.basis = [ [j] for j in xrange(len(self.cells))] #List keeping track of change of basis. If cells[j] represents a cycle, then basis[j] is a list of indices of all cells in this cycle
		else:
			self.basis = None #In case persistent homology has been computed before, and complex altered in between
		if countMarks:
			self.markCount = [int(cell.mark) for cell in self.cells] #Will start with 1 for cells already marked, 0 otherwise.
		

		oldFract=0
		for j in xrange(len(self.cells)): #For every cell
			
			if trackProgression:
				#Make a visualization of how far the algorithm has come
				newFract=(10*(j+1))/len(self.cells)
				if newFract > oldFract:
					print str(newFract)+'/10'
					oldFract=newFract

			c = self.cells[j] #Cell we are now computing boundary of
			
			if c.dim>maxDim+1:
				continue  #Not needed to calculate homology up to maxDim
			#Compute the boundary
			d = [self.index[f] for f in c.boundary()]
			#Remove pivot rows:
			toRemove=[]
			for k in d:
				#If k is pivot, remove k.
				if not nonPivot[k]:
					toRemove.append(k)
			for k in toRemove:
				d.remove(k)

					
			while len(d) > 0: #If not empty
				#Find max index of cells in boundary
				i = max(d)
				#Do nothing if i is not top element in previous column.
				if T[i] == None:
					break 
				#If i is already the top element in a previous column, add that column to d (with Z/2-coefficients) to remove i from d. This is column reduction of boundary matrix.
				if trackBasis:
					self.basis[j].extend(self.basis[T[i][1]]) #Change of basis
				if countMarks:
					self.markCount[j] += self.markCount[T[i][1]] 
				#assert self.cells[j].dim == self.cells[T[i][1]].dim #debug
				for k in T[i][0]:
					if k in d: 
						d.remove(k) # (1+1) % 2 = 0
					else:
						d.append(k) # (1+0) % 2 = 1
					
			#After getting out of this loop: Either d is empty or T[i] is empty
						
			#Now the matrix is column reduced in part corresponding to columns 0,...,j, and we can see what role cells[j] play in the persistent homology.
			
			if len(d)>0:
				#In this case, c kills cycle generated by cells[i].
				i = max(d)
				T[i] = (d,j)
				if self.cells[i].filtrationValue != c.filtrationValue: #Only add cycles that aren't killed immideately
					cycles.append( PersistenceCycle(self.cells[i], c, filtDir=self.filtDir) )
					
			else:
				#In this case, the boundary is empty, so c is in kernel of boundary map.
				nonPivot[j] = True
		
		#Have now found all cycles that are killed. Must find all cycles of dim<=maxDim that live for ever.
		for j in xrange(len(self.cells)):
			if nonPivot[j] and T[j]==None and self.cells[j].dim <= maxDim:
				#cells[j] generate a cycle, but is never killed.
				cycles.append( PersistenceCycle(self.cells[j], filtDir = self.filtDir) )
		
		
		self.diagrams = [ [cycle for cycle in cycles if cycle.dimension() == dim ] for dim in xrange(maxDim+1) ]
		self.maxDimComputed=maxDim
	
	
	
	
	
	
	def isFilteredComplex(self):
		"""
		Check if this really is a filtered complex:
			for every cell, the boundary must be included in the complex, and must appear in filtration before the cell
		
		return:	boolean saying if this is a filtered complex
		
		Warning: Do not use this method on large complexes. Will use VERY long time
		"""
		for c in self.cells:
			for face in c.boundary():
				#Every face of c must be in the complex
				ok=False
				for c2 in self.cells:
					if face==c2:
						ok=True
						face=c2 #Now face also has the filtration value
						break
				if not ok:
					return ok, str(face)+' is not in filtration, but is face of '+str(c) #face was not in complex, return false
				
				#And every face must appear in filtration before c
				if self.filt_cmp(c,face)<0:
					return False, str(c)+' appears before face '+str(face)
		return True
				
		
def filtFromTop_cmp(cell1, cell2):
	"""Comparator on cells giving a filtration from highest to lowest filtration value.
	cell1 < cell2 if filtVal1 > filtVal2
	if they have same filtration value, sort by dimension: cell1 < cell2 if dim1 < dim2
	"""
	#First compare filtrationValue, then dimension if filtrationValue is equal
	return cmp( (cell2.filtrationValue, cell1.dim), (cell1.filtrationValue, cell2.dim) )

def filtFromBottom_cmp(cell1, cell2):
	"""Comparator on cells giving a filtration from lowest to highest filtration value.
	cell1 < cell2 if filtVal1 < filtVal2
	if they have same filtration value, sort by dimension
	"""
	#First compare filtrationValue, then dimension if filtrationValue is equal
	return cmp( (cell1.filtrationValue, cell1.dim), (cell2.filtrationValue, cell2.dim) )
	
