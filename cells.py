"""Contains classes representing cells in filtration"""

import numpy as np
import collections


class CubicalCell:
	"""Class representing cubical cells.
	
	"""
	
	def __init__(self, vertices, filtrationValue=None, mark=False):
		"""Make a cubical cell with given vertices
		
		vertices: 			In general the vertices are ordered in an array of length 2 in all directions. I.e. shape(vertices)==[2,2,...,2]. Items in array must be hashable. In case dim == 0, this is just a hashable object. For example a tuple giving coordinates of point.
		filtrationValue:	A number giving the filtration value of the cell.
		"""
					
		
		if isinstance(vertices, np.ndarray):
			
			#Comment these away when all is working
			#assert all([l==2 for l in vertices.shape]), 'Vertice array has wrong shape: ' + str(vertices.shape)
			#assert all([isinstance(v, collections.Hashable) for v in vertices.flat]), 'Non-vertex object in vertices array: '+str(vertices)
			
			self.dim = vertices.ndim
			vertices.flags.writeable = False #Want this object to be immutable.
		elif isinstance(vertices, collections.Hashable): #If it is a zero-dimensional cube
			self.dim = 0
		else:
			raise Exception('Not valid input: '+str(vertices))			
		self.vertices=vertices
		self.filtrationValue = filtrationValue		
		self.mark=mark
	
	
	
	def __eq__(self, other):
		"""Two cells are equal if their arrays of vertices are equal. (So if you rotate a cube, you will get a different one.)
		
		Zero- and one-dimensional cubes can also be equal to simplices.
		
		Note that two cells can have different filtration values, but still be equal.
		"""
		
		if not isinstance(other, CubicalCell):
			if isinstance(other, Simplex):
				if self.dim==0 and other.dim==0:
					return (self.vertices,) == other.vertices # A point is a point
				if self.dim==1 and other.dim==1:
					return tuple(self.vertices) == other.vertices #A line is a line
			else: return False
		eq = self.vertices == other.vertices
		if not isinstance(eq, bool):
			#Assume that in this case, eq is an array of boolean values
			eq = eq.all()
		return eq
		
	def __hash__(self):
		if self.dim == 0:
			#In this case, vertices is just one vertex, which is hashable
			return hash((self.vertices,)) #Hashcode must correspond to hashcode of zero-dim 
		else:
			#In this case, vertices is an array of vertices. Return the hash of the tuple of all these.
			return hash(tuple(self.vertices.flat))
		
		
	def __ne__(self,other):
		return not self.__eq__(other)
		
	def __str__(self):
		return str(self.vertices)
	
	def __repr__(self):
		return 'Cubical cell: '+str(self.vertices)+', Filtration value: '+str(self.filtrationValue)
		
	def getDim(self):
		"""Get dimension of cell"""
		return self.dim
		
	def getVertices(self):
		"""Get vertices of cell.
		Case dim>0: Vertices of cell, sorted in numpy array of shape [2,2,...,2]. Dimension of cell is dimension of this array.
		Case dim=0:	Just one vertex, which can be any hashable object.
		"""
		return self.vertices
	
	def getFiltrationValue(self):
		"""Get filtration value of cell.
		"""
		return self.filtrationValue
	
	def isMarked(self):
		"""See whether cell is marked or not.
		"""
		return self.mark
	
		
	def boundary(self):
		"""Return the boundary of the cube. Orientation is ignored, so is only good for computations with Z/2-coefficients.
		
		A face of the cube has vertices given by fixing one index in the array self.vertices
		
		return: A list of dim-1 dimensional cubical cells, each corresponding to a face of this cube
		"""
		#A 0-cell has no boundary
		if self.dim == 0:
			return []
			
		b = []
		v = self.vertices
		for i in range(self.dim):
			#Add faces where i'th index is fixed
			v = v.swapaxes(0,i)
			b.extend([CubicalCell(f) for f in v])
		return b
		
	
				



class Simplex:
	"""Class representing simplices.
	"""

	
	def __init__(self, vertices, filtrationValue=None, mark=False):
		"""Make a simplex with given vertices
		
		vertices: 			a tuple of distinct, ordered vertices. A vertex must be a hashable object. 
		filtrationValue: 	will usually induce the sorting in a filtered complex
		mark:				boolean variable to be se to True if it is needed
		"""
					
		
		if type(vertices) == tuple:
			self.dim = len(vertices)-1
		elif isinstance(vertices, collections.Hashable): #If vertices is a hashable object (and not tuple), assume that it is a zero-dimensional simplex.
			self.dim = 0
			vertices = (vertices,)
		else:
			raise Exception('Not valid input: '+str(vertices))			
		self.vertices=vertices
		self.filtrationValue = filtrationValue
		
		self.mark=mark

		
	def __eq__(self, other):
		"""Two simplices are equal if they have the same vertices, in the same ordering.
		A zero simplex can be equal to a zero dimensional cube
		A one-simplex can be equal to a one dimensional cube
		
		Note that cells can have different filtration values, but still be equal.
		"""
		
		if not isinstance(other, Simplex):
			if isinstance(other, CubicalCell):
				if self.dim==0 and other.dim==0:
					return self.vertices == (other.vertices,) # A point is a point
				if self.dim==1 and other.dim==1:
					return self.vertices == tuple(other.vertices) #A line is a line
			else: return False
		
		return self.vertices == other.vertices
		
		
	def __hash__(self):
		return hash(self.vertices)
				
	def __ne__(self,other):
		return not self.__eq__(other)
		
	def __str__(self):
		return str(self.vertices)
	
	def __repr__(self):
		return 'Simplex: '+str(self.vertices)+', Filtration value: '+str(self.filtrationValue)
		
	def getDim(self):
		"""Get dimension of cell."""
		return self.dim
		
	def getVertices(self):
		"""Get vertices of simplex.
		return: Tuple of hashable objects of length dim-1.
		"""
		return self.vertices
	
	def getFiltrationValue(self):
		"""Get filtration value of cell.
		"""
		return self.filtrationValue
	
	def isMarked(self):
		"""See whether cell is marked or not.
		"""
		return self.mark
		
	def boundary(self):
		"""Return the boundary of the simplex. Orientation is ignored, so is only good for computations with Z/2-coefficients.
		
		return: A list of dim-1 dimensional simplices, each corresponding to a face of this simplex.
		"""
		#A 0-cell has no boundary
		if self.dim == 0:
			return []
			
		b = []
		for v in self.vertices:
			l = list(self.vertices)
			l.remove(v)
			b.append(Simplex(tuple(l)))
		return b
				
		
