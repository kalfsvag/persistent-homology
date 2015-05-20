"""
Module containing the class PersistenceCycle

Author: Kristian Alfsvag
"""

class PersistenceCycle:
	"""Class representing a cycle in persistent homology.

	representative:	cell generating cycle
	killingSimplex:	cell killing cycle
	birth: 			birth of cycle
	death: 			death of cycle
	"""
	
	def __eq__(self, other):
		"""Two cycles are equal if they have same dimension and persistence interval.
		"""
		if not isinstance(other, PersistenceCycle):
			return False
		else:
			return self.dimension() == other.dimension() and self.birth == other.birth and self.death == other.death
			
	def __hash__(self):
		"""Considers dimension, birth and death
		"""
		return hash( (self.dimension(), self.birth, self.death ) )
		
	
	def __init__(self, birthCell, deathCell = None, filtDir='b', mark=False): #, cells=None):
		"""
		birthCell:	first cell in cycle
		deathCell:	cell killing cycle
		cells:		list of all cells in cycle, if given.
		mark:		boolean value saying if cycle is marked.
		"""
		self.birthCell = birthCell
		self.deathCell = deathCell
		self.birth = birthCell.filtrationValue
		self.mark = mark
		#Find if we should use + or - infinity
		if filtDir == 't':
			inf = -float('inf')
		else:
			inf = float('inf')		
		
		if deathCell != None:
			self.death = deathCell.filtrationValue
		else:
			self.death = inf
			
	def __repr__(self):
		return 'Persistencecycle <dimension: '+repr(self.dimension())+', birth: '+repr(self.birth) + ', death: '+repr(self.death)+'>'

	def __str__(self):
		return 'dimension: '+repr(self.dimension())+', birth: '+repr(self.birth) + ', death: '+repr(self.death)

	def dimension(self):
		return self.birthCell.dim

	def lifelength(self):
		return abs(self.death - self.birth)
