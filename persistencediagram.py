"""
Module containing the class persistenceDiagram

Author: Kristian Alfsvag
"""

class PersistenceDiagram():
	"""
	
	"""
	
	def __init__(self, cycles, filtDir='t'):
		"""
		cycles:		List of cycles
		filtDir:	'b': Ordinary order on filtration values
					't': Reverse order on filtration values					
		"""
		self.cycles = cycles
		self.filtDir = filtDir
		#~ if filtDir == 'b':
			#~ self.reverse=False
			#~ self.comp = cmp
		#~ elif filtDir == 't':
			#~ self.reverse=True
			#~ self.comp = lambda x,y: cmp(y,x)
		#~ else:
		if not filtDir in ('t', 'b'):
			raise Exception('Invalid filtration direction: '+str(self.filtDir))

	
	def getBettiNumbers(self):
		"""Find the betti numbers of the persistent homology corresponding to the diagram.
		
		return: sing: All the singular filtration values, i.e. the filtration values where the persistent homology changes
				rank: rank[sing[i]] is the rank of the persistent homology in the filtration interval [sing[i], sing[i+1])
		"""
		import numpy as np
		if self.filtDir == 'b':
			reverse=False
			comp = cmp
			inf = float('inf')
		elif self.filtDir == 't':
			reverse=True
			comp = lambda x,y: cmp(y,x)
			inf = -float('inf')
			
		#The singular values are all filtration values where the homology changes.
		sing = []
		for cycle in self.cycles:
			sing.append(cycle.birth)
			if cycle.death != inf:
				sing.append(cycle.death)
		#Get rid of duplicates
		sing = list(set(sing))
		sing.sort(reverse=reverse)
			
		rank = [ [comp(val,cycle.birth)>=0 and comp(val,cycle.death)<0 for cycle in self.cycles].count(True) for val in sing ]
		
				
		
		return sing, rank		
		
		
	def __repr__(self):
		return repr(self.cycles)

	def __str__(self):
		return str(self.cycles)
			
