"""
An example illustrating how to compute homology.
"""
from filteredcomplex import FilteredComplex
from cells import Simplex

s = Simplex( (0,1,2) , 2) #A filled triangle with three vertices: 0,1,2 and filtration value 2.
simplices = [s]

#The following makes sure that a list of simplices will get all the boundaries added, so that it actually corresponds to a complex.
#All simplices have its dimension as filtration value.
i=0
while i < len(simplices):
	for simplex in simplices[i].boundary():
		simplex.filtrationValue = simplex.getDim()
		simplices.append(simplex)
	i+=1
simplices = list(set(simplices)) #Just a simple way to get rid of duplicates.

#Construct complex and compute homology up to dimension maxDim.
f = FilteredComplex(simplices)
maxDim = 1
f.computePersistentHomology(maxDim)

#Print the homology.
for d in range(maxDim+1):
	print 'Homology of dimension %d: \n %s.\n' %(d, f.getDiagrams(d) )




