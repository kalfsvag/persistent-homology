"""
Module for saving and loading.

Author: Kristian Alfsvag"""

path = '/Home/stud3/kal045/Masteroppgave/data/'

def load(name):
	"""Load file in path with given name."""
	import pickle
	
	f = file(path + name,'r')
	v = pickle.load(f)
	f.close()
	return v
	
def save(v, name):
	"""Save object v to path with given name."""
	import pickle
	
	f=file(path+name,'w')
	pickle.dump(v,f)
	f.close()
	return

