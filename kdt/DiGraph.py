import numpy as np
import scipy as sc
import scipy.sparse as sp
import pyCombBLAS as pcb
import Graph as gr

class DiGraph(gr.Graph):

	def toBool(self):
		if isinstance(self.spm, pcb.pySpParMat):
			self.spm = pcb.pySpParMatBool(self.spm)

	#print "in DiGraph"

	# NOTE:  for any vertex, out-edges are in the column and in-edges
	#	are in the row
	def __init__(self,*args):
		if len(args) == 0:
			self.spm = pcb.pySpParMat()
		elif len(args) == 1:	# no longer used
			[arg] = args
			if arg < 0:
				self.spm = pcb.pySpParMat()
			else:
				raise NotImplementedError, '1-argument case only accepts negative value'
		elif len(args) == 4:
			[i,j,v,nv] = args
			if type(v) == int or type(v) == long or type(v) == float:
				v = ParVec.broadcast(len(i),v)
			self.spm = pcb.pySpParMat(nv,nv,i.dpv,j.dpv,v.dpv)
		elif len(args) == 5:
			[i,j,v,nv1,nv2] = args
			if type(v) == int or type(v) == long or type(v) == float:
				v = ParVec.broadcast(len(i),v)
			if i.max() > nv1-1:
				raise KeyError, 'at least one first index greater than #vertices'
			if j.max() > nv2-1:
				raise KeyError, 'at least one second index greater than #vertices'
			self.spm = pcb.pySpParMat(nv1,nv2,i.dpv,j.dpv,v.dpv)
		else:
			raise NotImplementedError, "only 1, 4, and 5 argument cases supported"

	def __add__(self, other):
		if type(other) == int or type(other) == long or type(other) == float:
			raise NotImplementedError
		if self.nvert() != other.nvert():
			raise IndexError, 'Graphs must have equal numbers of vertices'
		elif isinstance(other, DiGraph):
			ret = self.copy()
			ret.spm += other.spm
			#ret.spm = pcb.EWiseApply(self.spm, other.spm, pcb.plus());  # only adds if both mats have nonnull elems!!
		return ret

	def __div__(self, other):
		if type(other) == int or type(other) == long or type(other) == float:
			ret = self.copy()
			ret.spm.Apply(pcb.bind2nd(pcb.divides(),other))
		elif self.nvert() != other.nvert():
			raise IndexError, 'Graphs must have equal numbers of vertices'
		elif isinstance(other,DiGraph):
			ret = self.copy()
			ret.spm = pcb.EWiseApply(self.spm, other.spm, pcb.divides())
		else:
			raise NotImplementedError
		return ret

	def __getitem__(self, key):
		#ToDo:  accept slices for key0/key1 besides ParVecs
		if type(key)==tuple:
			if len(key)==1:
				[key0] = key; key1 = -1
			elif len(key)==2:
				[key0, key1] = key
			else:
				raise KeyError, 'Too many indices'
		else:
			key0 = key;  key1 = key
		if type(key0) == int or type(key0) == long or type(key0) == float:
			tmp = ParVec(1)
			tmp[0] = key0
			key0 = tmp
		if type(key1) == int or type(key0) == long or type(key0) == float:
			tmp = ParVec(1)
			tmp[0] = key1
			key1 = tmp
		if type(key0)==slice and key0==slice(None,None,None):
			key0mn = 0; 
			key0tmp = self.nvert()
			if type(key0tmp) == tuple:
				key0mx = key0tmp[0] - 1
			else:
				key0mx = key0tmp - 1
		else:
			key0mn = int(key0.min()); key0mx = int(key0.max())
			if len(key0)!=(key0mx-key0mn+1) or not (key0==ParVec.range(key0mn,key0mx+1)).all():
				raise KeyError, 'Vector first index not a range'
		if type(key1)==slice and key1==slice(None,None,None):
			key1mn = 0 
			key1tmp = self.nvert()
			if type(key1tmp) == tuple:
				key1mx = key1tmp[1] - 1
			else:
				key1mx = key1tmp - 1
		else:
			key1mn = int(key1.min()); key1mx = int(key1.max())
			if len(key1)!=(key1mx-key1mn+1) or not (key1==ParVec.range(key1mn,key1mx+1)).all():
				raise KeyError, 'Vector second index not a range'
		[i, j, v] = self.toParVec()
		sel = ((i >= key0mn) & (i <= key0mx) & (j >= key1mn) & (j <= key1mx)).findInds()
		newi = i[sel] - key0mn
		newj = j[sel] - key1mn
		newv = v[sel]
		ret = DiGraph(newi, newj, newv, key0mx-key0mn+1, key1mx-key1mn+1)
		return ret

	def __iadd__(self, other):
		if type(other) == int or type(other) == long or type(other) == float:
			raise NotImplementedError
		if self.nvert() != other.nvert():
			raise IndexError, 'Graphs must have equal numbers of vertices'
		elif isinstance(other, DiGraph):
			#dead tmp = pcb.EWiseApply(self.spm, other.spm, pcb.plus())
			self.spm += other.spm
		return self

	def __imul__(self, other):
		if type(other) == int or type(other) == long or type(other) == float:
			self.spm.Apply(pcb.bind2nd(pcb.multiplies(),other))
		elif isinstance(other,DiGraph):
			self.spm = pcb.EWiseApply(self.spm,other.spm, pcb.multiplies())
		else:
			raise NotImplementedError
		return self

	def __mul__(self, other):
		if type(other) == int or type(other) == long or type(other) == float:
			ret = self.copy()
			ret.spm.Apply(pcb.bind2nd(pcb.multiplies(),other))
		elif self.nvert() != other.nvert():
			raise IndexError, 'Graphs must have equal numbers of vertices'
		elif isinstance(other,DiGraph):
			ret = self.copy()
			ret.spm = pcb.EWiseApply(self.spm,other.spm, pcb.multiplies())
		else:
			raise NotImplementedError
		return ret

	#ToDo:  put in method to modify _REPR_MAX
	_REPR_MAX = 100
	def __repr__(self):
		if self.nvert() == 0:
			return 'Null DiGraph object'
		if self.nvert()==1:
			[i, j, v] = self.toParVec()
			if len(v) > 0:
				print "%d %f" % (v[0], v[0])
			else:
				print "%d %f" % (0, 0.0)
		else:
			[i, j, v] = self.toParVec()
			if len(i) < self._REPR_MAX:
				print i,j,v
		return ' '

	def _SpMM(self, other):
		selfnv = self.nvert()
		if type(selfnv) == tuple:
			[selfnv1, selfnv2] = selfnv
		else:
			selfnv1 = selfnv; selfnv2 = selfnv
		othernv = other.nvert()
		if type(othernv) == tuple:
			[othernv1, othernv2] = othernv
		else:
			othernv1 = othernv; othernv2 = othernv
		if selfnv2 != othernv1:
			raise ValueError, '#in-vertices of first graph not equal to #out-vertices of the second graph '
		ret = DiGraph()
		ret.spm = self.spm.SpMM(other.spm)
		return ret

	def bool(self):
		#ToDo:  change for real Boolean matrices
		return DiGraph.ones(self)

	def copy(self):
		ret = DiGraph()
		ret.spm = self.spm.copy()
		return ret
		
	def degree(self, dir=gr.Out):
		if dir == gr.InOut:
			#ToDo:  can't do InOut if nonsquare graph
			tmp1 = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.plus(), pcb.ifthenelse(pcb.bind2nd(pcb.not_equal_to(), 0), pcb.set(1), pcb.set(0)))
			tmp2 = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.plus(), pcb.ifthenelse(pcb.bind2nd(pcb.not_equal_to(), 0), pcb.set(1), pcb.set(0)))
			return ParVec.toParVec(tmp1+tmp2)
		elif dir == gr.In:
			ret = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.plus(), pcb.ifthenelse(pcb.bind2nd(pcb.not_equal_to(), 0), pcb.set(1), pcb.set(0)))
			return ParVec.toParVec(ret)
		elif dir == gr.Out:
			ret = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.plus(), pcb.ifthenelse(pcb.bind2nd(pcb.not_equal_to(), 0), pcb.set(1), pcb.set(0)))
			return ParVec.toParVec(ret)
		else:
			raise KeyError, 'Invalid edge direction'

	# in-place, so no return value
	def removeSelfLoops(self):
		self.spm.removeSelfLoops()
		return

	@staticmethod
	def fullyConnected(n,m=None):
		if m == None:
			m = n
		i = (ParVec.range(n*m) % n).floor()
		j = (ParVec.range(n*m) / n).floor()
		v = ParVec.ones(n*m)
		ret = DiGraph(i,j,v,n,m)
		return ret

	def genGraph500Edges(self, scale):
		elapsedTime = self.spm.GenGraph500Edges(scale)
	 	return elapsedTime

	@staticmethod
	def load(fname):
		#FIX:  crashes if any out-of-bound indices in file; easy to
		#      fall into with file being 1-based and Py being 0-based
		ret = DiGraph()
		ret.spm = pcb.pySpParMat()
		ret.spm.load(fname)
		return ret

	def max(self, dir=gr.InOut):
		#ToDo:  is default to InOut best?
		if dir == gr.InOut:
			tmp1 = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.max())
			tmp2 = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.max())
			return ParVec.toParVec(tmp1+tmp2)
		elif dir == gr.In:
			ret = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.max())
			return ParVec.toParVec(ret)
		elif dir == gr.Out:
			ret = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.max())
			return ParVec.toParVec(ret)
		else:
			raise KeyError, 'Invalid edge direction'

	def min(self, dir=gr.InOut):
		#ToDo:  is default to InOut best?
		if dir == gr.InOut:
			tmp1 = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.min())
			tmp2 = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.min())
			return ParVec.toParVec(tmp1+tmp2)
		elif dir == gr.In:
			ret = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.min())
			return ParVec.toParVec(ret)
		elif dir == gr.Out:
			ret = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.min())
			return ParVec.toParVec(ret)

	def mulNot(self, other):
		if self.nvert() != other.nvert():
			raise IndexError, 'Graphs must have equal numbers of vertices'
		else:
			ret = DiGraph()
			ret.spm = pcb.EWiseApply(self.spm, other.spm, pcb.multiplies(), True)
		return ret

	#FIX:  good idea to have this return an int or a tuple?
	def nvert(self):
		nrow = self.spm.getnrow()
		ncol = self.spm.getncol()
		if nrow==ncol:
			ret = nrow
		else:
			ret = (nrow, ncol)
		return ret

	#in-place, so no return value
	def ones(self):
		self.spm.Apply(pcb.set(1))
		return

	#in-place, so no return value
	def reverseEdges(self):
		self.spm.Transpose()

	def save(self, fname):
		self.spm.save(fname)
		return

	#in-place, so no return value
	def scale(self, other, dir=gr.Out):
		#Note:  have to compare against gr.SpParVec, not (local) SpParVec
		if not isinstance(other,gr.SpParVec):
			raise KeyError, 'Invalid type for scale vector'
		selfnv = self.nvert()
		if type(selfnv) == tuple:
			[selfnv1, selfnv2] = selfnv
		else:
			selfnv1 = selfnv; selfnv2 = selfnv
		if dir == gr.Out:
			if selfnv2 != len(other):
				raise IndexError, 'graph.nvert()[1] != len(scale)'
			self.spm.ColWiseApply(other.spv, pcb.multiplies())
		elif dir == gr.In:
			if selfnv1 != len(other):
				raise IndexError, 'graph.nvert()[1] != len(scale)'
			self.T()
			self.spm.ColWiseApply(other.spv,pcb.multiplies())
			self.T()
		else:
			raise KeyError, 'Invalid edge direction'
		return

	#in-place, so no return value
	def set(self, value):
		self.spm.Apply(pcb.set(value))
		return

	def subgraph(self, *args):
		if len(args) == 1:
			[ndx1] = args
			ret = self[ndx1, ndx1]
		elif len(args) == 2:
			[ndx1, ndx2] = args
			ret = self[ndx1, ndx2]
		else:
			raise IndexError, 'Too many indices'
		return ret

	def sum(self, dir=gr.Out):
		if dir == gr.InOut:
			tmp1 = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.plus(), pcb.identity())
			tmp2 = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.plus(), pcb.identity())
			return ParVec.toParVec(tmp1+tmp2)
		elif dir == gr.In:
			ret = self.spm.Reduce(pcb.pySpParMat.Row(),pcb.plus(), pcb.identity())
			return ParVec.toParVec(ret)
		elif dir == gr.Out:
			ret = self.spm.Reduce(pcb.pySpParMat.Column(),pcb.plus(), pcb.identity())
			return ParVec.toParVec(ret)
		else:
			raise KeyError, 'Invalid edge direction'

	T = reverseEdges

	def toParVec(self):
		ne = self.nedge()
		reti = ParVec(ne)
		retj = ParVec(ne)
		retv = ParVec(ne)
		self.spm.Find(reti.dpv, retj.dpv, retv.dpv)
		#ToDo:  return nvert() of original graph, too
		return (reti, retj, retv)

	@staticmethod
	def twoDTorus(n):
		N = n*n
		nvec =   ((ParVec.range(N*4)%N) / n).floor()	 # [0,0,0,...., n-1,n-1,n-1]
		nvecil = ((ParVec.range(N*4)%N) % n).floor()	 # [0,1,...,n-1,0,1,...,n-2,n-1]
		north = gr.Graph._sub2ind((n,n),(nvecil-1) % n,nvec)	
		south = gr.Graph._sub2ind((n,n),(nvecil+1) % n,nvec)
		west = gr.Graph._sub2ind((n,n),nvecil, (nvec-1) % n)
		east = gr.Graph._sub2ind((n,n),nvecil, (nvec+1) % n)
		Ndx = ParVec.range(N*4)
		northNdx = Ndx < N
		southNdx = (Ndx >= N) & (Ndx < 2*N)
		westNdx = (Ndx >= 2*N) & (Ndx < 3*N)
		eastNdx = Ndx >= 3*N
		col = ParVec.zeros(N*4)
		col[northNdx] = north
		col[southNdx] = south
		col[westNdx] = west
		col[eastNdx] = east
		row = ParVec.range(N*4) % N
		ret = DiGraph(row, col, 1, N)
		return ret

	UFget = gr.UFget
	UFdownload = gr.UFdownload

	# ==================================================================
	#  "complex ops" implemented below here
	# ==================================================================


	#	creates a breadth-first search tree of a Graph from a starting
	#	set of vertices.  Returns a 1D array with the parent vertex of 
	#	each vertex in the tree; unreached vertices have parent == -1.
	#	sym arg denotes whether graph is symmetric; if not, need to transpose
	#
	def bfsTree(self, root, sym=False):
		"""
		calculates a breadth-first search tree from the edges in the
		passed DiGraph, starting from the root vertex.  "Breadth-first"
		in the sense that all vertices reachable in step i are added
		to the tree before any of the newly-reachable vertices' reachable
		vertices are explored.

		Input Arguments:
			root:  an integer denoting the root vertex for the tree
			sym:  a Boolean denoting whether the DiGraph is symmetric
			    (i.e., each edge from vertex i to vertex j has a
			    companion edge from j to i).  If the DiGraph is 
			    symmetric, the operation is faster.  The default is 
			    False.

		Input Arguments:
			parents:  a ParVec instance of length equal to the number
			    of vertices in the DiGraph, with each element denoting 
			    the vertex number of that vertex's parent in the tree.
			    The root is its own parent.  Unreachable vertices
			    have a parent of -1. 

		SEE ALSO: isBfsTree 
		"""
		if not sym:
			self.T()
		parents = pcb.pyDenseParVec(self.nvert(), -1)
		# NOTE:  values in fringe go from 1:n instead of 0:(n-1) so can
		# distinguish vertex0 from empty element
		fringe = pcb.pySpParVec(self.nvert())
		parents[root] = root
		fringe[root] = root
		while fringe.getnnz() > 0:
			#FIX:  setNumToInd -> SPV.range()
			fringe.setNumToInd()
			self.spm.SpMV_SelMax_inplace(fringe)
			pcb.EWiseMult_inplacefirst(fringe, parents, True, -1)
			parents[fringe] = 0
			parents += fringe
		if not sym:
			self.T()
		return ParVec.toParVec(parents)
	
	
		# returns tuples with elements
		# 0:  True/False of whether it is a BFS tree or not
		# 1:  levels of each vertex in the tree (root is 0, -1 if not reached)
	def isBfsTree(self, root, parents, sym=False):
		"""
		validates that a breadth-first search tree in the style created
		by bfsTree is correct.

		Input Arguments:
			root:  an integer denoting the root vertex for the tree
			parents:  a ParVec instance of length equal to the number
			    vertices in the DiGraph, with each element denoting 
			    the vertex number of that vertex's parent in the tree.
			    The root is its own parent.  Unreachable vertices
			    have a parent of -1. 
			sym:  a Boolean denoting whether the DiGraph is symmetric
			    (i.e., each edge from vertex i to vertex j has a
			    companion edge from j to i).  If the DiGraph is 
			    symmetric, the operation is faster.  The default is 
			    False.
		
		Output Arguments:
			ret:  The return value may be an integer (in the case of
			    an error detected) or a tuple (in the case of no
			    error detected).  If it's an integer, its value will
			    be the negative of the first test below that failed.
			    If it's a tuple, its first element will be the 
			    the integer 1 and its second element will be a 
			    ParVec of length equal to the number of vertices
			    in the DiGraph, with each element denoting the
			    level in the tree at which the vertex resides.  The
			    root resides in level 0, its direct neighbors in
			    level 1, and so forth.  Unreachable vertices have a
			    level value of -1.
		
		Tests:
			The tests implement some of the Graph500 (www.graph500.org) 
			specification. (Some of the Graph500 tests also depend on 
			input edges.)
			1:  The tree does not contain cycles,  that every vertex
			    with a parent is in the tree, and that the root is
			    not the destination of any tree edge.
			2:  Tree edges are between vertices whose levels differ 
			    by exactly 1.
		SEE ALSO: bfsTree 
		"""
		ret = 1		# assume valid
		nvertG = self.nvert()
	
		# calculate level in the tree for each vertex; root is at level 0
		# about the same calculation as bfsTree, but tracks levels too
		if not sym:
			self.reverseEdges()
		parents2 = ParVec.zeros(nvertG) - 1
		parents2[root] = root
		fringe = SpParVec(nvertG)
		fringe[root] = root
		levels = ParVec.zeros(nvertG) - 1
		levels[root] = 0
	
		level = 1
		while fringe.nnn() > 0:
			fringe.sprange()
			#FIX:  create PCB graph-level op
			self.spm.SpMV_SelMax_inplace(fringe.spv)
			#FIX:  create PCB graph-level op
			pcb.EWiseMult_inplacefirst(fringe.spv, parents2.dpv, True, -1)
			parents2[fringe] = fringe
			levels[fringe] = level
			level += 1
		if not sym:
			self.reverseEdges()

		# spec test #1
		# Confirm that the tree is a tree;  i.e., that it does not
		# have any cycles (visited more than once while building
		# the tree) and that every vertex with a parent is
		# in the tree. 

		# build a new graph from just tree edges
		tmp2 = parents != ParVec.range(nvertG)
		treeEdges = (parents != -1) & tmp2
		treeI = parents[treeEdges.findInds()]
		treeJ = ParVec.range(nvertG)[treeEdges.findInds()]
		# root cannot be destination of any tree edge
		if (treeJ == root).any():
			ret = -1
			return ret
		# note treeJ/TreeI reversed, so builtGT is transpose, as
		#   needed by SpMV
		builtGT = DiGraph(treeJ, treeI, 1, nvertG)
		visited = ParVec.zeros(nvertG)
		visited[root] = 1
		fringe = SpParVec(nvertG)
		fringe[root] = root
		cycle = False
		multiparents = False
		while fringe.nnn() > 0 and not cycle and not multiparents:
			fringe.spones()
			newfringe = SpParVec.toSpParVec(builtGT.spm.SpMV_PlusTimes(fringe.spv))
			if visited[newfringe.toParVec().findInds()].any():
				cycle = True
				break
			if (newfringe > 1).any():
				multiparents = True
			fringe = newfringe
			visited[fringe] = 1
		if cycle or multiparents:
			ret = -1
			return ret
		
		# spec test #2
		#    tree edges should be between verts whose levels differ by 1
		
		if (levels[treeI]-levels[treeJ] != -1).any():
			ret = -2
			return ret
	
		return (ret, levels)
	
	# returns a Boolean vector of which vertices are neighbors
	def neighbors(self, source, nhop=1, sym=False):
		if not sym:
			self.T()
		dest = ParVec(self.nvert(),0)
		fringe = SpParVec(self.nvert())
		fringe[source] = 1
		for i in range(nhop):
			self.spm.SpMV_SelMax_inplace(fringe.spv)
			dest[fringe.toParVec()] = 1
		if not sym:
			self.T()
		return dest
		
	# returns:
	#   - source:  a vector of the source vertex for each new vertex
	#   - dest:  a Boolean vector of the new vertices
	#ToDo:  nhop argument?
	def pathsHop(self, source, sym=False):
		if not sym:
			self.T()
		retDest = ParVec(self.nvert(),0)
		retSource = ParVec(self.nvert(),0)
		fringe = SpParVec(self.nvert())
		retDest[fringe] = 1
		fringe.sprange()
		self.spm.SpMV_SelMax_inplace(fringe.spv)
		retDest[fringe] = 1
		retSource[fringe] = fringe
		if not sym:
			self.T()
		return (retSource, retDest)
		
	def centrality(self, alg, **kwargs):
		"""
		where 'alg' can be one of 
		    'exactBC':  exact betweenness centrality
		    'approxBC':  approximate betweenness centrality

		Each algorithm may have algorithm-specific arguments as follows:
		    'exactBC':  
		        normalize=True:  normalizes the values by dividing by 
		                (nVert-1)*(nVert-2)
		    'approxBC':
			sample=0.05:  the fraction of the vertices to use as sources 
				and destinations;  sample=1.0 is the same as exactBC
		        normalize=True:  normalizes the values by multiplying by 
				nVerts / (nVertsCalculated * (nVerts-1) * (nVerts-2))
		The return value is a ParVec with length equal to the number of
		vertices in the DiGraph, with each element of the ParVec containing
		the centrality value of the vertex.
		"""
		if alg=='exactBC':
			cent = DiGraph._approxBC(self, sample=1.0, **kwargs)
			#cent = DiGraph._bc(self, 1.0, self.nvert())
		elif alg=='approxBC':
			cent = DiGraph._approxBC(self, **kwargs)
		elif alg=='kBC':
			raise NotImplementedError, "k-betweenness centrality unimplemented"
		elif alg=='degree':
			raise NotImplementedError, "degree centrality unimplemented"
		else:
			raise KeyError, "unknown centrality algorithm (%s)" % alg
	
		return cent
	
	
	def _approxBC(self, sample=0.05, normalize=True, nProcs=pcb._nprocs(), BCdebug=0):
		A = self.copy()
		self.ones()
		#Aint = self.ones()	# not needed;  Gs only int for now
		N = A.nvert()
		bc = ParVec(N)
		#nProcs = pcb._nprocs()
		nVertToCalc = int(self.nvert() * sample)
		# batchSize = #rows/cols that will fit in memory simultaneously.
		# bcu has a value in every element, even though it's literally
		# a sparse matrix (DiGraph).  So batchsize is calculated as
		#   nrow = memory size / (memory/row)
		#   memory size (in edges)
		#        = 2GB * 0.1 (other vars) / 18 (bytes/edge) * nProcs
		#   memory/row (in edges)
		#        = self.nvert()
		physMemPCore = 2e9; memFract = 0.1; bytesPEdge = 18
		batchSize = int(2e9 * memFract / bytesPEdge * nProcs / N)
		nBatches = int(sc.ceil(float(nVertToCalc) / float(batchSize)))
		nPossBatches = int(sc.ceil(float(N) / float(batchSize)))
		if sample == 1.0:
			startVs = range(0,nVertToCalc,batchSize)
			endVs = range(batchSize, nVertToCalc, batchSize)
			if nVertToCalc % batchSize != 0:
				endVs.append(nVertToCalc)
			numVs = [y-x for [x,y] in zip(startVs,endVs)]
		else:
			startVs = (sc.random.randint(0,nPossBatches,nBatches)*batchSize).tolist()
			numVs = [min(x+batchSize,N)-x for x in startVs]

		if BCdebug:
			print "batchSz=%d, nBatches=%d, nPossBatches=%d" % (batchSize, nBatches, nPossBatches)
		for [startV, numV] in zip(startVs, numVs):
			if BCdebug:
				print "startV=%d, numV=%d" % (startV, numV)
			bfs = []		
			batch = ParVec.range(startV, startV+numV)
			curSize = len(batch)
			nsp = DiGraph(ParVec.range(curSize), batch, 1, curSize, N)
			fringe = A[batch,ParVec.range(N)]
			depth = 0
			while fringe.nedge() > 0:
				depth = depth+1
				nsp = nsp+fringe
				tmp = fringe.copy()
				tmp.bool()
				bfs.append(tmp)
				tmp = fringe._SpMM(A)
				fringe = tmp.mulNot(nsp)
	
			bcu = DiGraph.fullyConnected(curSize,N)
			# compute the bc update for all vertices except the sources
			for depth in range(depth-1,0,-1):
				# compute the weights to be applied based on the child values
				w = bfs[depth] / nsp * bcu
				# Apply the child value weights and sum them up over the parents
				# then apply the weights based on parent values
				w.T()
				w = A._SpMM(w)
				w.T()
				w *= bfs[depth-1]
				w *= nsp
				bcu += w
	
			# update the bc with the bc update
			bc = bc + bcu.sum(Out)	# column sums
	
		# subtract off the additional values added in by precomputation
		bc = bc - nVertToCalc
		if normalize:
			nVertSampled = sum(numVs)
			bc = bc * (float(N)/float(nVertSampled*(N-1)*(N-2)))
		return bc
	
	def cluster(self, alg, **kwargs):
	#		ToDo:  Normalize option?
		"""
		Deferred implementation for KDT v0.1
			
		"""
		raise NotImplementedError, "clustering not implemented for v0.1"
		if alg=='Markov' or alg=='markov':
			clus = _markov(self, **kwargs)
	
		elif alg=='kNN' or alg=='knn':
			raise NotImplementedError, "k-nearest neighbors clustering not implemented"
	
		else:
			raise KeyError, "unknown clustering algorithm (%s)" % alg
	
		return clus

class ParVec(gr.ParVec):
	pass

class SpParVec(gr.SpParVec):
	pass

		

master = gr.master
sendFeedback = gr.sendFeedback
InOut = gr.InOut
In = gr.In
Out = gr.Out
