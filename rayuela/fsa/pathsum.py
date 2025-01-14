from collections import defaultdict
import numpy as np
from numpy import linalg as LA
from frozendict import frozendict

from rayuela.base.datastructures import PriorityQueue
from rayuela.base.semiring import Real, Semiring
from rayuela.fsa.push import push_with_potential

from rayuela.fsa.state import State
from rayuela.fsa.scc import SCC

class Strategy:
	VITERBI = 1
	BELLMANFORD = 2
	DIJKSTRA = 3
	LEHMANN = 4
	JOHNSON = 5
	FIXPOINT = 6
	DECOMPOSED_LEHMANN = 7

class Pathsum:

	def __init__(self, fsa):

		# basic FSA stuff
		self.fsa = fsa
		self.R = fsa.R
		self.N = self.fsa.num_states

		# state dictionary
		self.I = {}
		for n, q in enumerate(self.fsa.Q):
			self.I[q] = n

		# lift into the semiring
		self.W = self.lift()

	def _convert(self):
		mat = np.zeros((self.N, self.N))
		for n in range(self.N):
			for m in range(self.N):
				mat[n, m] = self.W[n, m].score
		return mat

	def max_eval(self):
		# computes the largest eigenvalue
		mat = self._convert()
		if len(mat) == 0:
			return 0.0
		vals = []
		for val in LA.eigvals(mat):
			vals.append(np.abs(val))
		return np.max(vals)

	def lift(self):
		""" creates the weight matrix from the automaton """
		W = self.R.zeros(self.N, self.N)
		for p in self.fsa.Q:
			for a, q, w in self.fsa.arcs(p):
				W[self.I[p], self.I[q]] += w
		return W

	def pathsum(self, strategy):
		if strategy == Strategy.DIJKSTRA:
			assert self.R.superior, "Dijkstra's requires a superior semiring"
			return self.dijkstra_early()

		elif strategy == Strategy.VITERBI:
			assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
			return self.viterbi_pathsum()

		elif strategy == Strategy.BELLMANFORD:
			assert self.R.idempotent, "Bellman-Ford requires an idempotent semiring"
			return self.bellmanford_pathsum()

		elif strategy == Strategy.JOHNSON:
			assert self.R.idempotent, "Johnson's requires an idempotent semiring"
			return self.johnson_pathsum()

		elif strategy == Strategy.LEHMANN:
			return self.lehmann_pathsum()

		elif strategy == Strategy.FIXPOINT:
			return self.fixpoint_pathsum()

		elif strategy == Strategy.DECOMPOSED_LEHMANN:
			return self.decomposed_lehmann_pathsum()

		else:
			raise NotImplementedError

	def forward(self, strategy):

		if strategy == Strategy.DIJKSTRA:
			assert self.R.superior, "Dijkstra's requires a superior semiring"
			return self.dijkstra_fwd()

		if strategy == Strategy.VITERBI:
			assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
			return self.viterbi_fwd()

		elif strategy == Strategy.BELLMANFORD:
			assert self.R.idempotent, "Bellman-Ford requires an idempotent semiring"
			return self.bellmanford_fwd()

		elif strategy == Strategy.JOHNSON:
			assert self.R.idempotent, "Johnson's requires an idempotent semiring"
			return self.johnson_fwd()

		elif strategy == Strategy.LEHMANN:
			return self.lehmann_fwd()

		elif strategy == Strategy.FIXPOINT:
			return self.fixpoint_fwd()

		else:
			raise NotImplementedError

	def backward(self, strategy):
		if strategy == Strategy.VITERBI:
			assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
			return self.viterbi_bwd()

		elif strategy == Strategy.BELLMANFORD:
			assert self.R.idempotent, "Bellman-Ford requires an idempotent semiring"
			return self.bellmanford_bwd()

		elif strategy == Strategy.JOHNSON:
			assert self.R.idempotent, "Johnson's requires an idempotent semiring"
			return self.johnson_bwd()

		elif strategy == Strategy.LEHMANN:
			return self.lehmann_bwd()

		elif strategy == Strategy.FIXPOINT:
			return self.fixpoint_bwd()

		else:
			raise NotImplementedError

	def allpairs(self, strategy):
		if strategy == Strategy.JOHNSON:
			assert self.R.idempotent, "Johnson's requires an idempotent semiring"
			return self.johnson()

		elif strategy == Strategy.LEHMANN:
			return self.lehmann()

		elif strategy == Strategy.FIXPOINT:
			raise self.fixpoint()

		else:
			raise NotImplementedError

	def allpairs_pathsum(self, W):
		pathsum = self.R.zero
		for p in self.fsa.Q:
			for q in self.fsa.Q:
				pathsum += self.fsa.λ[p] * W[p, q] * self.fsa.ρ[q]
		return pathsum

	def allpairs_fwd(self, W):
		α = self.R.chart()
		for p in self.fsa.Q:
			for q in self.fsa.Q:
				α[q] += self.fsa.λ[p] * W[p, q]
		return frozendict(α)

	def allpairs_bwd(self, W):
		pass
		𝜷 = self.R.chart()
		W = self.lehmann()
		for p in self.fsa.Q:
			for q in self.fsa.Q:
				𝜷[p] += W[p, q] * self.fsa.ρ[q]
		return frozendict(𝜷)

	def viterbi_pathsum(self, forward=False):
		pathsum = self.R.zero
		if forward:
			alpha = self.viterbi_fwd()
			for q in self.fsa.Q:
				pathsum += alpha[q] * self.fsa.ρ[q]
		else:
			𝜷 = self.viterbi_bwd()
			for q in self.fsa.Q:
				pathsum += self.fsa.λ[q] * 𝜷[q]
		return pathsum

	def viterbi_fwd(self) -> "defaultdict[State, Semiring]":
		# Homework 2: Question 2
		assert self.fsa.acyclic
		alpha = self.R.chart()
		# base
		for q, w in self.fsa.I:
			alpha[q] = w

		# recursion
		for p in self.fsa.toposort():
			for _, q, w, in self.fsa.arcs(p):
				alpha[q] += alpha[p] * w

		return alpha

	def viterbi_bwd(self) -> "defaultdict[State, Semiring]":
		""" The Viterbi algorithm run backwards"""

		assert self.fsa.acyclic

		# chart
		𝜷 = self.R.chart()

		# base case (paths of length 0)
		for q, w in self.fsa.F:
			𝜷[q] = w

		# recursion
		for p in self.fsa.toposort(rev=True):
			for _, q, w in self.fsa.arcs(p):
				𝜷[p] += w * 𝜷[q]

		return 𝜷

	def dijkstra_early(self):
		""" Dijkstra's algorithm with early stopping."""
		raise NotImplementedError


	def dijkstra_fwd(self, I=None):
		""" Dijkstra's algorithm without early stopping. """

		assert self.fsa.R.superior

		# initialization
		α = self.R.chart()
		agenda = PriorityQueue(R=self.fsa.R)
		popped = set([])

		# base case
		if I is None:
			for q, w in self.fsa.I:
				agenda.push(q, w)
		else:
			for q in I:
				agenda.push(q, self.R.one)

		# main loop
		while agenda:
			i, v = agenda.pop()
			popped.add(i)
			α[i] += v

			for _, j, w in self.fsa.arcs(i):
				if j not in popped:
					agenda.push(j, v * w)

		return α

	def _lehmann(self, zero=True):
		"""
		Lehmann's (1977) algorithm.
		"""

		# initialization
		V = self.W.copy()
		U = self.W.copy()


		# basic iteration
		for j in range(self.N):
			V, U = U, V
			V = self.R.zeros(self.N, self.N)
			for i in range(self.N):
				for k in range(self.N):
					# i ➙ j ⇝ j ➙ k
					V[i,k] = U[i,k] + U[i,j] * U[j,j].star() * U[j,k]

		# post-processing (paths of length zero)
		if zero:
			for i in range(self.N):
				V[i,i] += self.R.one


		return V

	def lehmann(self, zero=True):

		V = self._lehmann(zero=zero)

		W = {}
		for p in self.fsa.Q:
			for q in self.fsa.Q:
				if p in self.I and q in self.I:
					W[p, q] = V[self.I[p], self.I[q]]
				elif p == q and zero:
					W[p, q] = self.R.one
				else:
					W[p, q] = self.R.zero

		return frozendict(W)

	def lehmann_pathsum(self): return self.allpairs_pathsum(self.lehmann())
	def lehmann_fwd(self): return self.allpairs_fwd(self.lehmann())
	def lehmann_bwd(self): return self.allpairs_bwd(self.lehmann())

	def local_lehmann(self, component: set):
		N = len(component)
		I = {}  # state to index
		for i, q in enumerate(component):
			I[q] = i

		# Initialization
		W = self.R.zeros(N, N)
		for p in self.fsa.Q:
			for a, q, w in self.fsa.arcs(p):
				if p in component and q in component:
					W[I[p], I[q]] += w

		for j in range(N):
			V = self.R.zeros(N, N)
			for i in range(N):
				for k in range(N):
					# i ➙ j ⇝ j ➙ k
					V[i, k] = W[i, k] + W[i, j] * W[j, j].star() * W[j, k]
			W = V

		# Paths of zero length
		for i in range(N):
			W[i, i] += self.R.one

		return {(p, q): W[I[p], I[q]] for p in component for q in component}

	def decomposed_lehmann_bwd(self):
		beta = self.R.chart()
		# base
		for q, w in self.fsa.F:
			beta[q] = w

		scc_decomp = SCC(self.fsa)
		for scc in reversed(scc_decomp.scc()):
			# Run inter-component Viterbi backward algorithm
			for p in scc:
				for a, q, w in self.fsa.arcs(p):
					if q not in scc:
						beta[p] += w * beta[q]

			# Run in-component Lehmann's algorithm
			W = self.local_lehmann(scc)

			# Run in-component Viterbi backward algorithm
			gamma = self.R.chart()
			for p in scc:
				for q in scc:
					gamma[p] += W[p, q] * beta[q]
			for q in scc:
				beta[q] = gamma[q]

		return beta

	def decomposed_lehmann_pathsum(self) -> Semiring:
		# Homework 3: Question 4
		beta = self.decomposed_lehmann_bwd()

		pathsum = self.R.zero
		for q in self.fsa.Q:
			pathsum += self.fsa.λ[q] * beta[q]

		return pathsum

	def bellmanford_pathsum(self) -> Semiring:
		pathsum = self.R.zero
		𝜷 = self.bellmanford_bwd()
		for q in self.fsa.Q:
			pathsum += self.fsa.λ[q] * 𝜷[q]
		return pathsum

	def bellmanford_fwd(self, I=None) -> "frozendict[State, Semiring]":
		alpha = self.R.chart()  # distance from source
		# base case
		if I is None:
			for q, w in self.fsa.I:
				alpha[q] = w
		else:
			for q in I:
				alpha[q] = self.R.one
		# Relax edges N-1 times
		for i in range(self.N - 1):
			for p in self.fsa.Q:
				for a, q, w in self.fsa.arcs(p):
					alpha[q] += alpha[p] * w
		# Check negative-cycle
		for p in self.fsa.Q:
			for a, q, w in self.fsa.arcs(p):
				if alpha[q] + alpha[p] * w != alpha[q]:
					raise AttributeError("Graph contains a negative-weight cycle")

		return frozendict(alpha)

	def bellmanford_bwd(self, F=None) -> "frozendict[State, Semiring]":
		beta = self.R.chart()  # distance from source
		# base case
		if F is None:
			for q, w in self.fsa.F:
				beta[q] = w
		else:
			for q in F:
				beta[q] = self.R.one
		# Relax edges N-1 times
		for i in range(self.N - 1):
			for p in self.fsa.Q:
				for a, q, w in self.fsa.arcs(p):
					beta[p] += w * beta[q]
		# Check negative-cycle
		for p in self.fsa.Q:
			for a, q, w in self.fsa.arcs(p):
				if beta[p] + w * beta[q] != beta[p]:
					raise AttributeError("Graph contains a negative-weight cycle")

		return frozendict(beta)


	def johnson(self) -> "defaultdict[(State,State), Semiring]":
		# 1
		alpha = self.bellmanford_fwd(I=[q for q, w in self.fsa.I])
		# 2
		V = {q: ~w for q, w in alpha.items()}
		pfsa = push_with_potential(self.fsa, V, False)
		ps = Pathsum(pfsa)
		# 3
		W = self.R.chart()
		for p in self.fsa.Q:
			d = ps.dijkstra_fwd(I=[p])
			for q, w in d.items():
				W[p, q] = ~alpha[p] * w * alpha[q]
		return W

	def johnson_pathsum(self): return self.allpairs_pathsum(self.johnson())
	def johnson_fwd(self): return self.allpairs_fwd(self.johnson())
	def johnson_bwd(self): return self.allpairs_bwd(self.johnson())

	def fixpoint(self):
		raise NotImplementedError

	def fixpoint_pathsum(self): return self.allpairs_pathsum(self.fixpoint())
	def fixpoint_fwd(self): return self.allpairs_fwd(self.fixpoint())
	def fixpoint_bwd(self): return self.allpairs_bwd(self.fixpoint())
