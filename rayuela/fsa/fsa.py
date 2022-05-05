from __future__ import annotations

import copy
from collections import defaultdict as dd, deque
from itertools import product

from frozendict import frozendict

from rayuela.base.misc import epsilon_filter
from rayuela.base.partitions import PartitionRefinement
from rayuela.base.semiring import Boolean, Semiring, String, ProductSemiring
from rayuela.base.symbol import Sym, ε, ε_1, ε_2
from rayuela.fsa.pathsum import Pathsum, Strategy
from rayuela.fsa.state import State, PairState, PowerState
from rayuela.fsa.transformer import Transformer


class FSA:

	def __init__(self, R=Boolean):

		# DEFINITION
		# A weighted finite-state automaton is a 5-tuple <R, Σ, Q, δ, λ, ρ> where
		# • R is a semiring;
		# • Σ is an alphabet of symbols;
		# • Q is a finite set of states;
		# • δ is a finite relation Q × Σ × Q × R;
		# • λ is an initial weight function;
		# • ρ is a final weight function.

		# NOTATION CONVENTIONS
		# • single states (elements of Q) are denoted q
		# • multiple states not in sequence are denoted, p, q, r, ...
		# • multiple states in sequence are denoted i, j, k, ...
		# • symbols (elements of Σ) are denoted lowercase a, b, c, ...
		# • single weights (elements of R) are denoted w
		# • multiple weights (elements of R) are denoted u, v, w, ...

		# semiring
		self.R = R

		# alphabet of symbols
		self.Sigma = set([])

		# a finite set of states
		self.Q = set([])

		# transition function : Q × Σ × Q → R
		self.δ = dd(lambda : dd(lambda : dd(lambda : self.R.zero)))

		# initial weight function
		self.λ = R.chart()

		# final weight function
		self.ρ = R.chart()

	def add_state(self, q):
		if not isinstance(q, State):
			q = State(q)
		self.Q.add(q)

	def add_states(self, Q):
		for q in Q:
			self.add_state(q)

	def add_arc(self, i, a, j, w=None):
		if w is None: w = self.R.one

		if not isinstance(a, Sym): a = Sym(a)
		if not isinstance(w, self.R): w = self.R(w)

		self.add_states([i, j])
		self.Sigma.add(a)
		self.δ[i][a][j] += w

	def set_arc(self, i, a, j, w):
		self.add_states([i, j])
		self.Sigma.add(a)
		self.δ[i][a][j] = w

	def set_I(self, q, w=None):
		if w is None: w = self.R.one
		self.add_state(q)
		self.λ[q] = w

	def set_F(self, q, w=None):
		if w is None: w = self.R.one
		self.add_state(q)
		self.ρ[q] = w

	def add_I(self, q, w):
		self.add_state(q)
		self.λ[q] += w

	def add_F(self, q, w):
		self.add_state(q)
		self.ρ[q] += w

	def freeze(self):
		self.Sigma = frozenset(self.Sigma)
		self.Q = frozenset(self.Q)
		self.δ = frozendict(self.δ)
		self.λ = frozendict(self.λ)
		self.ρ = frozendict(self.ρ)

	@property
	def I(self):
		for q, w in self.λ.items():
			if w != self.R.zero:
				yield q, w

	@property
	def F(self):
		for q, w in self.ρ.items():
			if w != self.R.zero:
				yield q, w

	def arcs(self, i, no_eps=False):
		for a, T in self.δ[i].items():
			if no_eps and a == ε:
				continue
			for j, w in T.items():
				if w == self.R.zero:
					continue
				yield a, j, w

	def accept(self, string):
		""" determines whether a string is in the language """
		assert isinstance(string, str)

		fsa = FSA(R=self.R)
		for i, x in enumerate(list(string)):
			fsa.add_arc(State(i), Sym(x), State(i+1), self.R.one)

		fsa.set_I(State(0), self.R.one)
		fsa.add_F(State(len(string)), self.R.one)

		return self.intersect(fsa).pathsum()

	@property
	def num_states(self):
		return len(self.Q)

	def copy(self):
		""" deep copies the machine """
		return copy.deepcopy(self)

	def spawn(self, keep_init=False, keep_final=False):
		""" returns a new FSA in the same semiring """
		F = FSA(R=self.R)

		if keep_init:
			for q, w in self.I:
				F.set_I(q, w)
		if keep_final:
			for q, w in self.F:
				F.set_F(q, w)

		return F
	
	def push(self):
		return Transformer.push(self)

	def determinize(self) -> FSA:
		# Homework 4: Question 4
		det = self.spawn()
		QI = PowerState(dict(self.I))
		det.set_I(QI)
		stack = [QI]

		while len(stack) > 0:
			Q = stack.pop()
			for a, QP, wp in Transformer.powerarcs(self, Q):
				# New power state:
				is_new_state = QP not in det.Q
				det.add_arc(Q, a, QP, wp)
				if is_new_state:
					# Final?
					det.ρ[QP] = sum([rq * self.ρ[q] for q, rq in QP.residuals.items()], start=self.R.zero)
					stack.append(QP)

		return det

	def minimize(self, strategy=None) -> FSA:
		# Homework 5: Question 3
		assert self.deterministic

		P = {frozenset(q for q in self.Q if self.ρ[q] != self.R.zero),
			 frozenset(q for q in self.Q if self.ρ[q] == self.R.zero)}
		for a in self.Sigma:
			fa = {}
			for p in self.Q:
				for b, q, w in self.arcs(p):
					if b == a:
						fa[p] = q
						break
				else:
					raise NotImplementedError
			P = PartitionRefinement(fa, self.Q).hopcroft(P)

		return self.block_fsa_construction(P)

	def block_fsa_construction(self, P) -> FSA:
		fsa = self.spawn()
		mu = {}
		for Q in P:
			for q in Q:
				mu[q] = Q

		for Q in P:
			for p in Q:
				for a, q, w in self.arcs(p):
					fsa.add_arc(mu[p], a, mu[q], w)

		for q, w in self.I:
			fsa.λ[mu[q]] += w

		for q, w in self.F:
			fsa.ρ[mu[q]] += w

		return fsa

	def dfs(self):
		""" Depth-first search (Cormen et al. 2019; Section 22.3) """

		in_progress, finished = set([]), {}
		cyclic, counter = False, 0

		def _dfs(p):
			nonlocal in_progress
			nonlocal finished
			nonlocal cyclic
			nonlocal counter

			in_progress.add(p)

			for _, q, _ in self.arcs(p):
				if q in in_progress:
					cyclic = True
				elif q not in finished:
					_dfs(q)

			in_progress.remove(p)
			finished[p] = counter
			counter += 1

		for q, _ in self.I: _dfs(q)

		return cyclic, finished

	def finish(self, rev=False, acyclic_check=False):
		"""
		Returns the nodes in order of their finishing time.
		"""

		cyclic, finished = self.dfs()

		if acyclic_check:
			assert self.acyclic

		sort = {}
		for s, n in finished.items():
			sort[n] = s
		if rev:
			for n in sorted(list(sort.keys())):
				yield sort[n]
		else:
			for n in reversed(sorted(list(sort.keys()))):
				yield sort[n]

	def toposort(self, rev=False):
		return self.finish(rev=rev, acyclic_check=True)

	@property
	def acyclic(self) -> bool:
		cyclic, _ = self.dfs()
		return not cyclic

	@property
	def deterministic(self) -> bool:
		# Homework 1: Question 2
		for i in self.Q:
			for a in self.δ[i]:
				T = self.δ[i][a]
				n_nonzero = sum(1 for w in T.values() if w != self.R.zero)
				if n_nonzero > 1:
					return False
		return True

	@property
	def pushed(self) -> bool:
		# Homework 1: Question 2
		for i in self.Q:
			w_sum = self.R.zero
			for a in self.δ[i]:
				for w in self.δ[i][a].values():
					w_sum += w
			if w_sum + self.ρ[i] != self.R.one:
				return False
		return True

	def reverse(self) -> FSA:
		""" computes the reverse of the FSA """
		# Homework 1: Question 3
		rev = FSA(self.R)
		for i in self.Q:
			for a, j, w in self.arcs(i):
				rev.add_arc(j, a, i, w)
		for q, w in self.I:
			rev.add_F(q, w)
		for q, w in self.F:
			rev.add_I(q, w)
		return rev

	def accessible(self) -> set:
		""" computes the set of acessible states """
		# Homework 1: Question 3
		visited = set()
		# BFS from initial states
		queue = deque(q for q, w in self.I)
		while queue:
			i = queue.popleft()
			if i in visited: continue
			visited.add(i)
			for a, j, w in self.arcs(i):
				if j not in visited:
					queue.append(j)
		return visited

	def coaccessible(self) -> set:
		""" computes the set of co-acessible states """
		# Homework 1: Question 3
		return self.reverse().accessible()

	def trim(self) -> FSA:
		""" keeps only those states that are both accessible and co-accessible """
		# Homework 1: Question 3
		remaining_states = self.accessible().intersection(self.coaccessible())
		ret = FSA(self.R)
		for i in self.Q:
			if i in remaining_states:
				for a, j, w in self.arcs(i):
					if j in remaining_states:
						ret.add_arc(i, a, j, w)
		for q, w in self.I:
			if q in remaining_states:
				ret.set_I(q, w)
		for q, w in self.F:
			if q in remaining_states:
				ret.set_F(q, w)
		return ret

	def union(self, fsa) -> FSA:
		""" construct the union of the two FSAs """
		# Homework 1: Question 4
		assert self.R == fsa.R
		ret = FSA(self.R)
		for i in self.Q:
			for a, j, w in self.arcs(i):
				ret.add_arc(PairState(1, i.idx), a, PairState(1, j.idx), w)
		for i in fsa.Q:
			for a, j, w in fsa.arcs(i):
				ret.add_arc(PairState(2, i.idx), a, PairState(2, j.idx), w)
		for q, w in self.I:
			ret.set_I(PairState(1, q.idx), w)
		for q, w in self.F:
			ret.set_F(PairState(1, q.idx), w)
		for q, w in fsa.I:
			ret.set_I(PairState(2, q.idx), w)
		for q, w in fsa.F:
			ret.set_F(PairState(2, q.idx), w)
		return ret

	def concatenate(self, fsa) -> FSA:
		""" construct the concatenation of the two FSAs """
		# Homework 1: Question 4
		assert self.R == fsa.R
		ret = FSA(self.R)
		for i in self.Q:
			for a, j, w in self.arcs(i):
				ret.add_arc(PairState(1, i.idx), a, PairState(1, j.idx), w)
		for i in fsa.Q:
			for a, j, w in fsa.arcs(i):
				ret.add_arc(PairState(2, i.idx), a, PairState(2, j.idx), w)
		mid_state = State('*')
		ret.add_state(mid_state)
		for q, w in self.I:
			ret.set_I(PairState(1, q.idx), w)
		for q, w in self.F:
			ret.add_arc(PairState(1, q.idx), ε, mid_state, w)
		for q, w in fsa.I:
			ret.add_arc(mid_state, ε, PairState(2, q.idx), w)
		for q, w in fsa.F:
			ret.set_F(PairState(2, q.idx), w)
		return ret

	def kleene_closure(self) -> FSA:
		""" compute the Kleene closure of the FSA """
		# Homework 1: Question 4
		ret = FSA(self.R)
		initial_state = State('*I')
		final_state = State('*F')
		ret.add_state(initial_state)
		ret.add_state(final_state)
		ret.set_I(initial_state)
		ret.set_F(final_state)
		for i in self.Q:
			for a, j, w in self.arcs(i):
				ret.add_arc(i, a, j, w)
		for q, w in self.I:
			ret.add_arc(initial_state, ε, q, w)
		for q, wf in self.F:
			ret.add_arc(q, ε, final_state, wf)
			for p, wi in self.I:
				ret.add_arc(q, ε, p, wf * wi)
		ret.add_arc(initial_state, ε, final_state)
		return ret

	def pathsum(self, strategy=Strategy.LEHMANN):
		if self.acyclic:
			strategy = Strategy.VITERBI
		pathsum = Pathsum(self)
		return pathsum.pathsum(strategy)

	def forward(self, strategy=Strategy.LEHMANN):
		if self.acyclic:
			strategy = Strategy.VITERBI
		return Pathsum(self).forward(strategy)

	def backward(self, strategy=Strategy.LEHMANN):
		if self.acyclic:
			strategy = Strategy.VITERBI
		return Pathsum(self).backward(strategy)

	def edge_marginals(self) -> dd[dd[dd[Semiring]]]:
		""" computes the edge marginals μ(q→q') """
		# Homework 2: Question 2
		assert self.acyclic
		ps = Pathsum(self)
		alpha, beta = ps.viterbi_fwd(), ps.viterbi_bwd()

		res = dd(lambda: dd(lambda: dd(lambda: self.R.zero)))
		for p in self.Q:
			for a, q, w in self.arcs(p):
				res[p][a][q] = alpha[p] * w * beta[q]

		return res

	def intersect(self, fsa):
		"""
		on-the-fly weighted intersection
		"""

		# the two machines need to be in the same semiring
		assert self.R == fsa.R

		# add initial states
		product_fsa = FSA(R=self.R)
		for (q1, w1), (q2, w2) in product(self.I, fsa.I):
			product_fsa.add_I(PairState(q1, q2), w=w1 * w2)

		self_initials = {q: w for q, w in self.I}
		fsa_initials = {q: w for q, w in fsa.I}

		visited = set([(i1, i2, State('0')) for i1, i2 in product(self_initials, fsa_initials)])
		stack = [(i1, i2, State('0')) for i1, i2 in product(self_initials, fsa_initials)]

		self_finals = {q: w for q, w in self.F}
		fsa_finals = {q: w for q, w in fsa.F}

		while stack:
			q1, q2, qf = stack.pop()

			E1 = [(a if a != ε else ε_2, j, w) for (a, j, w) in self.arcs(q1)] + \
                            [(ε_1, q1, self.R.one)]
			E2 = [(a if a != ε else ε_1, j, w) for (a, j, w) in fsa.arcs(q2)] + \
                            [(ε_2, q2, self.R.one)]

			M = [((a1, j1, w1), (a2, j2, w2))
				 for (a1, j1, w1), (a2, j2, w2) in product(E1, E2)
				 if epsilon_filter(a1, a2, qf) != State('⊥')]

			for (a1, j1, w1), (a2, j2, w2) in M:

				product_fsa.set_arc(
					PairState(q1, q2), a1,
					PairState(j1, j2), w=w1*w2)

				_qf = epsilon_filter(a1, a2, qf)
				if (j1, j2, _qf) not in visited:
					stack.append((j1, j2, _qf))
					visited.add((j1, j2, _qf))

			# final state handling
			if q1 in self_finals and q2 in fsa_finals:
				product_fsa.add_F(
					PairState(q1, q2), w=self_finals[q1] * fsa_finals[q2])

		return product_fsa

	def equivalent(self, fsa) -> bool:
		""" Tests equivalence. """

		# Homework 5: Question 4
		raise NotImplementedError

	def coaccessible_intersection(self, fsa) -> FSA:
		# Homework 2: Question 3
		return self.reverse().intersect(fsa.reverse()).reverse()

	def tikz(self, max_per_row=4):

		tikz_string = []
		previous_ids, positioning = [], ''
		rows = {}

		initial = {q: w for q, w in self.I}
		final = {q: w for q, w in self.F}

		for jj, q in enumerate(self.Q):
			options = 'state'
			additional = ''

			if q in initial:
				options += ', initial'
				additional = f' / {initial[q]}'
			if q in final:
				options += ', accepting'
				additional = f' / {final[q]}'

			if jj >= max_per_row:
				positioning = f'below = of {previous_ids[jj - max_per_row]}'
			elif len(previous_ids) > 0:
				positioning = f'right = of {previous_ids[-1]}'
			previous_ids.append(f'q{q.idx}')
			rows[q] = jj // max_per_row

			tikz_string.append(f'\\node[{options}] (q{q.idx}) [{positioning}] {{ ${q.idx}{additional}$ }}; \n')

		tikz_string.append('\\draw')

		seen_pairs, drawn_pairs = set(), set()

		for jj, q in enumerate(self.Q):
			target_edge_labels = dict()
			for a, j, w in self.arcs(q):
				if j not in target_edge_labels:
					target_edge_labels[j] = f'{a}/{w}'
				else:
					target_edge_labels[j] += f'\\\\{a}/{w}'
				seen_pairs.add(frozenset([q, j]))

			for ii, (target, label) in enumerate(target_edge_labels.items()):

				edge_options = 'align=left'
				if q == target:
					edge_options += ', loop above'
				elif frozenset([q, target]) not in seen_pairs:
					edge_options += 'a, bove'
				elif frozenset([q, target]) not in drawn_pairs:
					if rows[q] == rows[target]:
						edge_options += ', bend left, above'
					else:
						edge_options += ', bend left, right'
				else:
					if rows[q] == rows[target]:
						edge_options += ', bend left, below'
					else:
						edge_options += ', bend left, right'
				end = '\n'
				if jj == self.num_states - 1 and ii == len(target_edge_labels) - 1:
					end = '; \n'
				tikz_string.append(f'(q{q.idx}) edge[{edge_options}] node{{ ${label}$ }} (q{target.idx}) {end}')
				drawn_pairs.add(frozenset([q, j]))

		if not len(list(self.arcs(list(self.Q)[-1]))) > 0:
			tikz_string.append(';')

		return ''.join(tikz_string)

	def __truediv__(self, other):
		return self.intersect(other)

	def __add__(self, other):
		return self.concatenate(other)

	def __sub__(self, other):
		return self.difference(other)

	def __repr__(self):
		return f'WFSA({self.num_states} states, {self.R})'

	def __str__(self):
		""" ascii visualize """

		output = []
		for q, w in self.I:
			output.append(f"initial state:\t{q.idx}\t{w}")
		for q, w in self.F:
			output.append(f"final state:\t{q.idx}\t{w}")
		for p in self.Q:
			for a, q, w in self.arcs(p):
				output.append(f"{p}\t----{a}/{w}---->\t{q}")
		return "\n".join(output)

	def _repr_html_(self):
		"""
		When returned from a Jupyter cell, this will generate the FST visualization
		Based on: https://github.com/matthewfl/openfst-wrapper
		"""
		from uuid import uuid4
		import json
		from collections import defaultdict
		ret = []
		if self.num_states == 0:
			return '<code>Empty FST</code>'

		if self.num_states > 64:
			return f'FST too large to draw graphic, use fst.ascii_visualize()<br /><code>FST(num_states={self.num_states})</code>'

		finals = {q for q, _ in self.F}
		initials = {q for q, _ in self.I}

		# print initial
		for q, w in self.I:
			if q in finals:
				label = f'{str(q)} / [{str(w)} / {str(self.ρ[q])}]'
				color = 'af8dc3'
			else:
				label = f'{str(q)} / {str(w)}'
				color = '66c2a5'

			ret.append(
				f'g.setNode("{repr(q)}", {{ label: {json.dumps(label)} , shape: "circle" }});\n')
				# f'g.setNode("{repr(q)}", {{ label: {json.dumps(hash(label) // 1e8)} , shape: "circle" }});\n')

			ret.append(f'g.node("{repr(q)}").style = "fill: #{color}"; \n')

		# print normal
		for q in (self.Q - finals) - initials:

			label = str(q)

			ret.append(
				f'g.setNode("{repr(q)}", {{ label: {json.dumps(label)} , shape: "circle" }});\n')
				# f'g.setNode("{repr(q)}", {{ label: {json.dumps(hash(label) // 1e8)} , shape: "circle" }});\n')
			ret.append(f'g.node("{repr(q)}").style = "fill: #8da0cb"; \n')

		# print final
		for q, w in self.F:
			# already added
			if q in initials:
				continue

			if w == self.R.zero:
				continue
			label = f'{str(q)} / {str(w)}'

			ret.append(
				f'g.setNode("{repr(q)}", {{ label: {json.dumps(label)} , shape: "circle" }});\n')
				# f'g.setNode("{repr(q)}", {{ label: {json.dumps(hash(label) // 1e8)} , shape: "circle" }});\n')
			ret.append(f'g.node("{repr(q)}").style = "fill: #fc8d62"; \n')

		for q in self.Q:
			to = defaultdict(list)
			for a, j, w in self.arcs(q):
				if self.R is ProductSemiring and isinstance(w.score[0], String):
					# the imporant special case of encoding transducers
					label = f'{str(a)}:{str(w)}'
				else:
					label = f'{str(a)} / {str(w)}'
				to[j].append(label)

			for dest, values in to.items():
				if len(values) > 4:
					values = values[0:3] + ['. . .']
				label = '\n'.join(values)
				ret.append(
					f'g.setEdge("{repr(q)}", "{repr(dest)}", {{ arrowhead: "vee", label: {json.dumps(label)} }});\n')

		# if the machine is too big, do not attempt to make the web browser display it
		# otherwise it ends up crashing and stuff...
		if len(ret) > 256:
			return f'FST too large to draw graphic, use fst.ascii_visualize()<br /><code>FST(num_states={self.num_states})</code>'

		ret2 = ['''
		<script>
		try {
		require.config({
		paths: {
		"d3": "https://cdnjs.cloudflare.com/ajax/libs/d3/4.13.0/d3",
		"dagreD3": "https://cdnjs.cloudflare.com/ajax/libs/dagre-d3/0.6.1/dagre-d3.min"
		}
		});
		} catch {
		  ["https://cdnjs.cloudflare.com/ajax/libs/d3/4.13.0/d3.js",
		   "https://cdnjs.cloudflare.com/ajax/libs/dagre-d3/0.6.1/dagre-d3.min.js"].forEach(function (src) {
			var tag = document.createElement('script');
			tag.src = src;
			document.body.appendChild(tag);
		  })
		}
		try {
		requirejs(['d3', 'dagreD3'], function() {});
		} catch (e) {}
		try {
		require(['d3', 'dagreD3'], function() {});
		} catch (e) {}
		</script>
		<style>
		.node rect,
		.node circle,
		.node ellipse {
		stroke: #333;
		fill: #fff;
		stroke-width: 1px;
		}

		.edgePath path {
		stroke: #333;
		fill: #333;
		stroke-width: 1.5px;
		}
		</style>
		''']

		obj = 'fst_' + uuid4().hex
		ret2.append(
			f'<center><svg width="850" height="600" id="{obj}"><g/></svg></center>')
		ret2.append('''
		<script>
		(function render_d3() {
		var d3, dagreD3;
		try { // requirejs is broken on external domains
		  d3 = require('d3');
		  dagreD3 = require('dagreD3');
		} catch (e) {
		  // for google colab
		  if(typeof window.d3 !== "undefined" && typeof window.dagreD3 !== "undefined") {
			d3 = window.d3;
			dagreD3 = window.dagreD3;
		  } else { // not loaded yet, so wait and try again
			setTimeout(render_d3, 50);
			return;
		  }
		}
		//alert("loaded");
		var g = new dagreD3.graphlib.Graph().setGraph({ 'rankdir': 'LR' });
		''')
		ret2.append(''.join(ret))

		ret2.append(f'var svg = d3.select("#{obj}"); \n')
		ret2.append(f'''
		var inner = svg.select("g");

		// Set up zoom support
		var zoom = d3.zoom().scaleExtent([0.3, 5]).on("zoom", function() {{
		inner.attr("transform", d3.event.transform);
		}});
		svg.call(zoom);

		// Create the renderer
		var render = new dagreD3.render();

		// Run the renderer. This is what draws the final graph.
		render(inner, g);

		// Center the graph
		var initialScale = 0.75;
		svg.call(zoom.transform, d3.zoomIdentity.translate(
		    (svg.attr("width") - g.graph().width * initialScale) / 2, 20).scale(initialScale));

		svg.attr('height', g.graph().height * initialScale + 50);
		}})();

		</script>
		''')

		return ''.join(ret2)


if __name__ == '__main__':
	from rayuela.base.semiring import Tropical
	fsa = FSA(Tropical)
	fsa.set_I(State(0))
	fsa.add_arc(State(0), Sym('a'), State(1), Tropical(1))
	fsa.add_arc(State(0), Sym('a'), State(2), Tropical(2))
	fsa.add_arc(State(1), Sym('b'), State(1), Tropical(3))
	fsa.add_arc(State(2), Sym('b'), State(2), Tropical(3))
	fsa.add_arc(State(1), Sym('c'), State(3), Tropical(5))
	fsa.add_arc(State(2), Sym('d'), State(3), Tropical(6))
	fsa.add_arc(State(2), Sym('d'), State(1), Tropical(-2))
	# fsa.add_arc(State(3), Sym('d'), State(0), Tropical(-6))
	fsa.set_F(State(3))
	ps = Pathsum(fsa)
	print(dict(ps.dijkstra_fwd()), '\n', dict(ps.bellmanford_fwd()))

	lps = fsa.pathsum(Strategy.LEHMANN)
	jps = fsa.pathsum(Strategy.JOHNSON)
	print(lps, jps)
	assert lps == jps

	# det = fsa.determinize()
	# assert fsa.pathsum() == det.pathsum()
