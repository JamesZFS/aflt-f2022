import json
import frozendict


class State:

    def __init__(self, idx, label=None):
        self._idx = idx
        self._label = label

    @property
    def idx(self):
        return self._idx

    @property
    def label(self):
        return self._label

    def set_label(self, label):
        self._label = label

    def copy(self):
        return State(self.idx)

    def __repr__(self):
        if self.label is not None:
            return f'{self.label}'
        return f'{self.idx}'

    def __str__(self):
        if self.label is not None:
            return f'{self.label}'
        return str(self.idx)

    def __hash__(self):
        return hash(self.idx)

    def __eq__(self, other):
        return isinstance(other, State) and self.idx == other.idx


class PowerState(State):
	"""
	A state that is an element of the powerset of states.
	This is useful for determinization.
	"""

	def __init__(self, residuals):
		super().__init__(frozenset({p for p, _ in residuals.items()}))
		self.residuals = residuals

	def __repr__(self):
		return "PowerState(" + str(self) + ")"

	def __str__(self):
		if self.label is not None:
			return f'{self.label}'
		contents = []
		if self.residuals is None:
			for state in self.idx:
				contents.append(f"{state.idx}")
		else:
			for state in self.idx:
				contents.append(f"{state.idx}/{str(self.residuals[state])}")

		return "{" + ", ".join(contents) + "}"

	def __hash__(self):
		return hash((self.idx, frozendict(self.residuals)))


class MinimizeState(State):
	"""
	A state that is an element of the powerset of states.
	This is useful for minimization.
	"""

	def __init__(self, states):
		super().__init__(frozenset(states))

	def __repr__(self):
		if self.label is not None:
			return f'{self.label}'
		contents = []
		for state in self.idx:
			contents.append(f"{state}")
		return "{" + ",".join(contents) + "}"

	def __str__(self):
		return self.__repr__()

	def __hash__(self):
		return hash(self.idx)


class PairState(State):

    def __init__(self, p, q):
        super().__init__((p, q))

    @property
    def state1(self):
        return self.idx[0]

    @property
    def state2(self):
        return self.idx[1]

    def __repr__(self):
        return f"( {str(self.state1)}, {str(self.state2)} )"

    def __str__(self):
        return self.__repr__()

    def __iter__(self):
        return iter((self.idx[0], self.idx[1]))


class PowerState(State):
    def __init__(self, residual_map: dict, label=None):
        super(PowerState, self).__init__(dict(residual_map), label)

    @property
    def states(self):
        for q in self._idx.keys():
            yield q

    def __iter__(self):
        return iter(self._idx.items())

    def residual(self, q):
        return self._idx[q]

    def __repr__(self):
        if self.label is not None:
            return f'{self.label}'
        return "{" + f"{[(p, w.score) for p, w in self._idx.items()]}"[1:-1] + "}"

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash(frozenset((q, w) for q, w in self._idx.items()))


if __name__ == '__main__':
    d = {1: 0.1, 2: 0.3}
    s = PowerState(d)
    for q, w in s:
        print(q, w)
