from collections import defaultdict as dd
from itertools import chain, product
from sys import float_repr_style
from frozendict import frozendict

from rayuela.base.symbol import ε
from rayuela.fsa.state import MinimizeState, PowerState
from rayuela.fsa.pathsum import Pathsum, Strategy
from rayuela.fsa.push import push_with_potential


class Transformer:
    @staticmethod
    def trim(fsa):
        raise NotImplementedError

    @staticmethod
    def _powerarcs(fsa, Q):
        """ This helper method group outgoing arcs for determinization. """

        symbol2arcs, unnormalized_residuals = dd(set), fsa.R.chart()

        for q, old_residual in Q.residuals.items():
            for a, p, w in fsa.arcs(q):
                symbol2arcs[a].add(p)
                unnormalized_residuals[(a, p)] += old_residual * w

        for a, ps in symbol2arcs.items():
            normalizer = sum([unnormalized_residuals[(a, p)] for p in ps], start=fsa.R.zero)
            residuals = {p : ~normalizer * unnormalized_residuals[(a, p)] for p in ps}

            yield a, PowerState(residuals), normalizer

    @staticmethod
    def push(fsa):
        from rayuela.fsa.pathsum import Strategy
        W = Pathsum(fsa).backward(Strategy.LEHMANN)
        return push_with_potential(fsa, W)

    @staticmethod
    def _eps_partition(fsa):
        """ partition fsa into two (one with eps arcs and one with all others) """

        E = fsa.spawn()
        N = fsa.spawn(keep_init=True, keep_final=True)

        for q in fsa.Q:
            E.add_state(q)
            N.add_state(q)

        for i in fsa.Q:
            for a, j, w in fsa.arcs(i):
                if a == ε:
                    E.add_arc(i, a, j, w)
                else:
                    N.add_arc(i, a, j, w)

        return N, E

    @staticmethod
    def epsremoval(fsa):

        # note that N keeps same initial and final weights
        N, E = Transformer._eps_partition(fsa)
        W = Pathsum(E).lehmann(zero=False)

        for i in fsa.Q:
            for a, j, w in fsa.arcs(i, no_eps=True):
                for k in fsa.Q:
                    N.add_arc(i, a, k, w * W[j, k])

        # additional initial states
        for i, j in product(fsa.Q, repeat=2):
            N.add_I(j, fsa.λ[i] * W[i, j])


        return N

