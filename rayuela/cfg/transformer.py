import itertools as it
from itertools import chain, combinations

from rayuela.base.semiring import Boolean, Real, Derivation
from rayuela.base.symbol import Sym, ε

from rayuela.fsa.state import State
from rayuela.fsa.pathsum import Pathsum

from rayuela.cfg.nonterminal import NT, S, Slash, Other
from rayuela.cfg.production import Production
from rayuela.cfg.cfg import CFG
from rayuela.cfg.treesum import Treesum

def unary(p):
    # X → Y
    if len(p.body) == 1 and isinstance(p.body[0], NT):
        return True
    return False

def preterminal(p):
    # X → a
    (head, body) = p
    if len(body) == 1 and isinstance(p.body[0], Sym):
        return True
    return False

def binarized(p):
    # X → Y Z
    (head, body) = p
    if len(body) == 2 and isinstance(p.body[0], NT) and isinstance(p.body[1], NT):
        return True
    return False


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


class Transformer:

    def __init__(self):
        self.counter = 0

    def _gen_nt(self):
        self.counter += 1
        return NT(f"@{self.counter}")

    def booleanize(self, cfg):
        one = Boolean(True)
        ncfg = CFG(R=Boolean)
        ncfg.S = cfg.S
        for p, w in cfg.P:
            if w != cfg.R.zero:
                ncfg.add(one, p.head, *p.body)
        return ncfg

    def _fold(self, cfg, p, w, I):

        # basic sanity checks
        for (i, j) in I:
            assert i >= 0 and j >= i and j < len(p.body)

        # new productions
        P, heads = [], []
        for (i, j) in I:
            head = self._gen_nt()
            heads.append(head)
            body = p.body[i:j+1]
            P.append(((head, body), cfg.R.one))

        # new "head" production
        body = tuple()
        start = 0
        for (end, n), head in zip(I, heads):
            body += p.body[start:end] + (head,)
            start = n+1
        body += p.body[start:]
        P.append(((p.head, body), w))

        return P

    def fold(self, cfg, p, w, I):
        ncfg = cfg.spawn()
        add = ncfg.add

        for (q, w) in cfg.P:
            if p != q:
                add(w, q.head, *q.body)

        for (head, body), w, in self._fold(cfg, p, w, I):
            add(w, head, *body)

        ncfg.make_unary_fsa()
        return ncfg

    def cnf(self, cfg):

        # remove terminals
        ncfg = self.separate_terminals(cfg)

        # remove nullary rules
        ncfg = self.nullaryremove(ncfg)

        # remove unary rules
        ncfg = self.unaryremove(ncfg)

        # binarize
        ncfg = self.binarize(ncfg)

        return ncfg.trim()

    def unaryremove(self, cfg) -> CFG:
        # Assignment 6
        raise NotImplementedError

    def nullaryremove(self, cfg) -> CFG:
        # Assignment 6
        cfg_null = cfg.spawn()
        for (X, body), w in cfg.P:
            assert len(body) == 1 or len(body) == 2
            if len(body) == 2 or (len(body) == 1 and body[0] == ε):
                cfg_null.add(w, X, *body)

        ts_null = Treesum(cfg_null).table()

        cfg_new = cfg.spawn()
        for p, w in cfg.P:
            (X, body) = p
            if len(body) == 2:
                Y, Z = body
                # X_null -> Y_null Z_null
                cfg_new.add(w * ts_null[Z], X, Y)  # X_notnull -> Y_notnull Z_null
                cfg_new.add(w * ts_null[Y], X, Z)  # X_notnull -> Y_null Z_notnull
                cfg_new.add(w, X, Y, Z)  # X_notnull -> Y_notnull Z_notnull
            else:
                x = body[0]
                if x != ε:
                    cfg_new.add(w, X, x)
        cfg_new.add(ts_null[S], S, ε)

        return cfg_new

    def separate_terminals(self, cfg) -> CFG:
        # Assignment 7
        raise NotImplementedError

    def binarize(self, cfg) -> CFG:
        # Assignment 7
        raise NotImplementedError

if __name__ == '__main__':
    cfg = CFG(Real)
    X = NT("X")
    Y = NT("Y")
    Z = NT("Z")
    x = Sym("x")
    y = Sym("y")
    cfg.add(Real(2), S, X, Y)
    cfg.add(Real(0.5), S, ε)
    cfg.add(Real(0.33), X, X, Z)
    cfg.add(Real(1), X, x)
    cfg.add(Real(3), X, ε)
    cfg.add(Real(2), Y, y)
    cfg.add(Real(4), Y, ε)
    cfg.add(Real(2), Z, y)
    # cfg.add(Real(2), Z, ε)

    print(cfg)
    print(cfg.treesum())
    ncfg = Transformer().nullaryremove(cfg)

    print(ncfg)
    print(ncfg.treesum())
    assert ncfg.treesum() == cfg.treesum()
