class SCC:

    def __init__(self, fsa):
        self.fsa = fsa

    def scc(self):
        """
        Computes the SCCs of the FSA.
        Currently uses Kosaraju's algorithm.

        Guarantees SCCs come back in topological order.
        """
        return self._kosaraju()

    def _kosaraju(self) -> "list[frozenset]":
        """
        Kosaraju's algorithm [https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm]
        Runs in O(E + V) time.
        Returns in the SCCs in topologically sorted order.
        """
        # Homework 3: Question 4
        visited = set()
        stack = []

        def dfs(u):
            if u in visited:
                return
            visited.add(u)
            for a, v, w in self.fsa.arcs(u):
                dfs(v)
            stack.append(u)

        for q in self.fsa.Q:
            dfs(q)

        component = {}  # state to its component
        rev_fsa = self.fsa.reverse()

        def rev_dfs(u, root):
            if u in component:
                return
            component[u] = root
            for a, v, w in rev_fsa.arcs(u):
                rev_dfs(v, root)

        while len(stack) > 0:
            q = stack.pop()
            rev_dfs(q, q)

        sccs = {root: set() for root in component.values()}  # component to states
        for q, root in component.items():
            sccs[root].add(q)

        # toposort on sccs
        G = {c: [] for c in sccs.keys()}  # adjacency graph of components
        in_deg = {c: 0 for c in sccs.keys()}
        for u in self.fsa.Q:
            for a, v, w in self.fsa.arcs(u):
                cu, cv = component[u], component[v]
                if cu != cv:
                    G[cu].append(cv)
                    in_deg[cv] += 1

        stack = [c for c, d in in_deg.items() if d == 0]  # zero indeg components
        result = []
        while len(stack) > 0:
            cu = stack.pop()
            result.append(sccs[cu])
            for cv in G[cu]:
                in_deg[cv] -= 1
                if in_deg[cv] == 0:
                    stack.append(cv)

        return result
