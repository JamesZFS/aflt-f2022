def push_with_potential(fsa, V, sanity_check=True):
    """
    Mohri (2001)'s weight pushing algorithm. See Eqs 1, 2, 3.
    Link: https://www.isca-speech.org/archive_v0/archive_papers/eurospeech_2001/e01_1603.pdf.
    """

    pfsa = fsa.spawn()
    for i in fsa.Q:
        pfsa.set_I(i, fsa.λ[i] * V[i])
        pfsa.set_F(i, ~V[i] * fsa.ρ[i])
        for a, j, w in fsa.arcs(i):
            pfsa.add_arc(i, a, j, ~V[i] * w * V[j])

    if sanity_check:
        assert pfsa.pushed  # sanity check
    return pfsa