# -*- coding: utf-8 -*-
"""
PHYSICS anchor for the closed-form off-shell reduction (main11 methodology).

Every intermediate reduction check (p_r = rational*sqrt(S); int p_r dr = sum a_k U_k; kernel
= A/sqrt(S); IBP into W_jk) is math-vs-math -- correct but self-referential, exactly the chain
that in the -alpha*P_r bug was internally consistent yet WRONG. The only check that matters is
against the PHYSICS: the true canonical flow of the original Hamiltonian.

This module builds delta phi from the CLOSED FORM Delta S = sum_k a_k(M,a,E,J) U_k (+third kind)
and verifies frozen + eps*delta phi against the TRUE original-variable flow: slope ~ 2, i.e. the
closed-form coefficients reproduce the real dynamics to O(eps^2), not merely the previous
quadrature. a_3=E^2, a_2=2E^2 M, a_1=4E^2 M^2-(E^2-1)J^2, a_0=2M[4E^2 M^2-J^2+2(1-E^2)Ja].
"""

