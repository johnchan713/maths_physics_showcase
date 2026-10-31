# Q49: the original axial tail annihilates the graph flux

Q48 (1f7db80) bounds a generic smooth projected source by repeated
integration by parts. Its high-order archive is sufficient but unnecessary
for the **original** prefast source. The reason is the exact equation in
[the initial scalar checkpoint, (3)–(4)](../core_initial_scalar_seed_budget/README.md).
This argument uses the unchanged operator, shift, width and complete real
primary. No later fast/slow difference is replaced by a zero-mean model.

## Identity before estimation

In U=R^d write V=partial_U+M v_r.grad_y. At fixed Z,T on one band,

    D_e F = d R^(d-1-e) V(R^e F).

A_e integrates f_e times its input over the WHOLE shifted radial line,
then multiplies by k_e chi_U (Q48); f_e=R^(e-d+1)/d. Thus for a smooth
compact flux F in the positive radial shell the integrand for A_e D_e F
is precisely d/da [(R^e F)(U+a,y+Ma v_r)]. Extend the weighted flux
by zero before translating. Both endpoint traces vanish, so A_e D_e F=0
pointwise, including on the moving cutoff. This left annihilation is
distinct from D_e T_e=I-A_e. It requires compact flux traces and the same
shift in BOTH operators; it does not follow from zero Haar alone.

For the original zero-mean velocities the axial source is
Qcirc[D1(Wrz-Sigma_z)+epsilon partial_Z(Wzz+pm)]. Qcirc commutes
with the graph derivative and shifted maps; Sigma_z has no auxiliary
dependence. The complete primary covariance has smooth compact radial
extension under the original support/edge premises. Consequently

    A1 F0_z = epsilon A1 partial_Z [Qcirc Wzz + T0 Qcirc gr].

Here pm=T0(gr-rho P) before projection, gr=-D1 Wrr-epsilon partial_Z
Wrz+Wtt/R. Only the auxiliary-independent compensator disappears.
The derivative acts on the entire moving pressure, including its cutoff;
it is NOT commuted through T0. The displayed remainder is generally
nonzero. In particular the argument does not set the initial fast axial
velocity correction to zero, or eliminate Q47's full pressure differences.

## Weighted inequality and interface

Use Q48's original hull recovered from the initial budget. Let n be the
requested row order, N=n+3, v=12N+18, D=n+4, W the original normalized
R covariance polynomial at N (physical order one). With the same g,r
graph and reciprocal bounds as Q48, gr through n+1 is bounded at order
1-kappa by (Dg+D+2r)W: D pays the one graph/source derivative and g
the frozen graph coefficients. All covariance entries come from the
same original primary. Put Pm=KT0(n+1,v)(Dg+D+2r)W using the full
shifted weighted pressure including both coordinate conversions.
Taking one normalized Z derivative costs at most n+1 <= D.
Using epsilon<=1 and allowing the common weaker order 1-kappa for Wzz,
the input to A1, including the exterior epsilon, has coefficient
2D(W+Pm) at order **2-kappa**. The two pays for auxiliary projection.

For the output use Q48's fixed-gap transition argument without any
Fourier inversion. Differentiating at fixed shift and applying the raw
product bound three times gives

    Kraw = 2^(3n) (Umax-Umin) N_n(k1) N_n(chi_U) N_n(f1).
    K = Greverse(n) Z_ax Kraw E(v) Fforward(n).

Here E unweights the original source, Z_ax restores zeta only on the
original axial transition, and the Q43 factors convert the normalized R
input and output. All derivative supports of chi_U stay there. Thus

    ||A1 F0_z||_(R,n,w0) <= epsilon^(2-kappa) 2 D K (W+Pm).

The weaker Q48 rate 2-2kappa and Q40 row rate 2-4kappa follow on the
same accepted domain. There is no inverse graph-floor factor, extra
temporal inverse or large integration count. The initial pressure shape
is needed through n+2 for the existing KT interface, axial shape through
n+1. For postslow order m, the original tail alone now needs covariance
through m+13 rather than m+200012, and pressure shape through m+12.
This does not reduce Q47's separate m+14 current-map shape request or
the other initial/signed-state ledgers. The factory supplies the tail
directly from the same initial budget to the unchanged Q47 rows. Its
extra factor two on a supplied tail remains a safe overestimate.

Next audit the lower original primary and fixed-profile inputs on the
remaining finite ledgers and bound the full postslow rows, including
their strict rate margins. No unconditional physical coefficient bound,
complete residual gain, infinite realization or blow-up is established.
The exact finite-mode checks below diagnose the identity; the continuum
argument is the compact boundary evaluation above.
