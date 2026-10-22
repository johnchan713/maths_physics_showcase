# Higher-pole operators and fresh post-sourced mean size

Parent: `702886d7be7e46fc038e544987a1dbb25e5d1417`. This checkpoint
closes the parent's **finite weighted-operator domain gap and the next
fresh mean SIZE composition**, conditionally on the original retained
profile, field, and support bounds. It does not close the next nonzero
residual. Historical frozen evidence is unchanged.

## 1. A fixed higher-pole operator constant, with the actual gap

Use the original shell length L, width w=C_setup^-20, and actual cutoff
plateau gap eta. No geometry is changed. Set a=w^2, b=4,
r=min(1,L/4,eta), and delta=min(1,y_a,y_b). For finite pole q, the
same two-edge primitive proof uses

    d=max(0,ceil((q-3)/2)),
    D(c)=(1/(2c))*sum_{j=0}^d binom(d,j) j! (r^2/c)^j,
    C_q=max{D(a) exp(b/(L-r)^2),
            D(b) exp(a/(L-r)^2),
            L E_q exp((a+b)/r^2)}.                     (1)

Here E_q is the existing exact edge absorption constant satisfying
sqrt(zeta)*delta^-q <= E_q. In the inherited construction it can be
chosen max(1,N!/(w^2/2)^N), N=ceil(q/2). Formula (1) is an upper
bound, not the optimal norm. All fixed differentiated multiplier,
cutoff, Jacobian, and chart factors remain in the original profile
envelopes at each requested jet order. Both one-sided kernels remain.

The **specific retained positive recipes** in
`core_initial_haar_row_budget` and `core_signed_prefix_source_budget`
give, at each n,

    KT_old >= L1,   KA_old >= Z*LA,
    Kp_old >= L0,   Kmom_old >= LM.                    (2)

L1,L0,LM are their pole-independent positive profile sums, including the
radial integration weights for LM; Z*LA includes the fixed interior
weight floor and actual A_e output profile. Indeed their old factors
1+C_2p, E_2p, the sum of pressure factors at 2p and 4p, and the sum of
moment edge factors at 2p and 4p are all at least one. Consequently the
SAME maps at the new proof domains obey

    KT_new(6p)       <= KT_old*(1+C_6p),
    KA_new(6p)       <= KA_old*E_6p,
    Kp_new(6p,12p)   <= Kp_old*(2+C_6p+C_12p),
    Kmom_new(6p,12p) <= Kmom_old*(E_6p+E_12p).         (3)

Equation (2) is a consequence of the retained recipe construction. It
is **not** a consequence of knowing an arbitrary old operator norm.
Callbacks supplied without that provenance do not certify (2). This
checkpoint derives (3) under that existing construction contract; it
does not validate the original actual profile enclosures. Interior
slow outputs keep their existing weight-floor/profile bound: delta^q
<=1, so no new shell-edge singularity is introduced on their closed
reserved hull. The pressure bump and full integration hull are unchanged.

Define a single finite rational constant by keeping its exact power recipe:

    K=3+E_6p+E_12p
      +sum_{q in {6p,12p}} sum_{three branches of (1)}
           coefficient(q,branch)*3^ceil(exponent(q,branch)).  (4)

Since exp(t)<=3^ceil(t) for t>=0, K dominates all four multipliers
in (3). It is fixed after original p,w,L,eta are fixed. The formula is
valid even when the selected width makes materialization impractical.
The diagnostic exponent guard in Python is not an analytical gate.
There is no claim of uniformity as width or gap shrinks, or as n or p
varies without bound.

## 2. Positive polynomials without a new accepted-band condition

Write x=S+K **only as an estimate argument**. The actual physical S,
accepted band, carrier, length, pulse supports, and cutoffs remain as
before. Every old positive polynomial A satisfies A(S)<=A(x), K<=x.
Thus each new operator polynomial can be bounded by x*A_old(x).
All unchanged input estimates can also be evaluated at x. Composition
with the existing positive-polynomial routines gives a finite explicit
polynomial F_n(S+K) in the original S. Equivalently,

    [S^j] F_n(S+K)=sum_{i>=j} f_i binom(i,j) K^(i-j).   (5)

This is an algebraic majorant; no physical call at scale x is made.
In particular S>=K is unnecessary. This avoids a new smallness or
starting-band gate and does not change the original inverse length proof.
For reproducible arithmetic the current-wave and accumulated-state polynomials are further
bounded by its sparse hull: sum c_j*x^j <= (sum c_j)*x^d for x>=1,
d=max{j:c_j>0}. This exact positive inequality retains every term's
coefficient mass, and never samples or removes a physical interaction.

## 3. Fresh fast, fresh slow, then pressure once again

Let B be the parent's original `before` prefix and Delta U the complete
generated sourced increment, with all signed modes, same cutoff curl,
and full unprojected ODE pressure retained by the parent. The next mean
calculation starts from

    M_initial=B.final,   U_current<=hull(B.wave.U+Delta U). (6)

The entire U_current is held fixed during these mean maps. Its covariance
majorant is U_current^2; (6) retains old-wave crosses and increment
squares, including the second harmonic. We do not replace (6) by the
increment alone, nor reset the old slow means to zero.

Apply the original complete mean equations and map formulas (3)--(6)
of `core_signed_prefix_source_budget`, replacing only the four operator
envelopes as in (3)--(5):

1. Reconstruct full gr, P, pm, Etheta, Ez, Jtheta, Jz using (6).
2. Add the fresh fast velocities, including actual gd-A1(gd) in the
   axial construction and its KA1 cost.
3. Reconstruct full gr, P, pm and all three scalars on that postfast
   prefix; apply the slow map to these freshly reconstructed quantities.
4. Add the fresh slow velocities to the accumulated state, then
   reconstruct full pressure and rows again on the final prefix.

The pressure memory is the full conservative source, compensated by
the actual bump and full current scalar; the old A0/A1 velocities,
complete covariance, viscosity, graph and cylindrical terms remain
inside its majorant. The inherited moving-compensator identities are
not replaced by selected transverse fluxes. The parent's exact tensor
and full moment controls are separately replayed.

For each fixed n this supplies conditional explicit finite-polynomial
SIZE bounds for the next actual mean construction: tangential/axial
order mu=1-kappa, radial order mu+1, sufficient final edge pole 6p.
The full final mean rows/pressure may require pole 12p, which (3)
covers. These are the same coarse physical orders used by the inherited
size closure; **no scalar-defect or full residual improvement follows**.
The order-zero final-row composition reads original pre-source mean
jets through order 31. Finiteness of the composed loss is distinct from
validity of those original input enclosures.

## 4. Scope and the next inequality

Removed independent gap: higher-pole weighted operator SIZE coverage
and fresh post-wave fast/slow mean SIZE are now derived from the original
recipe-defined envelopes and generated complete wave. Added analytical
assumptions: none relative to that retained contract. Still conditional:
valid original primitive/old-wave/source bounds, actual recipe provenance,
common-chart support and coefficients, and homogeneous-frame/root
identification where needed. A manufactured fixture is not evidence
that these premises hold for every physical scale.

The inherited nonzero `source_recipe` explicitly treats one real
frequency with signed harmonics +/-1. It cannot be reused unchanged
after (6), whose harmonics are +/-1,+/-2. To continue, derive and bound
the exact regular source at all output modes +/-1,...,+/-4 (plus the
fresh zero-mode equations) using the new means, full reconstructed
pressure, actual angular factors |m|<=2 and m^2<=4, viscosity, cutoff
tails and every interaction. In particular the needed inequality is

    ||F_next,m||_{n, regular weighted}
       <= epsilon^(proposed improved exponent) Q_{n,m}(S)

on the unchanged physical band, after identifying the exact principal
cancellations and keeping the complete pressure-gradient terms. Neither
the improved exponent nor this inequality is proved here. Generic
positive SIZE polynomials do not supply it. The new class therefore
rejects `source`, `raw_source`, and `inverse_source` rather than forwarding
the wrong frequency grammar. There is no unconditional pulse theorem,
complete residual gain, or blow-up result.

## 5. Evidence

`controls.py` compares (1) to the existing independent materializer at
both new poles (including zero and fractional p), checks exact finite
shift domination and binomial translation even at S<K, and composes the
fresh full-pressure maps on a manufactured original-recipe fixture.
It tests stale-pressure/zero-reset exclusions, A1 retention, the finite
jet dependency, and rejection of the old next-source interface.
These exact rational/polynomial diagnostics check implementation and
proof bookkeeping, not continuum enclosure or blow-up. `audit.py`
binds historical evidence and declared sources; output must be fresh.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
