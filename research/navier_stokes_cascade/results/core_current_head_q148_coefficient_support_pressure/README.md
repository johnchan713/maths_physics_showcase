# Q148: finalized-profile coefficient supports and reserved pressure

Parent: Q147, `45efb625b15fd4039b21086895a4c2fcecf1986f`.
The live branch, Q147 evidence and CI, the frozen coefficient-support recheck,
the positive-order repair, support-saturated tuple, retained Borel schedule and
physical pulse transfer were read first. Historical evidence, roots, supports,
pressure datum, scales, frequency and the `C^-20` width remain unchanged.

Q147 supplies the exact same-profile leading input which the earlier finite
coefficient theorem required: the accepted compact repair, all five global
functional identities and the unchanged pressure gauge. This checkpoint binds
that theorem to the finalized profile. It also states the exact inequality
which still prevents the fixed-order result from becoming a uniform
all-physical-scale pulse certificate.

## 1. The Q147 profile is the coefficient theorem's input

Q147 proves, as identities of `eta`,

\[
 M_\infty=R_\infty=J_\infty=S_\infty=0,
 \qquad \Pi_{\rm axis}+C_{p,\infty}=0.                 \tag{1}
\]

Here `R` is the renormalized leading heat moment. Positive coefficient
orders use the ordinary compact rows. The accepted leading profile, its axis
datum and its future heat tail are unchanged. Thus the hypotheses of the
frozen finite coefficient-support theorem now refer to Q147's actual final
profile; no support end or local zero is used to infer (1), and pressure is
not reset at an interface.

For every fixed `n>=1`, finalize all lower orders before forming the five
ordinary discrepancies. Current-order fields occur only in `(0,n)` and
`(n,0)`; their squares have degree `2n`. With the unchanged I3 plateau, the
five-row correction is the linear solve

\[
 c_n=-B^{-1}z_n,
 \qquad z_n=(M,(M-J)/\lambda,I,-S,C_p),                \tag{2}
\]

after the documented normalizations. The rational positive-mass determinant
argument gives `||B^-1||<800`; this is not a sampled inverse.

If `mu_i` are the five factorial angular discrepancy norms,
`K=(sqrt(2)rho_3^(3/2))^-1` and

\[
 W_k=e_*^{-1}\begin{cases}2,&k=0,\\4,&k=1,\\5,&k\ge2,
 \end{cases}
\]

then the actual fixed-order estimate is

\[
 |c_n|_k\le800\max\left\{
 {\mu_1\over\rho_3},
 {{\mu_1/\rho_3+KW_k\mu_4}\over\lambda},
 K\mu_2,{W_k\mu_5\over\rho_3},W_k\mu_3\right\}.     \tag{3}
\]

Equation (3) retains the reciprocal amplitude jets and the raw
`lambda^-1` cost. Together with the inherited finite inner Volterra lemma,
it gives finite mixed seminorms for each selected `n,k`; it is not uniform in
coefficient order, derivative order or physical scale.

## 2. Exact positive-order support ledger

The first total makes `M_n=0` after `Xvel`, including every available angular
jet. The flux identity then makes the whole positive-order velocity zero
there. The potential-cutoff commutator is proportional to `M_n`, so it also
vanishes on the reserved I4 interval. The conservative totals retain the
previous-order axial viscosity and give the same exact right edges:

| Order | Velocity | Pressure | Stress |
|---|---|---|---|
| `1` | `Xvel` | `Xv` | `Xb` |
| `2` | `Xvel` | `Xv` | `Xv` |
| `n>=3` | `Xvel` | `Xvel` | `Xvel` |

Thus the reserved mean velocity is exactly the leading velocity, even in
potential-cutoff transitions. Pressure is reconstructed from its complete
source and is not a velocity component. On I4 it remains

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},\qquad
 \Pi_2=K_2(\eta),\qquad \Pi_n=0\quad(n\ge3).           \tag{4}
\]

Both `K1` and `K2` contain later-source memory. A local zero pressure
derivative cannot set either function to zero.

## 3. The sharp reserved-pressure obstruction survives

Let `P=q^(2A)p_B-Pi0`. For any `Xa<Xb` in I4, both unknown angular constants
cancel:

\[
 P(X_b)-P(X_a)=\chi(c_1q)q^{2h}{\mu^2\over4}
 (X_a^{-1}-X_b^{-1}).                                  \tag{5}
\]

Therefore

\[
 \sup_{X\in I4}|P(X)|\ge
 { |\chi(c_1q)|q^{2h}\mu^2\over8\rho_3}
 (e^{-6}-e^{-11}).                                     \tag{6}
\]

The factor `1/8` is sharp even if arbitrary X-independent cancellation is
allowed. This is an analytical two-point inequality, independent of `K1`
and `K2`. Exact-rational diagnostics replay the algebra and exponential
enclosures; they do not evaluate a continuum supremum. Complete covariance
means must retain (4), its physical time/axial derivatives, leading radial
velocity, full-wave axial velocity and curl remainders.

## 4. Precise inequality blocking uniform physical pulses

The fixed-order estimate (3) alone supplies no bound on the actual full-hull
tuple seminorms `M_(n,m)` as `n,m` vary. The retained Borel construction
requires, for the same implemented cutoff and schedule,

\[
 C_{n,m}=\max\{1,2K_{m+1}(1+2nh)^mM_{n,m}\},\qquad
 C_{n,m}a_n^{-nh}\le2^{-n},\quad m\le n,               \tag{7}
\]

at every retained order. Q147 closes the leading exact identities needed by
the induction, but neither (1) nor fixed-order finiteness proves (7). In
particular, no supremum in `n,m`, common inverse/gap domain or all-scale
support estimate follows from (3). This is the first unresolved inequality
before the already audited primary-pulse transfer can be applied to this
same finalized background.

The finite positive-order coefficient and reserved-mean support obligation is
therefore complete on Q147's profile. Uniform actual-pulse bounds remain
conditional on (7) and the associated same-profile inverse, gap and support
premises. Complete covariance means, the full viscosity/cutoff/all-interaction
residual and blow-up remain unproved.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json` here.
The frozen record binds both Q147 and the historical support evidence without
editing either one.
