# Q205: stagewise generated pulses and the diagonal quantifier

Parent: Q204, `fd91fff61a1ce514283afbd267870883c14d352c`.
The live branch, Q204, its active CI, the sourced-inverse foundation, the
physical source/tail join, the certified complete-potential diagonal and the
Q203 endpoint/localization replay were read first. Historical evidence, the
finalized physical family, pressure gauge, roots, supports, accepted bands,
cutoff sequence and unchanged `C^-20` width are preserved.

Q204 correctly records that no constant uniform over the complete generated
hierarchy has been proved. It leaves ambiguous whether such a constant is
needed for the already invoked diagonal. This checkpoint resolves that
quantifier. Every fixed generated cycle has an all-scale value, derivative,
pressure and closed-support estimate. Its constants may grow arbitrarily with
the cycle. The certified diagonal was designed for exactly this situation: it
chooses the cycle cutoff only after the finite current constants are known.
Thus a hierarchy-uniform constant is not a missing premise of Q203's
conditional locally finite residual conclusion. It remains a genuine missing
premise for a stronger common bound on the uncut raw generated pulses.

All continuum conclusions below remain conditional on the named provider
estimates. The exact inequalities, countermodel and evidence hashes check the
quantifiers and composition; they do not prove those provider estimates.

## 1. Fixed-cycle actual pulse bounds

For a fixed finite harmonic set, cycle `j` and physical derivative counts
`a,b`, the sourced inverse gives

\[
 |D^{a,b}\delta u_j|\le C_{j,a,b}
 Q^{-A+h\alpha_j-a s_x-b s_t},                       \tag{1}
\]
\[
 |D^{a,b}\delta p_j|\le C'_{j,a,b}
 Q^{-2A+h(\alpha_j+1/2)-a s_x-b s_t},                \tag{2}
\]
\[
 |D^{a,b}\delta\mathcal A_j|\le C''_{j,a,b}
 Q^{1/2-A+h(\alpha_j+1/2)-a s_x-b s_t}.              \tag{3}
\]

The constants include the actual source, frame, graph and edge costs. They
are independent of the running physical band inside that fixed provider, but
may depend on `j,a,b` and the fixed profiles. The pressure rate in (2) cannot
be improved to the potential rate: for constant normal `n=e3` and a normal
source `f=e3 F(v)`, the exact solution has zero tangent response but
`pi=iF/(km)`. Both moving-normal and source-normal pressure terms therefore
remain.

Full-path zero data give zero response and differentiated zero extension on
the transverse, slow and shell boundaries. A compact source can leave a
nonzero exit, so the original cutoff must be applied to both potential and
pressure before curl. Curl, graph and frame derivatives do not enlarge the
resulting closed supports. These are actual fixed-cycle all-band estimates,
not a single estimate uniform over every generated cycle and not absolute
decay as `Q` tends to zero.

The finite source/tail recurrence separately yields, for every fixed
`j,m,N`,

\[
 |D_{phys}^m(T^{[j]}_{compact}+T^{[j]}_{Gaussian}
                 +T_{base/exterior})|
 \le C_{j,m,N}q^N.                                   \tag{4}
\]

It retains all six covariance entries, radial and axial velocity, the whole
moving pressure product, `Pi1/Pi2`, both same-base pressure tails, viscosity,
cutoffs and all interactions. Equation (4) is stagewise; it does not assert
`sup_j C_{j,m,N}<infinity`.

## 2. Precise obstruction to a common raw-hierarchy constant

The available premise has the quantifier order

\[
 \forall j\;\forall m\;\exists C_{j,m}<\infty.       \tag{5}
\]

It does not imply

\[
 \forall m\;\exists C_m<\infty\;\forall j.           \tag{6}
\]

This is not only a verbal gap. The representative estimate before the scale
cut is

\[
 |(A_j,b_je_\theta,\pi_j)|_m
 \le C_{j,m}(1+|\log q|)^{P_{j,m}}
 q^{h\delta j-L_m},\qquad
 \delta={1\over2}-3\kappa.                           \tag{7}
\]

The abstract scalar representatives

\[
 X_j(q)=2^{j^2}q^{h\delta j-L_m}
\]

with `C_{j,m}=2^(j^2)`, `P_{j,m}=0` satisfy exactly the stagewise premise
(5)--(7). At the fixed physical scale `q0=1/2`, their normalized size is

\[
 2^{j^2-h\delta j},                                  \tag{8}
\]

which is unbounded. Thus the displayed stagewise inequalities alone have a
model in which (6) fails. Neither a common starting band nor a uniform raw
generated-pulse constant follows from those premises. This is an obstruction
to that stronger inference, not a lower bound for the actual pulses or
residual.

## 3. Why the locally finite diagonal does not need (6)

For the fixed-gauge increment

\[
 \Delta u_j=\operatorname{curl}A_j+b_je_\theta,
 \qquad \pi_j=p^{[j]}-p^{[j-1]},                     \tag{9}
\]

the certified cut is

\[
 V_{j,a}=\operatorname{curl}(\chi(aq)A_j)
          +\chi(aq)b_je_\theta,
 \qquad P_{j,a}=\chi(aq)\pi_j.                      \tag{10}
\]

Stage-dependent potential or pressure gauges are forbidden because cutting
them creates new non-telescoping terms. After reserving one curl, (10) has

\[
 |(V_{j,a},P_{j,a})|_m
 \le \widetilde C_{j,m}(1+|\log q|)^{\widetilde P_{j,m}}
 q^{h\delta j-L_m},                                  \tag{11}
\]

uniformly in `a` on its transition. At cycle `j`, collect only the finitely
many requests `m<=j` and put `gamma_j=h delta j/2`. Once their finite
constants are known, choose `s_j=log a_j` above

\[
 \max\left\{s_{j-1}+1,\log(q_*^{-1})+1,1,
 {2(c+j)\over\gamma_j}+{P^2\over\gamma_j^2}\right\} \tag{12}
\]

for every current `c>=log C` and `P>=0`. Then, on the whole interval
`0<q<=a_j^-1`,

\[
 C(1+|\log q|)^Pq^{\gamma_j}\le2^{-j}.               \tag{13}
\]

Equation (12) exists for every finite list regardless of its growth. With
strictly increasing `a_j`, the sum of (10) is locally finite on `q>0`.
For a fixed jet and requested power one then fixes a sufficiently late
comparison cycle and lets `q` tend to zero. No uncontrolled substitution
`j=j(q)` and no `sup_j C_{j,m}` occur.

Therefore Q203's six-group residual and twelve-term localization retain time,
unit viscosity, complete pressure, both nonlinear cross orientations, the
square and every scale-cutoff derivative under the same conditional
providers. Stage-dependent constants suffice for that diagonal. They do not
establish a common raw-pulse norm.

## 4. Remaining limits and next obligation

This checkpoint corrects Q204's frontier but does not strengthen its pulse
decay. The primary and current sourced zero-jet upper-bound exponents remain
negative, the reserved-pressure floor remains positive, and the generated
fixed-cycle estimates remain conditional. No uniform uncut hierarchy bound,
absolute pulse smallness, unconditional complete residual improvement,
regular terminal velocity or pressure, competitor exclusion or blow-up is
proved.

The next obligation is to audit the same locally finite object at the
terminal origin: establish or obstruct one-sided velocity and pressure jets
in the fixed gauge, rather than importing the already closed smooth force
trace. That audit must retain the stagewise cutoff sequence, complete pressure
memory, axial velocity, all covariance means and the protected core bound.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q205_stagewise_pulse_diagonal/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q205_stagewise_pulse_diagonal/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q205_stagewise_pulse_diagonal/evidence.json \
  --output fresh.json
```
