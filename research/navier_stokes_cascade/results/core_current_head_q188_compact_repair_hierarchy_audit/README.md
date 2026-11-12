# Q188: exact compact modulation/repair hierarchy and the remaining continuum envelope

Parent: Q187, `06e77c9804e579c2d6acd18af03b0f0ac9797993`.
The live branch, Q187, frozen continuation/compact/repair/cone/source records
and current CI were read first.  Every historical parameter, root, support,
pressure gauge, frequency, force and the `C^-20` activation width is unchanged.

Q187 left the retained continuation, low-order loop/repair, leading strict
margins and regional source suprema as one broad independent-validation
frontier.  That description is accurate as a warning against unconditional
promotion, but it is stale as a provider inventory: the live repository
already contains all-fixed-jet continuation, the same-map all-order
modulation and repair, annular/exterior continuation, regional cone bounds,
global swirl/axial source bounds, the selected cutoff and explicit tilt.

This checkpoint independently redoes the exact low-order repair/frequency
arithmetic, retains all five pressure-memory rows, and joins that result to
the existing inner/compact/outer cone theorem.  It finds no remaining scalar
frequency obstruction.  It does not independently rederive the continuum
suprema of the selected profile from the initial equations.

## 1. Repair neighborhood and state estimate

Use the selected hierarchy

\[
 A=C^{4096}\ge2^{16},\qquad G=e^{A^{32}},\qquad H=e^{G^4},
\]

and the exact repair choices

\[
 C_{state}=D=H^{16},\quad C_{corr}=A^{128},\quad
 r=A^{-128},\quad \epsilon=G^{-64},\quad
 N=1+\lfloor H^{32}\rfloor .                         \tag{1}
\]

If the coefficient norm is `s<=r`, the actual width-`.1` bumps give

\[
 d=1280As,\qquad q=4\,10^6As .                       \tag{2}
\]

The positive-field estimate follows from
`2560<A^126`, and `1280<A^127` gives `d<1`.  The exact shear calculation
reduces to

\[
 16\,10^6+5120A^2<A^4,                               \tag{3}
\]

while the four pressure/shear coordinates obey

\[
 \|\Delta(a,b,p_1,p_2)\|_\infty
 \le327680A^9s<A^{11}s<A^{128}s .                    \tag{4}
\]

At `A=2^16`, (3)--(4) are strict integer inequalities; their slack only
increases with `A`.  No giant scale or binary floating value is materialized.

The compact cone has common gap `G^-4`, first-coordinate floor `G^-3` and
box radius at most `4G^3`.  The Lipschitz comparison reduces exactly to

\[
 G^{30}\ge512\,4^{10},                               \tag{5}
\]

which follows from `G>A` and the checked integer inequality
`A^30>512*4^10`.

## 2. The selected frequency satisfies every repair gate

The five-row inverse and quadratic constants remain

\[
 \beta=1000,\qquad q_*=10^5 .                         \tag{6}
\]

The accepted root must meet the state, coefficient-radius, contraction,
incoming-memory and modulation-transfer gates

\[
\begin{aligned}
 N&>2H^{16}(1+2000A^{128})G^{64},\\
 N&>4000H^{16}A^{128},\\
 N&>8\,10^{11}H^{16},\\
 N&>512A^6H,\\
 N&\ge4AH.                                             \tag{7}
\end{aligned}
\]

For the first inequality, its logarithm is bounded by

\[
 16G^4+\log(4002)+128\log A+64A^{32}
 <16G^4+193G<32G^4.                                  \tag{8}
\]

The executable proves (8) without floating logarithms:
`4002<2^12`, `log 2<1`, `G>A`, and `193<16A^3` suffice.
Similarly `4000<2^12` and `16*10^11<2^41` reduce the second and third
lines to `12+128A<16A^4` and `41<16A^4`.  The last two lines follow from
`512A^6<A^31` and `4A<A^31`.  Since (1) has `N>H^32`, every line of (7)
holds strictly.

Thus the same accepted root has coefficient norm below `r/2`, contraction
factor below `1/2`, and combined modulation/repair state error below
`epsilon/2`.  Equations (5) and (7) preserve the compact cone, and the
regional theorem joins it to the protected inner collar and compensated
exterior.  This is a theorem composition relative to the displayed actual
profile bounds, not a new numerical profile sample.

## 3. Complete pressure memory is part of the join

The forward rows are

\[
 (M,I,J,S,C_p),\qquad \Pi=\Pi_{axis}+C_p .             \tag{9}
\]

The loop calculation keeps the incoming error
`m=32 A^4 H/N` in every row and leaves `Pi_axis` unchanged.  The exact repair
solves the normalized five-row equation

\[
 Bc+Q(c,c)=z_N(\eta)                                  \tag{10}
\]

with the `lambda^-1` target loss and the pressure terms in `p1,p2` retained.
At the right endpoint all five rows agree with the original profile as
functions of `eta`.  Because the axis datum is the same, (9) then gives the
same pressure and integrated stress coordinates at every later radius.
Restoring four rows would leave a nonzero pressure defect and is not enough.

The all-order theorem identifies these maps with the same accepted C1 root;
it does not select a new root or frequency for each derivative order.  The
annular/exterior and regional providers use the same global moment ledger,
so the pressure memory in this checkpoint is continuous across the cone
join rather than reset at an interface.

## 4. Correction and precise remaining obstruction

The following items in Q187's provider inventory are already supplied on the
live branch, relative to their stated inputs:

* prescribed continuation through shape has all fixed jets;
* the moving phase and exact repair use one finite `N` and one accepted root
  at every fixed order;
* the annular/exterior continuation supplies the same-root global moments;
* global swirl and axial source envelopes and the signed pressure moments
  retain shifted viscosity and future pressure memory;
* the coefficient cutoff and tilt are explicit.

Consequently there is no missing low-order repair or frequency inequality
between Q187's symbolic selector and those providers.  The surviving
independent-validation boundary is the actual continuum input to this
arithmetic.  On one and the same selected full-pressure profile it is the
conjunction

\[
 \max\{\text{required mixed jets through order 2},
        \|(a,b,p_1,p_2)\|_{C^1},E^{-1},H_x^{-1},\lambda^{-1}\}
 \le A,                                                \tag{11}
\]

and, on the subsets where the last two terms are required,

\[
 \min\{a,P_c-2,P_c-v,
 2(P_c-v)^2-(v-2)J_c^2\}\ge A^{-1}.                  \tag{12}
\]

The named analytical records derive (11)--(12), and this checkpoint binds
their exact outputs.  It does not independently reconstruct every one of
their suprema from the defining ODEs and roots.  That independent
same-profile reconstruction is the next justified obligation.  Failure of
either (11) or (12) blocks the cone conclusion; passing exact-arithmetic and
source-binding tests cannot repair such a failure.

No prescribed historical first band, unconditional pulse family, complete
residual improvement, admissible force, competitor exclusion or
Navier--Stokes blow-up is proved here.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q188_compact_repair_hierarchy_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q188_compact_repair_hierarchy_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q188_compact_repair_hierarchy_audit/evidence.json \
  --output /tmp/q188_compact_repair_hierarchy.json
```
