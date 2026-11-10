# Q153: physical source, tail and complete-potential diagonal join

Parent: Q152, `30cdd67580073347b99bb9f4d84a75af326cbd1e`.
The live branch, Q152 evidence and CI, the regular-source envelope and
foundation, physical-tail closures, finite-stage induction and certified
complete-potential diagonal were read first. Historical evidence, the
finalized profile, physical pulse family, roots, pressure gauge, supports,
accepted bands and `C^-20` width are unchanged.

Q152 left the same-family generated-source Gaussian/endpoint estimate and
all excluded physical tails as its next gate. Those estimates have already
been proved as conditional, family-local implications, but had not been
bound to Q149--Q152's inherited physical family. This checkpoint performs
that binding. It also corrects Q152's tail classification: spatial cutoff
terms remain in the complete regular operator, while scale cutoffs enter
only in the complete-potential diagonal. Counting either as an excluded
per-stage tail would omit or double count an exact residual term.

## 1. Full-path source quotient on the physical family

Keep Q152's same-family orders

\[
 B=1-\kappa,\qquad C=B+\tfrac12,\qquad
 D=C-\kappa=\tfrac32-2\kappa .                         \tag{1}
\]

Q149 supplies the physical primary potentials and pressures on their whole
lifted paths, with closed support and smooth zero extension. Q150 and Q152
supply the complete current covariance, pressure and means: all six tensor
entries, axial velocity, the `A1` response, both exceptional background
pressure tails and the derivative of the whole moving pressure moment are
retained. The exact physical stress identity in Q152 closes the selected
row and gives the next regular nonzero source order `D` with pole `9p`.

The inherited source-envelope theorem is local to these current inputs; it
does not require the separately rebuilt canonical entrance. For a fixed
physical label let `P>0` be its original Gaussian on the closed path. If
`K_j` bounds raw mixed jets of `log P`, define the Bell majorants

\[
 H_0=1,\qquad
 H_j=\sum_{i=1}^j {j-1\choose i-1}K_iH_{j-i}.           \tag{2}
\]

If `A_j` bounds the same raw jets of the complete source divided only by
its displayed `P` factor, Leibniz gives

\[
 G_j=\sum_{i=0}^j {j\choose i}H_iA_{j-i},\qquad
 |D^j(f_m/P)|\le
 \epsilon^D G_j(S)\sqrt\zeta\,\delta^{-q_j}.           \tag{3}
\]

Every coefficient in (2)--(3) is finite for a fixed requested order. The
pole includes the current `9p`, the finite Gaussian-log pole and every
same-domain chart loss; no lower-order archive is promoted to a higher
order. There is no division by a tiny numerical lower bound for `P`.

The complete source grammar proves the premise for (3): the linear source
keeps radial advection, slow time, pressure gradients, vector and mixed
viscosity and cylindrical connections; wave-wave terms use the exact
complete-curl divergence identity; mean-wave terms use the actual current
means. Every nonzero term retains a physical wave factor. Distinct labels
have zero derivative products on the closed enlarged supports. Pure normal
forcing remains in

\[
 \pi_m={i\{(K^Tn-n')\cdot t_m+n\cdot f_m\}
             \over km|n|^2},                            \tag{4}
\]

so source-normal pressure is not lost by tangential projection.

The potentials and pressures are cut before the complete curl and global
zero extension. Linear regular expressions contain the same smooth cutoff;
nonlinear expressions use complete cut velocities. Therefore every mixed
entrance jet of `f_m/P` vanishes, while `P(0)>0`, and the source extends
smoothly by zero at both pulse endpoints. This discharges Q152 (10) for
every fixed output jet requested by the compact graphs, conditional on the
named current continuum estimates.

## 2. Exact tail allocation and physical high-jet bound

For a complete uncut wave `u=curl A`, pressure `pi` and cutoff `psi`, the
excluded temporal commutator remains exactly

\[
 T_\psi=(\psi_t+B\cdot\nabla\psi)u
 -2(\nabla u)\nabla\psi-(\Delta\psi)u
 +\pi\nabla\psi+\mathcal L_B(\nabla\psi\mathbin\times A). \tag{5}
\]

The per-stage split is

\[
 \mathcal R^{[j]}=\mathcal R^{[j]}_{\rm regular}
 +T^{[j]}_{\rm compact}+T^{[j]}_{\rm Gaussian}
 +T_{\rm base/exterior}.                                \tag{6}
\]

Here `Tcompact` contains the same-state `A0/A1` residuals. Spatial and slow
cutoff derivatives are inside `Rregular`, together with pressure gradients,
viscosity and all interactions. Actual compact velocities remain in both
orientations of every nonlinear product. No additive tail is forwarded to
the sourced inverse.

For a fixed physical jet `m`, requested power `q^N`,
`Lambda=3/2+2h`, and normalized compact input order `a`, the retained raw
estimate asks for

\[
 p(a)=\max\left(1,
 \left\lfloor{(N+\Lambda+2m)/h-a\over\kappa}\right\rfloor+1\right),
 \qquad J(a)=m+2p(a)+3.                                  \tag{7}
\]

Then

\[
 h(a+p(a)\kappa)-\Lambda-2m-N>0                         \tag{8}
\]

strictly. At the first sourced physical state `a=B` for the radial `A0`
input and `a=D` for the axial `A1` input. For the illustrative retained
values `h=1/200`, `m=N=0`, `kappa=10^-5`, this gives

| branch | integrations | raw source order | safe first-real ceiling |
|---|---:|---:|---:|
| radial `A0` | 30,100,002 | 60,200,007 | 60,400,023 |
| axial `A1` | 30,050,003 | 60,100,009 | 60,300,013 |

These are finite symbolic requests, not evaluated arrays or numerical
seminorms. The acyclic changed-state graph retains complete covariance,
pressure and axial velocity. The stage recurrence expands every later fixed
request back to this first-real graph in finitely many steps and archives
each required pole separately. Thus for every fixed stage `j,m,N`,

\[
 |D_{\rm phys}^m(T^{[j]}_{\rm compact}
                 +T^{[j]}_{\rm Gaussian}
                 +T_{\rm base/exterior})|
 \le C_{j,m,N}q^N .                                     \tag{9}
\]

The Gaussian square-completion estimate and same-domain polynomial
absorption provide the temporal part of (9); the compact kernel estimates
provide the `A0/A1` part. Constants can depend arbitrarily on the fixed
stage, jet, requested power and fixed profiles. No uniform-in-stage or
all-orders norm is inferred.

## 3. What happens to scale cutoffs

Scale cutoffs are not an extra summand in (6). For the fixed-gauge stage
increment

\[
 \Delta u_j=\operatorname{curl}A_j+b_je_\theta,
 \qquad \pi_j=p^{[j]}-p^{[j-1]},                       \tag{10}
\]

the complete-potential diagonal uses

\[
 V_{j,a}=\operatorname{curl}(\chi(aq)A_j)
          +\chi(aq)b_je_\theta,
 \qquad P_{j,a}=\chi(aq)\pi_j.                         \tag{11}
\]

Consequently every scale-cutoff curl, pressure gradient, time/viscous term,
both nonlinear crosses and the new square occur in the exact residual
difference. The certified diagonal chooses each cutoff only after its finite
stage constants are known. It needs no stage-uniform constant and is locally
finite on `q>0`.

Q152 supplies precisely the corrected physical seed
`B_0=D=3/2-2kappa`. Sections 1--2 supply the source/support/tail hypotheses
that were missing from its application to this family. The inherited gain

\[
 \delta_{\rm stage}=\tfrac12-3\kappa>0                 \tag{12}
\]

and finite compact recurrence now bind conditionally to the same physical
family. The certified complete-potential theorem can therefore be invoked
as a conditional implication: it constructs one locally finite real
divergence-free field whose complete residual is `O(q^N)` in every fixed
Cartesian jet, while retaining the complete pressure, viscosity, cutoffs
and interactions. This is not a new independent proof of its continuum
premises or an unconditional residual theorem.

## 4. Remaining boundary and reproduction

The next obligation is the one stated by the certified diagonal: bind this
same physical locally finite field, cutoff sequence, support and fixed
pressure gauge to the closed-time endpoint atlas and complete one-sided
jets. Localization, zero initialization, compact smooth forcing, whole-space
energy, global comparison and the core lower-bound argument must then be
checked on that same object.

No giant profile or high derivative is evaluated here. Passing tests checks
the Bell recurrence, strict physical request arithmetic, exact allocation
and frozen provenance; it is not a continuum proof. No unconditional pulse
hierarchy, admissible force, complete global residual theorem or blow-up is
claimed.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q153_physical_source_tail_diagonal_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q153_physical_source_tail_diagonal_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q153_physical_source_tail_diagonal_join/evidence.json \
  --output q153_replay.json
```
