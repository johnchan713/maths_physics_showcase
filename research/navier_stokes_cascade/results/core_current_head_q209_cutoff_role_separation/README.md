# Q209: separate background and stage cutoffs; correct Q208's identity blocker

Parent: Q208, `43813f5f730fee634fe65e3fb3b5ee73bb146edc`.
The live branch, Q208 and its active CI, Q205's complete-potential diagonal,
the canonical realization, the Borel prescription closure, Q120, Q192 and
the weighted profile/gap foundations were read first. Historical evidence,
the selected profile, pressure gauge, roots, supports, cutoff sequences and
unchanged `C^-20` activation width remain fixed.

Q208 correctly rejected an `exp(-1/x)` substitution for Q120's squared-flat
background cutoff, but then treated identification of Q207's cutoff with
Q120's function as the next blocker. That conflates two cutoff roles. This
checkpoint corrects the frontier without modifying Q208 or selecting a new
cutoff.

## 1. Two different locally finite constructions

The coefficient-background Borel sum uses coefficient index `n`:

\[
 q^{2nh}\chi_{bg}(a_nq)(E_n,U_n,W_n^r,J_n^r).       \tag{1}
\]

The executable background audit imports the squared-flat step

\[
 g(x)=e^{-1/x^2},\qquad
 \sigma(x)={g(x)\over g(x)+g(1-x)},
 \qquad \chi_{bg}(s)=1-\sigma(2s-1),                \tag{2}
\]

with flat extensions. Q120 proves its explicit Euler bound
`K_m<=82*4097^m(m!)^3`. The later Borel-prescription closure also proves
that the background diagonal works for any one fixed smooth cutoff using
that cutoff's own finite seminorms; equality to the diagnostic squared-flat
function is not needed for the rule-defined existence theorem.

The generated-cycle diagonal instead uses stage index `j`:

\[
 V_{j,a}=\operatorname{curl}(\chi_{st}(aq)A_j)
          +\chi_{st}(aq)b_je_\theta,
 \qquad P_{j,a}=\chi_{st}(aq)\pi_j.                 \tag{3}
\]

The Q205 provider says only: fix one smooth `chi_st`, equal to one for
`s<=1/2` and zero for `s>=1`. The canonical source keeps `chi(a_j*q)`
symbolic. It does not name (2). An older manufactured cutoff diagnostic uses
`exp(-1/(s-1/2))` and `exp(-1/(1-s))`, a different flatness profile, and its
own record explicitly says the polynomial fixture is not the actual
supported giant field.

Using the same letter `chi` does not identify (1) and (3). Their indices,
scale rules and analytical roles differ. Therefore Q120's numerical bit
budget cannot be imported into Q207's stage diagonal, and Q205 does not
require such an import.

## 2. The stage cutoff constant exists without a formula

Let `E=s partial_s`. For every integer `m>=0`,

\[
 E^m\chi_{st}(s)
 =\sum_{k=0}^m S_2(m,k)s^k\chi_{st}^{(k)}(s),       \tag{4}
\]

where `S_2` is the Stirling number of the second kind. Every positive
ordinary derivative is supported in the compact transition `[1/2,1]`.
Since the chosen cutoff is smooth, each term in (4) has a finite supremum.
For every fixed finite request `R`, therefore,

\[
 K_R^{st}:=\max_{0\le m\le R}\|E^m\chi_{st}\|_\infty<\infty. \tag{5}
\]

This is a continuum consequence of smoothness and compactness, not a sampled
or executable estimate. It gives no numerical value and no uniform control as
`R` tends to infinity. At Q207's finite stage-602 ceiling one may use (5)
with `R` (or conservatively `R+1`) without evaluating its
`36,479,350,301` entries.

## 3. Why an abstract finite constant suffices

After differentiating (3), Q205 has for each fixed stage and jet a finite
constant `C` and logarithmic power `P`. It chooses `s_j=log a_j` only after
the finite list of current constants is known, above

\[
 \frac{2(\log C+j)}{\gamma_j}+\frac{P^2}{\gamma_j^2},
 \qquad \gamma_j={h\delta j\over2}>0.               \tag{6}
\]

The cutoff seminorm (5) is simply part of `C`. The elementary
`log(1+s)<=sqrt(s)` estimate then gives on the whole interval
`0<q<=a_j^-1`

\[
 C(1+|\log q|)^Pq^{\gamma_j}\le2^{-j}.              \tag{7}
\]

Thus arbitrary growth of `K_R^st` at finite requests is absorbed stage by
stage. No stage-uniform cutoff constant, Q120 bit budget, uncut series
convergence or numerical scale list is used.

This proves only the cutoff-shape portion of the conditional diagonal. It
does not evaluate the stage representatives or promote their continuum
providers.

## 4. Corrected remaining frontier

Q208's proposed identity
`chi_Q207=chi_Q120` is not the next mathematical obligation. The first
unresolved direct input remains a selected-profile and first-coefficient
weighted majorant on the actual common domain:

\[
 \max_{|I|\le R+1}\sup_{\mathcal D}
 w_I S^{-d_I}
 \left|D^I(F_0,U_0,\Pi_0,M_0,V_0,
                 E_1,U_1,W_1,T_1/\zeta)\right|
 \le A_R,                                           \tag{8}
\]

with `R=36,479,350,300` and `S=ell^2`. Complete pressure memory,
`Pi1/Pi2`, the whole moving pressure product, axial velocity and every
selected root/frame coefficient remain in the defining maps.

One must also prove on the same object a common positive `g_R` for all
selected Jacobian, normal, root, determinant, variance and support gaps, and
finite compact/slow inverse constants at every archived pressure and moment
pole. Q192 propagates a supplied primitive array and floors through finite
recurrences; hashes and tests do not prove (8) or those continuum minima.

Q209 therefore establishes a justified correction and a finite stage-cutoff
estimate, but no numerical high-jet cutoff norm, unconditional pulse
hierarchy, complete residual, admissible force, competitor exclusion or
Navier--Stokes blow-up. The next obligation is (8), the common gaps and the
inverse archive from the explicit selected maps.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q209_cutoff_role_separation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q209_cutoff_role_separation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q209_cutoff_role_separation/evidence.json \
  --output fresh.json
```
