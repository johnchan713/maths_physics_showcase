# Q192: auxiliary frame jets and corrected physical graph losses

Parent: Q191, `9bad836772511537fc9ee491f5d12406a2528bed`.
The live branch, Q189, Q191, the leading-frame and physical-pulse providers,
their frozen evidence, and the active Q189--Q191 CI runs were read first. The
selected Q189 profile, full pressure memory, supports, carrier choices and
unchanged `C^-20` activation width remain fixed.

Q191 left the all-fixed frame/graph jet inequality as the next obligation.
This checkpoint closes the rational auxiliary-frame transfer relative to the
named primitive profile-jet estimates and corrects an ambiguity in Q191:
physical radial graph derivatives have an `epsilon^-kappa` loss per derivative
before the physical `Q^-1/2` scaling. A pure polynomial-in-`S` estimate applies
to auxiliary frame jets, not to the already pulled-back physical graph jets.

## 1. Reduce the primitive comparison error

The actual phase/frame comparison has the form

\[
 {1\over L_s}+(1+L_s)
 \left(S^{-3}+\epsilon+\epsilon^2+k^{-1}\right).     \tag{1}
\]

The accepted later-band gate and pulse geometry are

\[
 S^2(\epsilon+\epsilon^2+k^{-1})\le1,
 \qquad 2r_0S\le L_s\le S,
 \qquad r_0={1\over6057072}.                         \tag{2}
\]

For `S>=1`, (2) gives

\[
 S^{-3}+\epsilon+\epsilon^2+k^{-1}\le {2\over S^2},
\]

and therefore

\[
 {1\over L_s}+(1+L_s)(\cdots)
 \le {1\over L_s}+{4\over S}
 \le {5\over L_s}.                                 \tag{3}
\]

After `j` auxiliary derivatives, the frozen representative is not
differentiated and the compact profile factors have a finite `C_j S^j`
bound. Thus the primitive actual/reference differences needed below are
`C_j S^j/L_s`. This is a conditional analytical use of the named profile-jet
provider, not a numerical reconstruction of a historical first band.

## 2. One denominator floor controls every fixed order

If `g>=gamma>0`, differentiating `g g^{-1}=1` gives the exact majorant
recurrence

\[
 H_0=\gamma^{-1},\qquad
 H_j\le {1\over\gamma}
 \sum_{i=1}^j {j\choose i}G_iH_{j-i}.               \tag{4}
\]

The product recurrence is

\[
 (FG)_j\le\sum_{i=0}^j{j\choose i}F_iG_{j-i}.       \tag{5}
\]

Equations (4)--(5), applied to the difference identity
`g^{-1}-g0^{-1}=-(g-g0)/(g g0)`, show that reciprocal differences retain the
same `1/L_s` factor. They require only the value floor `gamma`; higher orders
increase finite constants and powers of `S` but impose no new smallness gate.

The actual frame has the fixed floors

\[
 |n_{tan}|\ge B_*/2,\qquad |n|^2\ge B_*^2/4,
 \qquad |c_0|\ge c_*>0.                             \tag{6}
\]

Consequently normalization of `n_tan`, `s_a=n_r/|n_tan|`, the explicit
`Bleft`, `|n|^-2`, square roots and all products are closed by (4)--(6).
Because `A_n` contains `n'` and the frame error contains `B'`, order `j`
uses primitive jets through order `j+1`. A conservative resulting statement
is

\[
 |\partial^\alpha E_{frame}|
 +|\partial^\alpha(d-d_{ref})|
 \le {C_\alpha S^{2|\alpha|+4}\over L_s}.           \tag{7}
\]

The exponent in (7) is intentionally nonsharp. Its finiteness and the common
value domain matter; no order-uniform constant is asserted.

## 3. Full pressure denominator and numerator

The same reciprocal recurrence applies to the complete constrained pressure

\[
 \pi={i\over km|n|^2}
 \left(n\!\cdot K_{bg}t-n'\!\cdot t+n\!\cdot f\right). \tag{8}
\]

All three numerator channels in (8) remain present. In particular,
transversality does not erase either the moving-normal term or the
source-normal term. The specialization `f=0` applies only to the homogeneous
primary pulse and cannot be propagated to a sourced inverse.

## 4. Correction for physical graph derivatives

The exact coefficient operators include

\[
 D_R=\partial_R+M_i d_rR^{d_r-1}L_i,
 \qquad D_Z=\epsilon\partial_Z,
 \qquad D_T=-\epsilon\partial_T+c_iN_i,              \tag{9}
\]

with `M_i<=epsilon^-kappa S^-rho_g`. Repeated radial graph derivatives also
differentiate `R^(d_r-1)`, but the pulse stays on a compact annulus away from
`R=0`. Thus `a` spatial graph derivatives can cost
`epsilon^(-a kappa)` times a finite polynomial in `S`. After physical scaling,
the worst coefficient losses are

\[
 a\left({1\over2}+h\kappa\right)+b(1+h).            \tag{10}
\]

The harmonic phase budget is

\[
 a\left({1\over2}+{h\over2}\right)
 +b\left(1+{3h\over2}\right).                      \tag{11}
\]

Subtracting (10) from (11) leaves the exact margin

\[
 ah\left({1\over2}-\kappa\right)+{bh\over2}>0
 \quad\text{whenever }a+b>0,                       \tag{12}
\]

because `kappa=1/100000<1/2`. Hence the omitted graph loss in the literal
reading of Q191's displayed inequality is real but is strictly absorbed by
the already selected phase exponents. Q191's inequality should be read as an
auxiliary-frame statement; (10)--(12) are the physical pullback statement.

## 5. Remaining boundary

Q192 has not independently completed the field estimates. The next obligation
is to insert (7) and the corrected losses (10) into the differentiated Duhamel
and Faà di Bruno recurrences for the actual uncut amplitude, then assemble the
complete curl, pressure and potential jets while retaining every cutoff,
cylindrical connection, graph coefficient and support-endpoint term.

The named physical-pulse provider supplies that conclusion relative to its
inputs. This checkpoint does not claim an unconditional pulse hierarchy,
complete residual improvement, admissible force, comparison result or
Navier--Stokes blow-up. Executable recurrences and frozen hashes check the
ledger; passing tests are not a continuum proof.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q192_frame_graph_jet_transfer/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q192_frame_graph_jet_transfer/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q192_frame_graph_jet_transfer/evidence.json \
  --output q192-fresh.json
```
