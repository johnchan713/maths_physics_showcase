# Q146: specialize the complex caps and expose the exact-repair identity gap

Parent: Q145, `db174e5ac98dad3fbad4718a45c5442e3bbfb4f0`.
The live branch, Q145 evidence and current CI, Q72/Q73, Q127, Q131/Q132,
Q137/Q138 and the inherited regional cone join were read first. Historical
evidence, roots, supports, pressure datum, frequency, scales and the
`C^-20` width remain unchanged.

Q145 reduced the scalar hierarchy gate and left one same-profile complex-cap
statement open. The first result here is that the statement is already a
consequence of the Q73/Q127 formulas once their provenance is kept. The
second result is a correction to the proposed next step: Q138's stage-only
repair join does not certify that the supplied Q137 repair corrects the same
modulation memory. The full five-row target binding in Q132 remains open.

## 1. The actual Q73/Q127 cap tuple

On Q73's strip let `e,v` be its actual complex caps for `E,E^-1` at the
shape entry. Q72 derives them from the same core trace and all five incoming
moments before the prescribed continuation. Q73 proves

    L0=J6+e v,
    |ell|<=L0,
    |E|<=exp(T/10+L0+4)/C,       T=2000 B0.             (1)

The same displayed logarithm also gives

    |E^-1|<=C exp(T/10+L0+4).                            (2)

No lower bound on a complex real part is being inferred here; (2) follows
from the explicit nonvanishing exponential formula. Q127 proves for the
same earlier accepted root

    |E|<=4P,  |E^-1|<=20/P                              (3)

on the earlier annular patch and ideal exit. Q143 supplies the incoming
endpoint caps. Therefore one valid tuple is the maximum of the endpoint
caps and (1)--(3). Since `exp(x)<3^ceil(x)` for `x>=0`, it has the exact
outward representation

    n=ceil(200 B0+L0+4),
    Ecap=max(Eendpoint,3^n/C,4P),
    invEcap=max(invEendpoint,C 3^n,20/P).                (4)

The implementation evaluates (4) only for moderate controls. At the selected
hierarchy it records the expression symbolically; forming the integer
`3^n` is neither needed for the analytical implication nor feasible. Q73's
strip, Q143's endpoint and Q127's root all refer to the unchanged selected
profile, so (4) discharges Q145's equation (5). The real shear floor still
uses the sharper real root cap, not Q127's complex increment radius.

## 2. Stage compatibility is not repair identity

The mixed route is

    Q140 four fields -> Q139 primitives -> Q138 modulation,
    Q132 target/root -> Q137 physical repair -> Q138 join_repair.

The last helper checks only the two stage strings. It does not check that
Q132's target was formed from the same modulation. Q132 itself correctly
requires, for every row `r` in `{M,I,J,S,Pi}`, actual constants `T_r,L_r`
such that

    |partial_eta^m d_r^base| <= T_r L_r^m (m!)^3,        (5)

then adds the modulation-memory prefactor for that same row. Its normalized
target retains

    (m,(m-j)/lambda,i,-s,cp),                            (6)

including the inverse-`lambda` loss and the pressure row. The accepted
pointwise root and positive-field neighborhood must also be the original
ones. A zero baseline entry is allowed only when an actual preceding target
identity proves it; ending the modulation support is insufficient because
the cumulative moments remain.

No current machine record binds Q138's mixed modulation object to Q131's
five discrepancy functions and the five actual baseline pairs in (5).
Consequently adding Q137 and Q138 envelopes is a conditional norm statement,
not an exact-repair or cone certificate. The implementation now emits a
five-row repair request and deliberately leaves `same_modulation_memory`,
`baseline_defects`, `accepted_root`, `positive_field` and `cone_join` false.
It rejects a missing `Pi` row or changed pressure-memory identity.

## 3. Scope and next obligation

Exact controls check (4), including endpoint domination and the annular
caps, then feed the specialized tuple through Q144 and Q145. A separate
control verifies that all five Q132 rows and the inverse-`lambda` combination
are requested, and that a stage-compatible Q137 object alone cannot mark the
join complete. The parent closure and source hashes replay byte for byte.

These are analytical implications and manufactured exact controls, not an
evaluation of the selected enormous profile or a continuum proof. Next bind
Q131's five full-memory discrepancies to the same Q138 mixed modulation,
supply (5) from the actual pre-repair target identities, and recheck the
original root/positive-field neighborhood. Only then may Q137 be joined to
Q138 and compared with the already protected inner and outer value cones.
Coefficient and reserved-mean supports, uniform physical pulses, complete
wave-covariance axial/pressure means, the full viscous/cutoff/interaction
residual and blow-up remain unproved.
