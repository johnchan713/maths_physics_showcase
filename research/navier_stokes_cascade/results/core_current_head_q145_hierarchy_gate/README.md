# Q145: the hierarchy discharges the scalar gate, but not the profile binding

Parent: Q144, `a28bbd8bcc13f8cee03f7e1a2ac5e8c00d2a5506`.
The live branch, Q144 frozen evidence and current CI, Q64/Q65, the original
loop-modulation bounds, the matching-root frequency selection, Q131 and Q140
were read first. Historical evidence, fields, roots, supports, pressure datum,
frequency, scales and the width `t=C^-20` are unchanged.

Q144 derived a positive real shear floor and one fixed mixed recurrence on the
shape, axial restoration and earlier annular join, starting from all five
actual incoming moments and retaining both `Pi` and `Pi_eta` in `p2`. It then
listed

    G>=max(2,|mu|,|p2|,|mu p2|),
    Eexp>=exp(G),  V0<=N                              (1)

as regional estimates still to be proved. This checkpoint audits (1) against
the unchanged original hierarchy. The audit shows that (1) is not an
independent regional obligation. It also identifies the smaller obligation
that prevents the conditional Q140--Q138 chain from becoming an actual cone
join.

## 1. Exact hierarchy comparison

The original construction fixes

    A=C^4096 >= 2^16,
    G=exp(A^32),
    H=exp(G^4),
    N=1+floor(H^32).                                  (2)

The compact-envelope contract used by the original loop proof gives
`|p2|<=A` on the whole compact profile, with `p2` containing the incoming
pressure datum and all five moment rows. The signed-root calculation in that
same proof gives `|mu|<G` and `|mu p2|<G`, including zero pressure and the
zero-variance boundary. These are value bounds; no high derivative norm is
inserted into an exponential.

All remaining comparisons are immediate and symbolic. Since `A>=2`,

    log G=A^32>log A,  hence A<G,
    G^4>G,             hence exp(G)<H,
    H^32>=H,           hence N>H^32>=H.                (3)

Consequently the unchanged choices `Eexp=H` and `V0=H` discharge (1),
conditional on the original whole-profile value contract. No value of
`G,H,N` is rounded or materialized. This also explains why replacing `Eexp`
by a newly evaluated floating-point exponential would weaken, rather than
strengthen, the proof ledger.

The pressure conclusion uses Q144's full-pressure identity. A reset pressure
or a recurrence omitting either `Pi` term would not be covered by the
whole-profile contract.

## 2. What the audit does not bind

Q144's function still accepts the free tuple

    (Hstrip,Z,B0,J6,L0,Ecap,invEcap,R,P,length),        (4)

where `Hstrip` is an angular strip radius and is unrelated to the hierarchy
quantity `H` in (2). Its guards prove that a supplied tuple closes the mixed
recurrence. They do not prove that one tuple bounds the selected Q73/Q127
profile at every real center on this region. In particular the missing
same-profile statement is

    sup_(shape/restoration/earlier-annulus, |z-eta|<rho)
      max{|ell|/L0, |E|/Ecap, |E^-1|/invEcap} <= 1,      (5)

with the Q73 strip, Q143 endpoint and Q127 root referring to the same accepted
branch, and with the resulting four-field mixed envelope for `(a,b,p2,E)`
using one order-independent tuple. The `U` cap is then supplied by Q144's
explicit restoration recurrence rather than an additional free input. The
implementation records (5) as false; it does not turn a caller-supplied
Boolean into a certificate.

Thus the precise current obstruction is (5), followed by specialization of
Q144's free caps to that actual profile. Once it is proved, the scalar part of
Q140's gate may use (2)--(3), Q140 can feed Q139/Q138, and the same full
five-row repair can be compared with the inner and outer cone estimates.
The present checkpoint does not perform that join.

## 3. Checks and scope

Exact integer-exponent checks verify every implication in (3), the unchanged
frequency formula and the distinction between `Hstrip` and the hierarchy
`H`. A replay constructs Q144's recurrence, verifies its positive floor and
complete pressure memory, reduces its four scalar requests to the single
profile-binding obstruction, and rejects a changed identity or a pressure
reset. The frozen closure and source hashes are replayed byte for byte.

These checks audit implications between established contracts. They do not
independently re-prove the compact envelope, evaluate the selected enormous
profile, establish (5), or constitute a continuum proof. The actual compact
cone join, coefficient/reserved-mean supports, uniform physical pulse bounds,
complete wave-covariance axial/pressure means, full viscous/cutoff/interaction
residual and blow-up remain unproved.

