# Q83: carry the radial principal moment into the axial pressure row

Parent: `53f1b98642391752ca4a2e4477a85f8eeced8144` (Q82).
There is one common integral in Q56 and Q82, before their different output
multipliers. Estimating it on the radial scale puts the entire moving-cutoff
part of the axial principal moment one epsilon order above the target. We
also give the leading-profile requests for the two remaining estimates.
This is conditional on the original common-domain comparison, spectral,
chart, smooth support, inverse and lower-field contracts. The historical
records and the physical fields are unchanged.

## Identify the common integral before projection onto the output collar

Work labelwise in the SAME frozen common chart of Q56/Q57/Q82. Write
brackets for the whole-line shifted dU integral at fixed output (Uout,y),
Q for removal of the auxiliary Haar mode, I for the zero-mode-free temporal
inverse including its c_i^-1, and P for the original dimensionless principal
tensor in common units. Use U=R^d, D=R_U/R, F=V/R and C=c_chi=s mu. Slow
weights commute with Q and I but generally do NOT commute with radial
integration. C is independent of the INPUT radius; it is derived from the
complete pressure cutoff transition. Define

    J = [ D Q(Pθθ−Prr) + (2 F_U−4 D F) I Q Prθ ].                 (1)

Indeed f0=R_U and a1=6V/R²−2V_R/R give −f0 a1=2F_U−4DF.
Thus (epsilon/2) B0 J is the principal part of POSITIVE A0 Q gr in Q82;
the final radial row takes its negative. Remainders of Q82 are not in J.
Now define

    K = [ D R² Q(Pzz−(Prr+Pθθ)/2)
          + 2R G_U I Q Prz − R² F_U I Q Prθ ].                    (2)

Then Q57's full axial principal term is exactly

    (epsilon²/2) B1 {partial_Z K + C partial_Z J + C_Z J}.        (3)

To check (3), the diagonal weights in Q58 are DR² and D(C−R²/2);
its last shear weight is C(2F_U−4DF)−R²F_U. Expand in the same tensor
entries, with the same temporal inverse and translations. The exterior Z
derivative is taken at fixed shift and frozen label; B1 is outside it.
Both diagonal combinations vanish for an isotropic addition. The two
remaining shear terms in K must still be estimated. Equations (1)--(3)
also hold with the fixed-carrier comparison tensor and leading F,G if all
terms are changed consistently. They are not identities between pressures
of different stages: only the extracted ORIGINAL principal entries recur.
The post-signed axial and final radial remainders retain their own states.

## What the radial estimate actually implies for the axial row

Let t=3−5κ. Impose Q82's successful-branch principal shell hypotheses at
output n+1, including the conversions to unweighted native norms. The
Fourier proof bounds the integral before multiplying by B0: equivalently
use output multiplier one in that proof. It gives

    N_(n+1)((epsilon/2) J) <= epsilon^t P_(n+1)(S).                (4)

Here and below P_j(S) denotes a positive bounding polynomial, not the
tensor P. Q82's direct comparison error is order three at the requested
jets, so it also fits (4). On its other branch the direct integral is O₁
with the epsilon prefactor; inequality (2) of Q82 pays t−1<2 from the
same frozen polynomial gate. These uses require the polynomial input
bounds, not merely a scalar output value. No derivative hits a branch
selector. On the successful branch both proofs use the original inverse
and the same target. Q81 localizes their differences, not J or K.

Multiplication by C and differentiation once in Z now give

    N_n((epsilon²/2) B1 partial_Z(C J))
      <= epsilon^(t+1) 2^(2n+1) N_n(B1) N_(n+1)(C) P_(n+1)(S),  (5)

where the through-order norms are raw derivative maxima and epsilon is
frozen. The coefficient follows by two Leibniz bounds. The required
weighted output conversion is on B1's own transition as in Q57; no
unweighted-to-weighted conversion at a zero of B0 is used. In particular
(5) retains C_Z=mu s_Z. The estimate is for the shared principal term
only. Nothing in (5) makes the complete axial error order 4−5κ.

Merely knowing N_j(B0 J) does not give (4): a smooth J supported where
B0 vanishes can have arbitrarily large norm without changing the product.
This is an obstruction to that *inference*, not a constructed physical
pulse. We use the pre-multiplier shell proof and its whole-hull contract.
To reach t it now suffices to estimate the three terms in K via Q57 plus
the original axial curl and mean remainders. For this implication the very
large radial estimates need one extra output jet; it is optional to obtain
(5). Q57's separate four-weight calculation remains valid without it.

## Same-band leading inputs for this simultaneous estimate

For a shell with count p, write m for the output order and z=1 if the
integrand has the exterior Z slot, zero otherwise. Q59's squared-column
solve and Q60's affine-path recurrence use at most m+z+p+1 slow slots,
m+z+2p−r total slots and m+z axial/time slots (r=0,1 denotes I^r).
At fixed path time the generator can add one slow field derivative; positive
path orders add none. Q62's unstripped target can add one derivative of a
primitive, including the full terminal angular pressure trace. Q63's
subindex argument applies unchanged: each eta derivative in the chart
requires an axial/time slot. For the leading fields a sufficient request is

    i+j <= m+z+p+2,    j <= m+z+1,                              (6)

where i,j are native X,eta orders. The same count bounds the weights in
(1),(2): the geometrical powers cost fixed positive-hull constants, while
F_U and G_U take just one extra radial slot. Chart composition uses the
whole interval; zero-velocity support cannot shorten the leading integral.
For K take (m,z,p)=(n,1,99995), giving n+99998, j<=n+2. For (4) take
(n+1,0,199995), giving n+199998, j<=n+2. Thus **one** downward-closed
leading profile table with the latter bounds covers both principal channels
including the extra output jet in (5). It does not request hundreds of
thousands of angular derivatives of the implicit compact repair at n=0.
The terminal moment traces in Q62 need at most n+2 angular slots.

More explicitly: take Q63--Q78's finite mixed bounds on the shell with
these indices, and the original core/heat bounds wherever the positive
hull extends outside it. Q62 supplies the target, Q58's finite chart
composition supplies the weights, and Q59--Q60 give the fixed-index tensor
polynomial from the slow column inverse, ambient value bound, cutoff/atlas
and partition polynomials. All these recurrences are finite at fixed n;
their degrees may depend on n but not on the running band. Insert the
result into the two shell bounds on Q61's common comparison branch. This
removes an independent high principal-product norm **relative to** those
same-field provider contracts. It does not prove those primitive contracts
from a grid. The other branch uses the true low-order bounds directly as
in Q61/Q82. No extra q threshold is chosen for (5) or (6).

Next quantify Q82's nonlinear Rem through n+199996 and its lower-field
requests on the same domain, including the internal temporal inverses,
final angular equation, complete curl and mean updates. Neither J nor K
is the full pressure; both exceptional reserved orders and the moving
compensator remain. Unconditional actual-pulse estimates, complete
residual improvement and blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. The symbolic tests check the independent
original weights and moving derivative; they do not certify profile norms.
