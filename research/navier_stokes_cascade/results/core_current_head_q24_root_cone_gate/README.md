# Q24: exact finite-head column and positive-root gate

Parent: `47934df884d9a195ac23906485ddec98cd651726` (Q23).
The upper propagator completion does not supply the lower constants K0,beta
in `core_primary_root_budget` (1). This checkpoint reduces that remaining
finite-head obligation to three oriented determinant inequalities and gives
rigorous perturbation bounds sufficient to certify them. No actual selected
column enclosure is asserted. The precise missing inequalities replace an
appeal to compactness or positive pulse mass; the field and selector stay fixed.

## 1. Use the canonical averaged columns, not a pointwise polarization

Use the original unrooted principal flux matrix H from
`core_canonical_haar_column_jet_closure` (1), with BOTH its angular factor
1/2 and its transverse Haar measure. Write H=s0 C, s0=r0^2/sqrt(L)>0,
C=[c_plus,c_minus]. The columns are the actual canonical pulse integrals,
including the fixed temporal cutoff. They are not the complete curl-wave
covariance and cannot be replaced by a leading polarization at one point.

Let rho=|T|/zeta and tau=T/|T| on the open shell. Q2's factored stress
bound gives rho>=rho_min>0 and a smooth limiting direction tau at both
flat edges. The magnitude rho need NOT be bounded: at the outer edge it
can grow as distance^-3. Thus compactness must be applied to C and tau,
not to the generally unbounded T/zeta. The unrooted columns contain no
primary square root, so testing this gate is acyclic.

For a fixed ordering of the two signs define

    D=det(c_plus,c_minus),
    Aplus=det(tau,c_minus), Aminus=det(c_plus,tau).       (1)

If D has negative orientation, multiply ALL THREE determinants by the
same sign o=-1; otherwise use o=1. Do not take their absolute values
independently. On any connected domain where D never vanishes this
orientation is fixed; it can be selected separately for discrete labels.

Cramer's rule gives the EXACT normalized squared roots

    b_plus=rho Aplus/D, b_minus=rho Aminus/D,
    y_sigma=(zeta/s0)b_sigma.                          (2)

Therefore invertibility alone is insufficient: positive roots require
Aplus/D>0 and Aminus/D>0. At a zero numerator the original-root signed
map divides by zero even if C is perfectly invertible.

## 2. Three lower margins and explicit root constants

On the finite earlier-band set singled out by Q23, and on every closed
slow chart including shell-edge limits, a sufficient uniform certificate is

    oD>=d>0, oAplus>=a>0, oAminus>=a>0,
    |c_plus|,|c_minus|<=M, rho>=rho_min>0.               (3)

The Euclidean column upper bound M is available conditionally from Q23
and the canonical Haar integration recipe. The three LOWER bounds are
not consequences of it. Since |D|<=M^2, (2) gives

    b_sigma>=beta_head=rho_min*a/M^2,
    ||C^-1||_infinity<=K_head=2M/d.                    (4)

For the inverse estimate each row of adj(C) has absolute sum at most
sqrt(2)M<=2M. No inverse Gaussian minimum, radius replacement or extra
q exponent is introduced. Combine the head with an already certified
canonical direct-cone tail by taking beta=min(beta_head,beta_tail) and
K0=max(K_head,K_tail). If the head is empty, only the tail is needed.
These constants feed the EXISTING root jet recurrences, unchanged.

Pointwise strict inequalities on a finite union of compact closed domains
would give (3), but establishing them is the obligation. Finitely many
bands do not mean finitely many slow points. Sampling, exact positivity
at the frozen centers, or invertibility of the ODE fundamental matrix
cannot replace a bound over the whole charts. A nonsingular linear flow
can have singular selected covariance columns after integration.

## 3. An explicit enclosure test for actual column calculations

Suppose a closed subdomain has reference vectors c_plus0,c_minus0,tau0
and proved Euclidean error bounds e_plus,e_minus,e_tau for the actual
vectors. The reference tau0 need not be exactly unit. Put bounds
U_plus>=|c_plus0|, U_minus>=|c_minus0|, U_tau>=|tau0|.
Bilinearity of determinant and |det(x,y)|<=|x||y| give

    ED=e_plus U_minus+e_minus U_plus+e_plus e_minus,
    EAplus=e_tau U_minus+e_minus U_tau+e_tau e_minus,
    EAminus=e_plus U_tau+e_tau U_plus+e_plus e_tau.      (5)

Consequently the exact sufficient local tests are

    o det(c_plus0,c_minus0)-ED>0,
    o det(tau0,c_minus0)-EAplus>0,
    o det(c_plus0,tau0)-EAminus>0.                     (6)

Take d and a to be the minima of these positive remainders over a FINITE
cover of the head, and M the maximum of U_plus+e_plus,U_minus+e_minus.
Then (3)--(4) follow. Formula (5) retains the quadratic error products;
dropping them can falsely certify a boundary case. Use directed bounds
for computed integrals and analytic chart variation, not point samples.
The executable test uses rational l1 upper bounds for U to avoid square
roots; this is conservative for Euclidean error balls. It reports each
margin and refuses to return root constants if a margin is nonpositive.
A failed sufficient test is not proof that the actual canonical roots fail.

The actual values/errors in (6), or a direct proof of (3), are NOT supplied
by the frozen records inspected here. `canonical_value_domain` and
`core_primary_tail_recheck` supply conditional later direct-cone gates,
not certified arbitrary earlier-head columns. If the original rule-defined
accepted domain is proved already inside that actual direct-cone tail,
the head can instead be empty. That requires the same constants/domain
identification; Q23's enlarged analysis threshold cannot establish it.
Historical-versus-canonical ROOT COMPARISON is not needed for (2).

## 4. Exact controls explain the obstruction

Take constant pulse vectors t_plus=t_minus=(1,1,0) on a unit temporal
interval with unit cutoffs. Each radial component and each pulse's kinetic
mass is positive. Each principal radial/tangential flux column is (1,0)
up to the same positive Haar prefactor, so D=0. This is a manufactured
algebraic example, not the actual normal, entrance or Navier--Stokes field.
It proves that positivity of radial amplitude and Gaussian mass alone
cannot imply independence of the two sign columns.

Separately, C=I and tau=(3/5,-4/5) have D=1 but Aminus<0;
tau=(1,0) gives Aminus=0. Both controls have perfect inverse bounds.
They detect a false promotion from an upper/absolute determinant estimate
to a positive canonical root. Simultaneously reversing column orientation
changes all determinant signs and preserves the two Cramer ratios.

The full covariance has six entries and includes curl correlations and
squares. None is removed by solving the two principal flux rows (2).
Full pressure, Pi1/Pi2, current axial velocity, A0/A1 and the complete
moving c_rho P moment remain in their separate mean/source ledgers.

## 5. Result, next obligation and reproduction

Established analytically: (1)--(6), including the flat-edge normalization,
orientation, finite-head aggregation and explicit K0,beta consequences.
Established computationally: exact manufactured positive/negative controls,
error-budget checks and frozen provenance replay. Not established: the
actual three determinant margins on the original head, or identification
of that head as empty. These are precise blocking inequalities for the
canonical-root input of the first signed-prefix source bound.

Next derive whole-chart enclosures for the canonical column integrals and
target directions in (6), or prove original-domain membership in the actual
direct-cone tail. No unconditional source/pulse bound, complete residual
improvement, endpoint result or blow-up is claimed.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here.
