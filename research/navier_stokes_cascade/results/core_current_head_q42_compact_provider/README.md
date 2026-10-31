# Q42: explicit compact provider for the current fast inverse

Parent: `39897b3143afbd749a344be22ab925505ef624e2` (Q41).
On 31 October 2026 parent focused CI 41152311182 passed; broad CI
41152311236 was still running. Historical frozen evidence is unchanged.

Q40 accepted independent KT and KA constants. This checkpoint supplies
that interface using Q29's weighted shifted-integral theorem. It removes
those two independent operator inputs, conditionally on the original raw
multiplier jets, plateau gap, geometry, kernel and input-weight hypotheses.
It does not certify those lower hypotheses or the initial A1 tail.

For raw coefficient order j and EXACT requested pole v, use the Q29
notation f,g,k,chi,1-chi,chi_U,J,L,E,C,Z. The provider computes

    KT = 2^(3j) g (chi+(1-chi)) f J C,
    KA = 2^(3j) k chi_U f J L E Z.

The derivation is the fixed-characteristic-shift Leibniz estimate: four
factors cost 2^(3j); differentiation at fixed shift introduces no M.
After differentiation, the native-to-logarithmic integration Jacobian is
bounded by J. Directed integrals use C only on the branch admitted by the
cutoff; the full shifted moment uses L E, and the output transition uses
Z. The latter is the inverse-weight bound on the ACTUAL cutoff transition,
not a rectangular union of moving supports. There is no extra L in KT.
The complete KA term survives even for zero-Haar input. Q40 still uses
1+KA in its axial velocity and supplies the physical S factor exactly once.

The original time-inverse/S bound is retained as a separate lower input.
All constants are finite nonnegative polynomials with symbolic exponential
coefficients in the existing exact ring. Q29's rational diagnostic formula
is used as an independent test oracle, not as a restriction to rational
profile constants. This provider can be passed directly to Q40.FastFields
or Q41.build, using their same imported module family.

## Coverage and the cutoff derivative

Q29's chi_U norm at order j uses the cutoff through j+1. Q40 requires
KA at fast output m and KT at m+1. This provider returns both maps at each
requested order, so its sufficient cutoff ceiling is m+2. Q41 reserves
fast fields through m+5; its sufficient compact/cutoff ceilings are thus
m+6 and m+7. The extra cutoff derivative is mandatory for the combined
provider, though a specialized KT-only call could use fewer cutoff jets.
These counts do not replace the much larger integration-by-parts tail DAG.
Each request uses its exact pole. A declaration covering a larger pole
cannot silently replace the source-weight hypothesis at the requested pole;
constants need not be monotone in pole. Derivative and pole guards reject
insufficient coverage. They do not verify a submitted numerical bound.

## Precise unresolved inequalities

For the original selected profiles and every required band, the remaining
input statement is Q29's raw weighted hypothesis

    sup_{|I|<=j} |D_(U,Z,T,y)^I F| <= A_j zeta delta^(-v),

together with each declared raw multiplier seminorm bound and one uniform
actual plateau gap eta>0. Native/R normalized derivatives must be converted
before applying this hypothesis. The provider does not manufacture that
conversion, change the cutoff, or deduce a uniform gap from pointwise
compactness. Exact plateaux, flat zero extension and unchanged shift are
also required. If these hypotheses hold, Q29 proves the two operator
bounds above; passing the finite tests does not prove the hypotheses.

Separately, Q40 still needs the SAME initial axial A1 remainder bounded at
least at order D=2-4kappa, at the exact requested pole and through the whole
row ledger. Neither KT nor KA proves an epsilon-small integration-by-parts
tail: the inequality with M^(-p) needs the source's higher jets. No initial
A1 remainder is set to zero. Complete pressure memory, moving multipliers,
old means and all six covariance entries are untouched by this connection.

Tests compare exact constants with Q29, retain distinct requested poles,
reject missing cutoff jets and check the Q40 fast formula with polynomial
lower data. Frozen Q41 evidence and its transitive bindings are replayed.
Next connect the actual multiplier/chart contracts and initial A1 tail to
the existing all-order source graph, then assemble the complete final rows.
The zero strict fast-residual margin, uniform pulse estimates, full source,
viscosity/cutoff interactions and blow-up remain unresolved.
