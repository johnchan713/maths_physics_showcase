# Q103: current Haar history and the sourced postslow scalar reset

Parent: `7a8aee9a6c509d0e2950b7135a8fd745efbb3e47` (Q102).
Its CI runs 41465340032 and 41465339945 were still running when inspected
on 2 November 2026. Q101 field-comparison run 41450637089 had passed.
This checkpoint binds the frozen Haar slow-reset calculus to the CURRENT
Q102 fresh slow response. It derives conditional postslow scalar,
compensator and whole axial pressure-moment coefficients. It does not
derive the full pointwise pressure or final radial residual bound.

## 1. The old Haar mean contains two slow histories

Let S0 be the original initial slow field, S1 the first-real refreshed
slow increment inside Q86 FinalMeans, F0/F1 the corresponding complete
zero-Haar fast fields, and F2 the Q101 newly sourced fast increment.
The new wave changes covariance, not mean velocity. Immediately BEFORE
the Q102 fresh slow increment S2 is applied,

    O=S0+S1+F0+F1+F2,
    Haar(O)=S0+S1, Haar(O S2)=(S0+S1) S2.                (1)

S2 is auxiliary-independent. The original time/compact inverse identities
give zero Haar of ALL complete fast components, including their axial
A1 velocity corrections. They do not permit deleting fast correlations
in Q102's earlier postfast update, which has already formed those products
before Haar. Equation (1) applies only to products with the new slow field.

Retrieve S0 from the same initial.slow_fields(N,False), and S1 from the
same CurrentScalars.slow(N,False). A safe unweighted coefficient is

    O_H,N=243 [S0_N+S1_N], N=n+2.                        (2)

The 243 pads raw maxima into the five-variable factorial-through norm
and is harmless for already-normalized inputs. Q41.old_haar is correct
for its original FIRST-real reset, where only S0 preceded that slow
solve. Reusing it alone now would omit S1. No historical evidence or old
implementation is changed; this is a new stage binding.

S0 has retained tangential/radial grades 1,2 from the sharp initial scalar
constructor. S1 has grades alpha,alpha+1 with alpha=2-3kappa>1 for the
retained 0<kappa<1/4. Thus (2) still bounds the current Haar mean at
grades 1,2. A full old-mean bound at B=1-kappa would be weaker and would
lose the scalar gain below. A zero old Haar bound would be false.

## 2. Exact linear solve followed by the complete remainder

Q102 constructs S2 from the CURRENT postfast scalars, which retain the
old scalar state, complete sourced covariance change ONCE, and all new
fast correlations. The original five-map identities cancel those scalar
defects in the linear reserved solve. This is an analytical identity
of the underlying signed physical fields, not a cancellation performed on
positive coefficient dictionaries. The code derives a remainder majorant
from the resulting S2 norm program; it does not set actual pressure or
scalars to zero.

For old Haar norms (ox,ou,oz), new slow norms (x,u,z), the original radial
base b with its extra epsilon, and tau=(n+1)+(n+1)(n+2), define

    S=2(b+ox)x+x^2, H=(b+ox)z+oz*x+x*z,
    T=2ou*u+u^2, K=oz*u+ou*z+u*z, Zf=2oz*z+z^2.

On the WHOLE reserved potential hull [a,bh], put
M=(bh-a)(a^-1+1+bh+bh^2). The inherited sufficient remainder is

    P_after <= M[tau*x+S+T+(n+1)H],
    Jtheta_after <= M K,
    Jz_after <= M Zf+(M/2)[tau*x+S+T+(n+1)H].             (3)

The same original full-hull G=0, nondegenerate five-map, compensator,
velocity-moment and smooth compact-extension identities remain required.
No check of dictionary keys, finite samples or separate bump supports
proves them. The norm hull is explicitly checked to contain the entire
reserved potential hull. M here may use the reserved hull because EVERY
remainder integrand contains S2 or one of its derivatives. Q102's full
fast scalar integrals could not be restricted to that smaller hull.

Whole compact traces are essential for the radial vector viscosity:

    int partial_R(partial_R x+x/R) dR=0,
    int R^2 partial_R(partial_R x+x/R) dR=0.              (4)

These are integral identities with compact endpoint jets, not deletion of
pointwise viscosity. The remaining slow time and axial viscous terms are
bounded conservatively by tau*x at inputs through n+2, dropping beneficial
epsilon factors. Both full pressure contributions to Jz and all nonlinear
old/new slow correlations remain in (3).

## 3. Gain and the exact G=0 support gate

For C=3/2-kappa, S2 has grades C,C+1. Old Haar grades 1,2 and b's
extra epsilon give the possible worst remainder grades

    min(C+1,2C,C+2)=C+1=5/2-kappa.                      (5)

Here C>1, so 2C>C+1. This is a conditional NORMALIZED scalar gain over
Q102's postfast C bound. It is not complete residual improvement or
physical absolute decay. Original initial/Haar and native norm contracts
are retained, not established by summing the coefficient histories.

If G is not zero on the full potential hull, the omitted linear scalar
leaks include int R^2 G u_s dR and int 2R G z_s dR. At a generic bounded
G they have only grade C. An alternative to exact G=0 would have to bound
the FULL missing terms (including G x_s in the radial-axial flux) at
grade C+1 with all mixed jets, or prove an appropriate combined exact
cancellation. The current wrapper does not supply such an alternative.
A compact pure swirl comparison with G(R)=R and ZERO angular velocity
moment has a nonzero first leak while the G=0 remainder K is zero. For
any nonnegative smooth compact phi of nonzero width, take u=phi(R)(R-c),
c=(int R^3 phi)/(int R^2 phi). The leak is the positive
weighted variance int R^4 phi-(int R^3 phi)^2/(int R^2 phi). The axial
velocity moment is also zero since x=z=0. This demonstrates the missing
information in a generic G-free estimate, not canonical five-map
realizability or failure of the actual construction. A constant-G example
would be weaker because its angular leak vanishes under that velocity
moment constraint.

## 4. Pressure scalar consequences, with memory preserved

For the SAME moving original rho and c_rho=(int R^2 rho dR)/2, (3) gives

    ||rho P_after||_n <=243 weighted_rho(n) P_after,n,
    ||partial_Z(Jz_after+c_rho P_after)||_n
       <=(n+1)[Jz_after,n+1+243 c_rho(n+1) P_after,n+1].   (6)

The first is a weighted bound for the COMPENSATOR COMPONENT at grade
C+1. The second acquires the original epsilon when used in the whole
axial scalar row, giving grade C+2=7/2-kappa. The derivative hits the
whole product; c_rho,Z P is retained. The code obtains postslow scalars
at n+1 rather than freezing the profile or using only the latest wave
pressure increment.

Small P_after does not bound the auxiliary nonzero-Haar pressure. Full
pm remains the original shifted T0 applied to the complete FINAL radial
source minus rho P_after, with its local terms, Pi1/Pi2/base archives and
earlier pressure memory in their original stages. No statement that total
pressure, A0 source or the radial residual has grade C+1 follows from (6).
Their domains, interactions, cutoffs and viscosity must be recomputed on
the final mean/covariance state. Haar moments alone are insufficient.

## 5. Requests and next obligation

Postslow scalar output n requests S2 and both old slow histories through
n+2. Q102 then needs covariance n+11, Q101 rows n+9, Q94 mean factory
n+14, conservative target/profile n+28,n+41, and new cutoff shape n+11.
The whole axial output n requires reset n+1, raising these requests to
covariance n+12, rows n+10, mean factory n+15, target/profile n+29,n+42
and cutoff shape n+12. These are newly composed finite bounds; the older
compact-tail and original primitive DAGs remain separate and unexpanded.
Neither giant physical jets nor uniform all-scale constants were evaluated.

Next assemble the complete final means and actual final radial source,
rebuild its full pressure and compact remainder at the required poles,
and join the complete mean and nonzero-angular residual ledgers. The
actual native/trace/support premises, unconditional pulse estimates,
complete residual closure, infinite realization and blow-up remain unproved.
No original root, field, width C_setup^-20, accepted domain or frozen
historical record is changed. State objects remain unchanged: these are
derived conditional norm recipes for the indicated postslow stage.

Exact controls distinguish both slow histories from the stale first-real
Haar input, verify (3), test the full-hull G gate and the integral-only
viscosity cancellation, and retain the next moving pressure jet in (6).
Manufactured lower norm leaves test state, reserve, unweighted-stage and
profile requests; they do not certify canonical continuum fields. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json`.
