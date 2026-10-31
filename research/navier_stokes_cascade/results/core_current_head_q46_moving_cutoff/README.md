# Q46: moving cutoff jets in the native compact chart

Parent: `dfdc6a763cb05140b90f0e876ac913f2d652cc62` (Q45). On 31 October
2026 its focused CI 41170467084 passed; broad CI 41170467028 was running.
Q45 left three native raw cutoff norm inputs. For the fixed physical
choice these follow from one one-dimensional profile ledger. No new cutoff
is selected and no archived estimate or evidence is overwritten.

On the original shell write chi=hat_chi(X), X=U^(2/d)/(2s(Z,T)).
Let H_k bound sup |hat_chi^(k)|/k! on the closed shell including the
plateaux, and use Q44's u=min(1,Rmin^3). For j>=1 put

    L_j=2/(d j u^j)+(1/3) sum_(i=0)^j 32^i 64^(j-i).
    B_0=Xb,   B_j=sum_(i=1)^j (i/j) L_i B_(j-i).

L_j bounds the sum of absolute normalized mixed jets of log X of exact
degree j: the radial contribution is the derivative of (2/d)log U and
the other contribution is the established two-variable log s bound in
core_primary_target_cutoff_budget, on its accepted physical scale branch.
There are no auxiliary-angle derivatives. The formal Taylor identity
X(point+increment)=X(point) exp(increment of log X) implies the B
recurrence by coefficient comparison; products include radial/slow mixed
terms even though log X separates. Only finite Taylor polynomials through
the requested degree are needed. No analyticity of hat_chi is assumed.

Let V(t)=sum_(j>=1) B_j t^j. For each i, define the raw bound

    D_i=i! sum_(k=0)^i H_k [t^i] V(t)^k.

Taylor composition with the normalized outer seminorms bounds the sum of
absolute normalized jets at degree i by D_i/i!. Multiplication by i!
therefore bounds every raw mixed partial of that degree. In particular
three sufficient Q45 inputs at order n are

    chi:       max_(0<=i<=n) D_i
    1-chi:     max({1+H_0} union {D_i:1<=i<=n})
    chi_U:     max_(1<=i<=n+1) D_i.

The positive-degree restriction in the last expression excludes the
constant profile term. It does not restrict derivatives to radial ones.
The output derivatives of chi_U are mixed derivatives of total order at
most n+1 containing at least one U slot. Using all such degrees is safe.
Outside the transition the positive-order derivatives vanish by Q30's
plateau identities. Constant extensions retain the same bounds (H_0>=1
for an actual cutoff with outer plateau one). The Xb estimate is needed
only on the shell. Nothing is differentiated at a translated endpoint;
Q29's fixed-shift proof remains unchanged.

The factory uses the same Xb, d, width and R hull for these constants
and Q44/Q45. It requires H through compact order+1, hence through m+7
for Q41 outputs through m. The original time inverse remains separate.
Every H_k exists for the specified fixed smooth shape, but height and
plateau witnesses do not determine its value. For example, centering a
smooth zero-to-one step of width 1/N strictly inside (c0,c1) preserves
both witnesses and height one while its maximum slope grows like N.
These are comparison profiles, not changes to the original construction.
A numerical certificate needs independently justified seminorms of the
actual shape; finiteness at each fixed order is all used here.

This supplies the moving cutoff inputs for operator boundedness. Next
bound the initial axial A1 remainder from its same-state higher source
jets, preserving pressure and the full shifted moment. No compact tail
gain, uniform actual pulse, complete residual improvement or blow-up
follows from the cutoff estimate. Tests exercise finite jet algebra and
the full provider interfaces, not the physical high-order continuum field.
