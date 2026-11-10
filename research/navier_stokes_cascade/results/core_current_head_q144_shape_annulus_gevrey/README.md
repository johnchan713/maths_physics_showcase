# Q144: shape/annulus mixed recurrence and the missing gate

Parent: Q143, `b3d49e6478b3b8475f9d7f58749d426acf796559`.
The live branch, Q143 evidence and CI, Q63, Q67, Q70, Q73, Q127 and
Q140 were read first. Historical evidence, fields, roots, supports,
pressure datum, frequency, scales and the width `t=C^-20` are unchanged.

Q143 closed the post-taper interval through `Xi=110`. The next region is
the explicit shape, axial restoration and earlier annular join. Q70 gives
all radial orders only through angular order six, while Q127 gives all
angular orders of the same earlier root only at radial order zero. Neither
statement alone is the mixed four-field input required by Q140. This
checkpoint applies the product and composition recurrences directly to
the explicit formulas, retaining Q143's actual endpoint memory.

## A real shear floor that Q70 did not record

On the shape,

    (log E)_y=1/10+sigma'(y/T)(log f-ell)/T,
    a=1-2(log E)_y.

Q70 has `||sigma'||<=64`, `|ell|<=J6<=B0/100`, `|log f|<=6` and
`T=2000 B0`, with `B0>=1`. Therefore

    a >= 4/5-128(J6+6)/T
      >= 4/5-128(1/100+6)/2000
       =1298/3125 >2/5.                                (1)

The fixed axial restoration changes `U`, not `E`, so its floor is `4/5`.
On the earlier annular join put `q=sum_(l=1)^3 e_l beta_l`. Q70's accepted
real root cap is `b=5.72e-14`, `|beta_l|<=64`, and
`|beta_l'|<=L1=172167120000`. Hence

    |q|<=3b*64,  |q_y|<=3b*L1,
    a=4/5-2q_y/(1+q)
      >=144709422735784/195312499997855 >0.7409.       (2)

The common floor on this entire region is thus the value in (1). The
larger complex root cap from Q127 is deliberately not used in this real
acceptance inequality.

## Fixed mixed recurrence with complete memory

For the shape, retain the actual Q143 endpoint envelopes for `E,U` and
the complex order-zero cap `L0` for `ell`. The Q63 step satisfies a fixed
radial Gevrey-3 envelope after division by `T`; `ell` and `log f` have a
fixed angular base on the Q73 strip. Differentiate the displayed logarithm
and compose the exponential at the supplied order-zero value cap. No
high-order seminorm is exponentiated.

For the earlier join Q127 supplies a genuine complex radius `sigma` for
the same eta-dependent root, including every derivative of its varying
quadratic coefficient. Its coefficients obey a fixed angular envelope.
The Q70 width-`.08` radial bumps obey

    |beta^(r)| <=(162*3260)[8*3260/.08]^r(r!)^3,

including `r=0` after enlarging the prefactor. Products with the ideal
field therefore give one fixed mixed Gevrey-3 envelope for `E,U,1/E` on
the shape, restoration and earlier annulus.

Convert `E=sqrt(2X)F` using the rational caps
`sqrt(2X)<=2R`, `1/sqrt(2X)<=1`. Starting with the five ACTUAL Q143
moments at `Xi`, reintegrate the unchanged densities to `R`. The full Q67
state recurrence then bounds `a,b,p2,E`; `p2` still contains

    X{4 Aphys eta Pi-(1-eta^2)Pi_eta}/(L E),

and the `M,I,J,S,Pi` rows are neither reset nor copied from the reference.
The logarithmic length is certified by `length>=T+10` and
`R/Xi<=2^length`; no floating-point logarithm is used.

## What remains unproved

The mixed recurrence and positive floor are conditional on binding the
supplied Q73/Q127 complex order-zero caps to the same selected profile.
Before Q140--Q138 can be applied on this region one must also prove

    G>=max(2,|mu|,|p2|,|mu p2|),
    Eexp>=exp(G),  V0<=the unchanged integer N.             (3)

These are now the precise regional gate, rather than a missing shear
floor or missing pressure row. The tests check the exact floor arithmetic,
same-root radius, endpoint/memory dependence, pressure response, guards
and byte-identical evidence. They are diagnostics, not a continuum proof.
The whole compact maximum, final repair/cone join, uniform physical pulses,
complete covariance means, full residual improvement and blow-up remain
unproved.
