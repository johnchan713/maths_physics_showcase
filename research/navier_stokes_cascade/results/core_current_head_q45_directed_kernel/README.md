# Q45: weighted directed kernels on the fixed cutoff shell

Parent: `c857803f1f83eed6328f1a3f78220079e1533bfa` (Q44).
Q44's list of open cutoff inputs needs a distinction. Q30 already proves
the existence of fixed plateau witnesses for the specified physical cutoff
in X. That existence does not require a new band-by-band hypothesis.
Quantitative bounds still require witnesses for that same profile and
its raw derivative majorants. Historical evidence is unchanged.

Let Xa<c0<c1<Xb be its plateau witnesses, x=log(X/Xa),
L=log(Xb/Xa), a=width^2 in (0,1], b=4 and p the requested nonnegative
integer pole. Write

    w(x)=exp(-a/x^2-b/(L-x)^2) min(1,x,L-x)^(-p).

Use Q30's rational gap g (no larger than either logarithmic plateau
width), lo=2(Xb-Xa)/(Xb+Xa)<=L<=hi=(Xb-Xa)/Xa. These logarithm
inequalities follow by differentiating at ratio >=1. Put

    r=min(1,lo/4,g,a/(1+p+a+b)).

This is only a proof collar for an integral estimate. It does not change
the field width, the transition, the pole, or the accepted bands.

**Global height.** If delta=min(1,x,L-x)<1, suppress the farther
exponential and set t=1/delta. Then log w <= p log t - min(a,b)t^2
<= p t - min(a,b)t^2 <= p^2/(4 min(a,b)). For delta=1, w<=1.
Thus E=exp(p^2/(4 min(a,b))) bounds the complete weight. The possibly
large exponential stays symbolic; there is no factorial-sized list.

**Directed integral.** For 0<x<=r and 0<y<=x, delta(y)=y.
The function exp(-a/y^2)y^(-p) is nondecreasing there because
p r^2<=a<2a. The ratio of the other exponential factors is at most
exp(b/(lo-r)^2). Thus the left integral of w divided by w(x) is at
most r exp(b/(lo-r)^2). The reflected right collar yields
r exp(a/(lo-r)^2), since b>=a. On [r,L-r], 1/w(x)<=exp((a+b)/r^2);
either integral ratio is at most hi E exp((a+b)/r^2).
The sum of these three bounds is a valid common K. The left branch is
multiplied by 1-chi or one of its derivatives and is absent on the
right collar; the right branch is absent on the left collar. K is NOT
a global bound for a single uncut one-sided primitive: for example its
left integral tends to a positive mass at the right edge while w tends to
zero there. Preserving the cutoff branches is essential.

On the transition only, 1/w<=exp((a+b)/g^2). These four finite constants
(hi,E,K,transition inverse bound) supply Q44's length, edge, kernel and
interior slots at its exact pole. Q29 already accounts for the characteristic
Jacobian after differentiating at fixed shift; there is no second Jacobian
or logarithmic length to insert into K. The calculation works for every
fixed requested pole, without an epsilon loss or a uniform-in-p claim.

Pass Shell(...).remaining(cutoff) to Q44.build on the same physical hull.
The three raw cutoff norms and original time inverse remain inputs. The
rational witnesses must apply to the original profile; the code checks
their ordering, not physical plateau values. The shifted axial A1 moment
and complete pressure memory remain in the existing equations.

Next establish the raw cutoff jets of the fixed shape composed with the
moving X coordinate, then connect the initial A1 integration-by-parts tail
to higher jets of the same source. Kernel boundedness gives no small
tail by itself. No uniform actual pulse, full residual gain or blow-up
is established. Finite integration tests are diagnostics for the formulas;
the continuum estimate above uses the collar split.
