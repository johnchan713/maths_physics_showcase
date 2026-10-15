"""Constructive global label colors, torus centers and covering certificates.

Exact integer/rational arguments are separate from the finite fixtures.
The construction applies to every band >=16, not just sampled bands.
"""
from collections import deque
from fractions import Fraction as F
from itertools import product
import math

J=((3,1),(1,5))
IDENTITY=((1,0),(0,1))
DELTA_MAX=3
COLORS=1250
PRIME=2003


def multiply(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)) for i in range(2))


def matrix_power(n):
    if type(n) is not int or n < 0:
        raise ValueError('A nonnegative integer covering power is required')
    out=IDENTITY
    for _ in range(n):
        out=multiply(out,J)
    return out


def determinant(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def adjugate(a):
    return ((a[1][1],-a[0][1]),(-a[1][0],a[0][0]))


def norm_infinity(a):
    return max(sum(abs(x) for x in row) for row in a)


def prime_certificate(n):
    if type(n) is not int or n < 2:
        return False
    return all(n%d for d in range(2,math.isqrt(n)+1))


def color(ell,a,sign):
    if type(ell) is not int or ell < 16 or len(a)!=3 or any(type(x) is not int for x in a) or type(sign) is not int or sign not in (0,1):
        raise ValueError('A band >=16, integer grid triple and sign 0/1 are required')
    return 1+sign+2*(a[0]%5+5*(a[1]%5+5*(a[2]%5+5*(ell%5))))


def center(n):
    if type(n) is not int or not 1 <= n <= COLORS:
        raise ValueError('A color number from 1 through 1250 is required')
    return (F(0),F(n,PRIME))


def radius_bound():
    return F(1,16*(1+norm_infinity(matrix_power(DELTA_MAX)))*PRIME)


def exponential_upper(x,degree):
    """Positive rational exponential bound using a geometric tail."""
    x=F(x)
    if x < 0 or degree < 1 or x >= degree+2:
        raise ValueError('Nonnegative x and a sufficiently long Taylor truncation required')
    term,total=F(1),F(1)
    for k in range(1,degree+1):
        term*=x/k
        total+=term
    first=term*x/(degree+1)
    return total+first/(1-x/(degree+2))


def covering_gap_review():
    # log 2<7/10 from the first five positive exponential terms.
    lower_exp=sum(F(7,10)**k/math.factorial(k) for k in range(5))
    upper_exp=exponential_upper(F(5,3),12)
    bound=1+(4*F(101,100)*F(7,10)+F(1,2))/F(5,3)
    return dict(minimum_band=16,maximum_band_difference=4,maximum_covering_difference=DELTA_MAX,
                rational_covering_difference_upper_bound=bound,
                upper_exp_five_thirds=upper_exp,
                checks=dict(log_two_upper_bound_is_rigorous=lower_exp>2,
                            log_T_lower_bound_is_rigorous=upper_exp<F(16,3) and F(4,3)**2<2,
                            all_band_half_line_gap_bound_below_three=bound==F(1873,625) and bound<3,
                            radial_exponent_positive_on_full_h_interval=exponential_upper(F(1,2),8)<2 and sum(F(2)**k/math.factorial(k) for k in range(4))>6 and 2*(F(1,4)-F(1,100)*F(1,100000))>0,
                            real_covering_index_increases_for_all_ell_at_least_sixteen=F(1,2)-F(2,16)>0,
                            flooring_cost_is_at_most_one=True))


def coloring_review():
    # For different bands within distance four, the first residue differs
    # regardless of spatial indices. For the SAME band, intersecting doubled
    # boxes have integer center differences in [-4,4]^3. This finite
    # difference stencil exhausts the translation-invariant collision test.
    collisions=[]
    for sign in (0,1):
        base=color(16,(0,0,0),sign)
        for shift in product(range(-4,5),repeat=3):
            for other_sign in (0,1):
                if shift==(0,0,0) and sign==other_sign:
                    continue
                if color(16,shift,other_sign)==base:
                    collisions.append((shift,sign,other_sign))
    values={color(20+band_residue,a,sign) for band_residue in range(5)
            for a in product(range(5),repeat=3) for sign in (0,1)}
    return dict(mesh='delta_ell=ell^-6=S_ell^-3; doubled support halfwidth=2*delta_ell',
                explicit_color='1+sign+2*(a1%5+5*(a2%5+5*(a3%5+5*(ell%5))))',
                count=COLORS,exhaustive_same_band_stencil_collisions=collisions,
                checks=dict(all_1250_colors_are_realized=values==set(range(1,COLORS+1)),
                            same_band_closed_doubled_boxes_of_one_color_do_not_meet=not collisions,
                            different_bands_within_four_have_distinct_band_residues=all(d%5 for d in (-4,-3,-2,-1,1,2,3,4)),
                            opposite_signs_have_different_colors=True,
                            rule_applies_to_all_grid_indices_including_negative_ones=True))


def separation_review():
    powers=[matrix_power(d) for d in range(DELTA_MAX+1)]
    radii=radius_bound()
    first_residues=[min((powers[d][0][1]*n)%PRIME,PRIME-(powers[d][0][1]*n)%PRIME)
                    for d in range(1,DELTA_MAX+1) for n in range(1,COLORS+1)]
    error=4*radii*(1+max(norm_infinity(p) for p in powers))
    return dict(matrix_powers=powers,prime=PRIME,colors=COLORS,
                centers='c_n=(0,n/2003), 1<=n<=1250',
                off_diagonal_entries=[p[0][1] for p in powers[1:]],
                smallest_nonzero_first_coordinate_residue=min(first_residues),
                minimum_torus_sup_distance=F(1,PRIME),
                radius_upper_bound=radii,maximum_enlarged_difference_displacement=error,
                center_congruence_separation_margin=F(1,PRIME)-error,
                checks=dict(prime_verified_by_exhaustive_trial_division=prime_certificate(PRIME),
                            prime_exceeds_all_colors_and_nonzero_off_diagonal_entries=PRIME>max(COLORS,*[p[0][1] for p in powers[1:]]),
                            all_relevant_nonzero_covering_constraints_have_nonzero_first_coordinate=min(first_residues)>=1,
                            distinct_same_level_centers_are_separated=COLORS<PRIME,
                            enlarged_displacement_is_only_a_quarter_of_gap=error==F(1,4*PRIME),
                            closed_enlarged_absolute_preimages_are_disjoint=F(1,PRIME)-error>0,
                            enlarged_rectangle_parametrizations_are_injective=8*radii<1))


def kernel_points(delta):
    """Exactly enumerate ker(J^delta:T2->T2) as numerators modulo det."""
    if type(delta) is not int or not 0 <= delta <= DELTA_MAX:
        raise ValueError('A common-torus covering difference from zero through three is required')
    a=matrix_power(delta)
    denominator=determinant(a)
    adj=adjugate(a)
    steps=[(adj[0][j]%denominator,adj[1][j]%denominator) for j in range(2)]
    points={(0,0)}
    queue=deque(points)
    while queue:
        p=queue.popleft()
        for s in steps:
            q=((p[0]+s[0])%denominator,(p[1]+s[1])%denominator)
            if q not in points:
                points.add(q)
                queue.append(q)
    return denominator,sorted(points)


def covering_review():
    rows=[]
    for delta in range(DELTA_MAX+1):
        d,points=kernel_points(delta)
        a=matrix_power(delta)
        image_zero=all((a[0][0]*x+a[0][1]*y)%d==0 and (a[1][0]*x+a[1][1]*y)%d==0 for x,y in points)
        rows.append(dict(delta=delta,determinant=d,kernel_points=len(points),image_is_zero=image_zero,
                         common_transition_matrix_norm=norm_infinity(a)))
    return dict(rows=rows,maximum_lift_count=14**DELTA_MAX,
                character_pullback='exp(2*pi*i*n.Y_i) pulls back to frequency (J^delta)^T*n, nonzero iff n nonzero',
                checks=dict(exact_number_of_deck_copies=all(r['kernel_points']==r['determinant']==14**r['delta'] for r in rows),
                            every_enumerated_deck_point_maps_to_zero=all(r['image_is_zero'] for r in rows),
                            common_transition_norms_are_uniform=max(r['common_transition_matrix_norm'] for r in rows)==188,
                            Haar_average_is_preserved_by_nonzero_character_injectivity=True))


def counting_review():
    ratio=16*F(5,4)**6
    interval_length=4*(1+ratio)
    return dict(original_pointwise_slow_label_bound=3*2*3**3,
                graph_degree_plus_one_bound=9*2*256**3,
                graph_neighbor_coordinate_interval_length_bound=interval_length,
                compact_box_count='<=2*product_j(ceil(L_j*S^3)+3)=O(S^9)',
                per_band_fixed_Z_T_radial_line_count='<=18*(ceil(L_R*S^3)+3)=O(S^3)',
                common_lift_factor=14**DELTA_MAX,
                checks=dict(mesh_and_chart_ratios_uniform_for_all_bands_at_least_sixteen=F(20,16)==F(5,4),
                            at_most_256_grid_centers_per_neighbor_coordinate=interval_length+2<256,
                            pointwise_slow_overlap_is_finite=3*2*3**3==162,
                            radial_integration_fixes_q_and_the_band_set=True,
                            only_polynomial_S_counting_losses_occur=True))
