"""Exact derivative requests for the common raw principal moment."""
def request(m,z,p,r):
    if any(type(x) is not int for x in (m,z,p,r)) or m<0 or z not in (0,1) or r not in (0,1) or p<3:
        raise ValueError('Invalid shell request')
    return dict(total=m+z+2*p-r,slow=m+z+p+1,profile=m+z+p+2,eta=m+z+1)

def record():
    return dict(status='shared-principal-pressure-moment',
                radial_extra_output=request(1,0,199995,0),axial=request(0,1,99995,0),
                shared_order='4-5*kappa',full_axial_order_improved=False,
                pressure_reset=False,actual_pulse_residual_or_blowup=False)
