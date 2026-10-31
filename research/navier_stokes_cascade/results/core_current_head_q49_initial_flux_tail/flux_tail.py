"""Original A1 tail via exact compact graph-flux cancellation."""
from pathlib import Path
import importlib.util
from functools import lru_cache

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q49_native',HERE.parent/'core_current_head_q48_original_axial_tail/native_tail.py')
old=importlib.util.module_from_spec(s);s.loader.exec_module(old)
full, moving, radial=old.full,old.moving,old.radial


class InitialTail:
    def __init__(self,current,max_row,*,d,pressure_cutoff,axial_cutoff):
        if not isinstance(current,full.refresh.maps.NewPoleMaps):
            raise TypeError('Current state must retain the original primary')
        full.order(max_row)
        b=current.target.source.haar.initial.budget
        self.b,self.maximum=b,max_row
        N=max_row+3
        if b.pole<12*N+18:
            raise ValueError('Initial covariance archive does not cover the rows')
        a,z=old.initial_shell(b.hull)
        self.chart=c=b.primary.chart
        def bind(data,needed):
            c0,c1,shape=data
            if len(shape)<=needed:raise ValueError('Fixed cutoff jet archive too short')
            w=b.primary.primary.root.width
            return moving.Moving(d,c.Rmin,c.Rmax,z,w,shape),moving.directed.Shell(a,z,c0,c1,w)
        self.press=bind(pressure_cutoff,max_row+2)
        self.ax=bind(axial_cutoff,max_row+1)

    @lru_cache(None)
    def __call__(self,n,output_pole):
        full.order(n);full.order(output_pole)
        if n>self.maximum:raise ValueError('Outside initial tail ledger')
        N=n+3;v=12*N+18;c=self.chart
        W=full.convert(self.b.primary.covariance_R(N)['W'])
        r=sum(c.Rmin**(-k-1) for k in range(N+1))
        g=2+6*c.B*sum((k+1)*c.Rmin**(-k) for k in range(N+1))
        D=n+4
        pm=old.kt(n+1,v,*self.press)*((D*g+D+2*r)*W)
        jets,shell=self.ax;geo=jets.geo
        f=geo.multipliers(n,1)
        raw=(2**(3*n)*(geo.Umax-geo.Umin)*f['k']*jets(n,v)['bounds']['chi_U']*f['f'])
        K=(raw*radial.chart.transfer(n,1/geo.d,geo.Rmax,geo.Umin,True)
           *radial.chart.transfer(n,geo.d,geo.Umax,geo.Rmin,False)
           *shell.bounds(v)['edge']*shell.bounds(0)['interior'])
        return dict(order=n,pole=0,Ez=2*D*K*(W+pm))


def build(current,m,*,d,pressure_cutoff,axial_cutoff,c0,c1,shape,time_inverse):
    """Initial cutoffs are bound separately from Q47's current compact maps."""
    full.order(m)
    tail=InitialTail(current,m+10,d=d,pressure_cutoff=pressure_cutoff,axial_cutoff=axial_cutoff)
    return full.build(current,tail,m,d=d,c0=c0,c1=c1,shape=shape,time_inverse=time_inverse)


def record():
    return dict(original_flux_annihilated=True,original_tail_rate='2-kappa',
                original_covariance_order='n+3',initial_pressure_shape='n+2',
                initial_axial_shape='n+1',output_pole=0,
                actual_lower_field_bounds_or_full_residual_certified=False)
