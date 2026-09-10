"""Positive backward stress and a stable factorization at the outer edge.

These formulas require the exact total moments of the compensated reference.
They do not construct a regular axis or cancel the complete PDE residual.
"""
import mpmath as mp
from construction import Rule, step, step_prime, heat_coefficients


class Exterior:
    def __init__(self,h='.005',xtail=1000,amplitude='1e-10',order=32,panels=8,degree=32):
        self.h,self.xtail,self.amplitude = map(mp.mpf,(h,xtail,amplitude))
        self.co = mp.mpf('.001')
        self.rho = self.co*self.h
        self.rule = Rule(order,panels)
        self.b = heat_coefficients(self.h,degree)
        self.degree = degree
        v,w = mp.gauss_quadrature(order,'laguerre')
        self.laguerre = list(zip(v,w))

    def state(self,y,eta):
        y,eta = mp.mpf(y),mp.mpf(eta)
        X = self.xtail*mp.exp(y)
        r = mp.sqrt(2*X)
        Z = 2*(1-eta*eta)/X
        H = 1+self.h*mp.fsum(self.b[n]*Z**n for n in range(1,self.degree+1))
        Hp = self.h*mp.fsum(n*self.b[n]*Z**(n-1) for n in range(1,self.degree+1))
        K = self.amplitude*mp.exp(-(mp.mpf('.5')+self.h)*y)*H
        fo = 1-self.rho*step((3-y)/2)
        fp = self.rho*step_prime((y-1)/2)/2
        a = 2+2*self.h+2*Z*Hp/H-2*fp/fo
        return dict(X=X,r=r,Z=Z,H=H,Hp=Hp,K=K,fo=fo,fp=fp,a=a,
            radial_weight=2+2*self.h+2*Z*Hp/H)

    def direct(self,y,eta):
        y,eta = mp.mpf(y),mp.mpf(eta)
        if not mp.mpf('.5')<=y<3 or abs(eta)>1:
            raise ValueError('Interior terminal coordinates required')
        s = self.state(y,eta)
        L = 1-2*self.h*eta*eta
        pieces = [mp.mpf(0),mp.mpf(0),mp.mpf(0)]
        left = max(y,mp.mpf(1))
        for x,w in self.rule.nodes:
            v = left+(3-left)*x
            q = self.state(v,eta)
            density = q['fp']/self.rho
            weight = (3-left)*w*density
            pieces[0] += weight*q['r']*q['K']*q['radial_weight']
            pieces[1] += weight*q['r']**3*q['K']
            pieces[2] += weight*s['r']**2*mp.expm1(v-y)*q['K']**2*q['fo']
        theta = 2*s['K']/s['r']*s['fp']/self.rho+pieces[0]/s['r']**2+pieces[1]/(2*L*s['r']**2)
        axial = eta*pieces[2]/(L*s['r'])
        ratio = axial/theta
        return dict(y=y,eta=eta,theta_over_rho=theta,axial_over_rho=axial,ratio=ratio,
            ratio_upper=2*s['K'],a_minus_two=s['a']-2,
            directional_loss=(s['a']-2)*ratio**2,directional_gap=2-(s['a']-2)*ratio**2)

    @staticmethod
    def flat_kernel(u):
        """The positive coefficient in f'=rho exp(-4/u^2) u^-3 B(u)."""
        u = mp.mpf(u)
        if u==0: return 8*mp.e
        if not 0<u<2: raise ValueError('Flat collar coordinate must be below 2')
        c = 1-u/2
        A,B = mp.exp(-1/c**2),mp.exp(-4/u**2)
        return A*(8+u**3/c**3)/(A+B)**2

    def flat(self,delta,eta):
        """Remove exp(-4/delta^2) before integrating either stress component."""
        delta,eta = mp.mpf(delta),mp.mpf(eta)
        if not 0<delta<=mp.mpf('.5') or abs(eta)>1:
            raise ValueError('Require 0<delta<=.5 and |eta|<=1')
        y = 3-delta
        s = self.state(y,eta)
        L = 1-2*self.h*eta*eta
        integral1,integral2,integralz = mp.mpf(0),mp.mpf(0),mp.mpf(0)
        for w,weight in self.laguerre:
            root = mp.sqrt(1+delta**2*w/4)
            u = delta/root
            separation_over_delta3 = (w/4)/(root*(1+root))
            separation = delta**3*separation_over_delta3
            ratio = mp.expm1(separation)/separation if separation else 1
            q = self.state(3-u,eta)
            B = self.flat_kernel(u)
            integral1 += weight*B*q['r']*q['K']*q['radial_weight']
            integral2 += weight*B*q['r']**3*q['K']
            integralz += weight*B*ratio*separation_over_delta3*q['K']**2*q['fo']
        btheta = 2*self.rho*s['K']/s['r']*self.flat_kernel(delta)
        btheta += self.rho*delta**3/(8*s['r']**2)*(integral1+integral2/(2*L))
        bz = eta*self.rho*s['r']*integralz/(8*L)
        edge = self.state(3,eta)
        expected_theta = 16*self.rho*edge['K']*mp.e/edge['r']
        expected_z = eta*self.rho*edge['r']*edge['K']**2*mp.e/(8*L)
        return dict(delta=delta,eta=eta,btheta=btheta,bz=bz,ratio=delta**6*bz/btheta,
            ratio_over_delta6=bz/btheta,expected_btheta=expected_theta,expected_bz=expected_z,
            expected_ratio_over_delta6=eta*edge['X']*edge['K']/(64*L),
            positive_flat_factor=mp.exp(-4/delta**2),
            binary64_flat_factor=float(mp.exp(-4/delta**2)))

    def endpoint(self,eta):
        edge = self.state(3,eta)
        return dict(theta=0,axial=0,direction=[1,0],directional_gap=2,
            positive_angular_coefficient=16*self.rho*edge['K']*mp.e/edge['r'])
