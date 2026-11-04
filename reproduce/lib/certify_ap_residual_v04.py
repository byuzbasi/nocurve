#!/usr/bin/env python3
"""Exact dyadic-rational stationarity/readout check. No fitting or dependencies."""
import hashlib
import json
import math
import sys
from fractions import Fraction


def floats(a):
    values=[float.fromhex(v) for v in a]
    assert all(math.isfinite(v) for v in values)
    return values


def integers(values):
    pairs=[v.as_integer_ratio() for v in values]
    powers=[d.bit_length()-1 for _,d in pairs];e=max(powers,default=0)
    return [v<<(e-k) for (v,_),k in zip(pairs,powers)],e


def upper(value):
    if not value:return 0.0
    return math.nextafter(float(value),math.inf)


def certify(a):
    n,p=a['n'],a['p'];assert 1<=n<=256 and 1<=p<=96
    X,ex=integers(floats(a['X']));S,es=integers(floats(a['S']))
    Y,ey=integers(floats(a['y']));R,er=integers(floats(a['rhs']))
    hi,lo=floats(a['hi']),floats(a['lo']);pred=floats(a['prediction'])
    c=float.fromhex(a['multiplier']);tol=float.fromhex(a['tolerance'])
    assert c>0 and math.isfinite(c) and tol==1e-8
    assert len(X)==n*p and len(S)==p*p and len(Y)==n and len(R)==p
    assert len(hi)==len(lo)==p and len(pred)==n
    cn,cd=c.as_integer_ratio();ec=cd.bit_length()-1

    def normal(parts,symmetric=True):
        B,eb=integers(parts);B=[B[i]+B[p+i] for i in range(p)]
        # Column-major matrix inputs preserve the original R values exactly.
        F=[sum(X[i+j*n]*B[j] for j in range(p)) for i in range(n)];ef=ex+eb
        ev=max(ey,ef);v=[(Y[i]<<(ev-ey))-(F[i]<<(ev-ef)) for i in range(n)]
        legacy_pen=[cn*sum(S[j+k*p]*B[k] for k in range(p)) for j in range(p)];ep0=ec+es+eb
        pen=([cn*sum((S[j+k*p]+S[k+j*p])*B[k] for k in range(p)) for j in range(p)]
             if symmetric else legacy_pen);ep=ep0+int(symmetric)
        eg=max(ex+ev,ep)
        grad=[(sum(X[i+j*n]*v[i] for i in range(n))<<(eg-ex-ev))-
              (pen[j]<<(eg-ep)) for j in range(p)]
        # Retain the original normalization, including its stored S operator.
        ed=max(er,ep0)
        denominator=(1<<ed)+max(max(abs(v) for v in R)<<(ed-er),max(abs(v) for v in legacy_pen)<<(ed-ep0))
        numerator=max(abs(v) for v in grad)
        q=Fraction(numerator<<max(0,ed-eg),denominator<<max(0,eg-ed))
        return q,F,ef

    q,F,ef=normal(hi+lo);rounded,_,_=normal(hi+[0.0]*p)
    legacy,_,_=normal(hi+lo,symmetric=False)
    P,ep=integers(pred);ed=max(ep,ef)
    delta=max(abs((P[i]<<(ed-ep))-(F[i]<<(ed-ef))) for i in range(n))
    readout=Fraction(delta,1<<ed)
    readout_limit=Fraction.from_float(8*sys.float_info.epsilon)*(1+max(Fraction.from_float(abs(v)) for v in pred))
    normal_pass=q<=Fraction.from_float(tol);readout_pass=readout<=readout_limit
    return dict(valid=normal_pass and readout_pass,normal_pass_exact=normal_pass,
        normal_equation_residual_upper=upper(q),rounded_coefficient_normal_upper=upper(rounded),
        legacy_operator_normal_upper=upper(legacy),legacy_operator_pass_at_same_tolerance=legacy<=Fraction.from_float(tol),
        readout_error_upper=upper(readout),readout_pass_exact=readout_pass,
        arithmetic='exact Python integers and dyadic fractions of hexadecimal binary64 inputs',
        coefficient_representation='hi + lo as an exact real sum',tolerance=tol,n=n,p=p,
        penalty_representation='exact symmetric part; identical quadratic form; original S normalization retained',
        multiplier=c,optimization_executed=False)


if __name__=='__main__':
    raw=sys.stdin.buffer.read();a=json.loads(raw);r=certify(a)
    r['input_sha256']=hashlib.sha256(raw).hexdigest()
    print(json.dumps(r,allow_nan=False))
