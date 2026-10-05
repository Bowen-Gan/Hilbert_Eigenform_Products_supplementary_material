#!/usr/bin/env python3
"""Certify finite numerical endpoints in the manuscript's Eisenstein bounds.

Standard-library rational arithmetic only. Machin's identity and alternating
arctangent sums enclose pi; a positive exponential series and a geometric
remainder enclose exp(x). Every reported comparison is an exact rational
inequality. The cited discriminant/special-value bounds and the manuscript's
symbolic inductions remain mathematical inputs; finite sampling is not used
to prove an infinite range.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import argparse
import hashlib
import json
import sys

if not __debug__:
    raise SystemExit('Do not use optimized Python; assertions are required.')


def arctan_interval(reciprocal, terms=40):
    """Alternating-sum enclosure of atan(1/reciprocal)."""
    assert reciprocal > 1 and terms > 0
    x=Q(1,reciprocal)
    partial=sum(((-1)**j*x**(2*j+1)/Q(2*j+1)
                 for j in range(terms)),Q(0))
    next_sum=partial+(-1)**terms*x**(2*terms+1)/Q(2*terms+1)
    return min(partial,next_sum),max(partial,next_sum)


def pi_interval():
    a,b=arctan_interval(5)
    c,d=arctan_interval(239)
    # Machin's identity: pi=16 atan(1/5)-4 atan(1/239).
    lo,hi=16*a-4*d,16*b-4*c
    # Round outwards to a common decimal denominator. This remains an
    # exact rational enclosure and keeps high-power certificates compact.
    scale=10**40
    lower=lo.numerator*scale//lo.denominator
    upper=-((-hi.numerator*scale)//hi.denominator)
    return Q(lower,scale),Q(upper,scale)


def exp_interval(x, terms=100):
    """Enclose exp(x) for rational x>=0 with a geometric series tail."""
    x=Q(x)
    assert x>=0 and terms>=1 and x<terms+2
    total=Q(1);term=Q(1)
    for j in range(1,terms+1):
        term*=x/j;total+=term
    next_term=term*x/(terms+1)
    # Later ratios x/(j+1) are at most x/(terms+2)<1.
    return total,total+next_term/(1-x/(terms+2))


def root_discriminant_interval(degree):
    a,b=exp_interval(Q(83185,10000*degree))
    return Q(29099,1000)/b,Q(29099,1000)/a


def run():
    pi_lo,pi_hi=pi_interval()
    assert Q(3)<pi_lo<pi_hi<Q(22,7)
    assert pi_hi<Q(355,113)
    rows=[]
    def above(name, lower, threshold, section):
        lower,threshold=Q(lower),Q(threshold)
        assert lower>threshold,(name,str(lower),str(threshold))
        rows.append({'name':name,'direction':'>','certified_lower':str(lower),
                     'threshold':str(threshold),'section':section,'passed':True})
    def below(name, upper, threshold, section):
        upper,threshold=Q(upper),Q(threshold)
        assert upper<threshold,(name,str(upper),str(threshold))
        rows.append({'name':name,'direction':'<','certified_upper':str(upper),
                     'threshold':str(threshold),'section':section,'passed':True})
    equal='Equal-weight Eisenstein products in higher degree'
    distinct='Distinct-weight Eisenstein products in higher degree'
    quadratic='Quadratic Eisenstein products with larger narrow class number'
    common='General discriminant-bound endpoints used for Eisenstein--cuspidal products'

    # Quadratic EE: D>=12 and the two retained gamma-ratio factors.
    above('9720/pi^8',Q(9720)/pi_hi**8,1,quadratic)
    above('12*2^2/(4*pi^2)',Q(12)/pi_hi**2,1,quadratic)

    # Equal weights k>=3: all infinite ranges use the symbolic ratios in
    # the proof. These are their finite starting values and ratio constants.
    z3,z6=Q(5,4),Q(33,32)
    h3_lower=Q(15)/(16*z6*z3**2)
    assert h3_lower==Q(32,55)
    above('7*(32/55)^3',7*h3_lower**3,1,equal)
    above('(4^4/4!)*(32/55)^4',Q(4**4,factorial(4))*h3_lower**4,1,equal)
    above('64/55',Q(64,55),1,equal)
    above('Gamma(6)/((2*pi)^3*Gamma(3)*zeta6_upper*zeta3_upper)',
          Q(120)/(8*pi_hi**3*2*z6*z3),Q(1,6),equal)
    u3=Q(3**3,factorial(3))**6/6**3
    assert u3==Q(19683,512)
    above('U3=19683/512',u3,9,equal)
    above('64/6',Q(64,6),8,equal)
    above('7/pi',Q(7)/pi_hi,2,equal)

    # Equal weight two: the displayed root-discriminant and small-degree
    # decimal thresholds use a certified pi interval, not binary floats.
    c_hi=pi_hi**8/9720
    below('16*(pi^8/9720)^2',16*c_hi**2,Q(15248,1000),equal)
    above('29.099*exp(-8.3185/13)',root_discriminant_interval(13)[0],
          Q(15248,1000),equal)
    below('((4^3-1)*(pi^8/9720)^3)^2',(4**3-1)**2*c_hi**6,
          Q(343461,100),equal)
    below('((4^4-1)*(pi^8/9720)^4)^2',(4**4-1)**2*c_hi**8,
          Q(5362189,100),equal)
    sqrt725_lower=Q(2692582403,10**8)
    assert sqrt725_lower**2<725
    above('2*(4^4-1)*725^(3/2)*(3/pi^4)^4',
          2*(4**4-1)*725*sqrt725_lower*(Q(3)/pi_hi**4)**4,8,equal)
    t6=2*(4**6-1)*Q(6**6,factorial(6))**3*(Q(3)/pi_hi**4)**6
    above('T6',t6,Q(188,100),equal)
    above('2187/(16*pi^4)',Q(2187)/(16*pi_hi**4),1,equal)
    above('(17/8)^3',Q(17,8)**3,2**3+1,equal)
    below('17*pi^8/6480',17*pi_hi**8/6480,25,equal)
    assert 5**4<725
    above('29.099*exp(-8.3185/5)',root_discriminant_interval(5)[0],
          Q(5512,1000),equal)

    # Distinct weights k1>=4. The proof reduces to these finitely many
    # smallest pairs using factors>1, n>=3 and decreasing zeta(s).
    pi22=Q(22,7)
    r0=Q(7)/(4*pi22)
    z2=pi22**2/6;z4=pi22**4/90;z6_exact=pi22**6/945
    above('delta_cubic^3>=49>(7/2)^3',Q(49),Q(7,2)**3,distinct)
    above('29.099*exp(-8.3185/4)',root_discriminant_interval(4)[0],
          Q(7,2),distinct)
    above('2*r0 with pi replaced by22/7',2*r0,1,distinct)
    above('(6*r0^2/(zeta4*zeta2))^3',(6*r0**2/(z4*z2))**3,
          Q(11,10),distinct)
    above('(3*r0/(zeta4*zeta3_upper))^3',(3*r0/(z4*Q(5,4)))**3,
          Q(11,10),distinct)
    above('(20*r0^2/(zeta6*zeta4))^3',(20*r0**2/(z6_exact*z4))**3,
          100,distinct)
    above('(R-1)*S endpoint',Q(1,10)*100,1,distinct)

    # The (3,2) pair: even degrees>=6 start at6; degree4 has its own
    # exact local alternatives. We verify each numerical starting value.
    above('29.099*exp(-8.3185/6)',root_discriminant_interval(6)[0],7,distinct)
    above('(168/(5*pi^3))^6',(Q(168)/(5*pi22**3))**6,Q(3,2),distinct)
    above('147/(pi^2*zeta5_upper*zeta3_upper)',
          Q(147)/(pi22**2*Q(67,64)*Q(5,4)),10,distinct)
    assert 26**2<725
    above('16*725^2*26/(625*(22/7)^12)',Q(16*725**2*26)/(625*pi22**12),
          Q(1,3),distinct)
    above('725^2*(768/(335*pi^2))^4',725**2*(Q(768)/(335*pi22**2))**4,
          1000,distinct)
    above('725^3*(1152/(67*pi^5))^4',725**3*(Q(1152)/(67*pi22**5))**4,
          1000,distinct)
    for q in (2,4):
        assert Q(q**3+1,q-1)<=Q(65,3)
    below('65/3',Q(65,3),1000,distinct)
    assert (2**8+1,2**4+1,2**16+1)==(257,17,65537)
    assert (65537-17,65537-257)==(65520,65280)
    assert 65520-65280==240
    above('(240-65280/1000)/3',Q(240-Q(65280,1000),3),16,distinct)

    # Common finite exponential endpoints not covered by the older
    # special-value verifier. Their infinite degree extension is symbolic.
    above('29.099*exp(-8.3185/5)>51/10',root_discriminant_interval(5)[0],
          Q(51,10),common)
    assert Q(51,10)**4<725 and Q(51,10)**3<133
    assert Q(18,5)**3<49
    above('29.099*exp(-8.3185/8)',root_discriminant_interval(8)[0],
          Q(102,10),common)
    assert Q(7403,1000)**4>3003

    return {
        'schema':'eigenform-eisenstein-finite-bound-verification-v1',
        'status':'passed','completed':True,'passed':True,
        'python_version':sys.version.split()[0],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic':'exact rational interval arithmetic; no floating-point decisions',
        'pi_interval':{'lower':str(pi_lo),'upper':str(pi_hi),
                       'method':'Machin identity with40-term alternating arctangent enclosures'},
        'exp_interval_method':'100-term positive Taylor sum; geometric rational tail',
        'comparison_count':len(rows),'comparisons':rows,
        'external_mathematical_dependencies':[
            'Machin identity, alternating-series remainder theorem, and exponential Taylor expansion.',
            'Critical Hecke L-value bounds, nonvanishing and Galois equivariance.',
            'Product coefficient identities, algebraic-integrality/norm comparisons, and local Hecke recurrences.',
            'Minkowski discriminant bound and published minimum discriminants49and725.',
            'Unconditional Odlyzko--Takeuchi bound delta>29.099 exp(-8.3185/degree).',
            'Roblot class-number tables and ordinary Hilbert class-field discriminant identity.',
            'Zeta(s) decreases for s>1; integral-test upper bounds; exact zeta values at2,4,6.',
            'The manuscript proves integer-weight and integer-degree inductions by symbolic successive ratios.',
            'Parity and local residue-degree arguments reducing the distinct-weight endpoints.',
        ],
        'scope_limitations':[
            'No numerical sampling is treated as a proof of continuous monotonicity or an infinite range.',
            'The field tables, class numbers, analytic functional equations and classification theorem are not proved by this verifier.',
            'Special-value tables and EC conjugate cutoffs are checked by verify_special_values.py; lattice and cusp-space arithmetic have separate programs.',
        ],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        result=run()
    except Exception as error:
        result={'schema':'eigenform-eisenstein-finite-bound-verification-v1',
                'status':'failed','completed':False,'passed':False,
                'error_type':type(error).__name__,'error':str(error)}
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(result,indent=2)+'\n')
        print('FAIL: '+type(error).__name__+': '+str(error),file=sys.stderr)
        return 1
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('EISENSTEIN_FINITE_BOUNDS_PASS: comparisons='+str(result['comparison_count']))
    print('SCOPE: finite rational endpoints verified; analytic theorems and symbolic inductions remain explicit inputs.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
