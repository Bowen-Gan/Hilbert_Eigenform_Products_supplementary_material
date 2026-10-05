#!/usr/bin/env python3
"""Exact finite certificate for Case 11 (Python 3, standard library only).

Run: python3 rq_certificate.py --output rerun/rq_certificate.json

Computational inputs are the discriminant inequalities proved in the paper.
The mathematical inputs used to interpret the certificates are Minkowski's
bound, the real-quadratic Bernoulli formula, and the discriminant/conductor
formula for a biquadratic genus field.  No floating-point arithmetic, tables,
class-number oracle, guessed special value, or numerical reconstruction occurs.

An ideal is represented by integer generators in the integral basis (1,w).
The gcd of its 2-by-2 minors is its lattice index.  Equality (gamma)=I is
certified by gamma in I and |Norm(gamma)|=[O:I].  This gives an ordinary class
group upper bound from all prime ideals of norm <= sqrt(D)/2.  For each of
D=40,60,85,165 the unramified biquadratic extension gives the lower bound 2.
"""

from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
import json
import sys

if not __debug__:
    raise SystemExit('Assertions must be enabled: do not use -O/-OO.')


def squarefree(n):
    return n > 0 and all(n % (p*p) for p in range(2, isqrt(n)+1))


def fundamental(D):
    return D > 1 and ((D % 4 == 1 and squarefree(D)) or
                     (D % 4 == 0 and D//4 % 4 in (2, 3)
                      and squarefree(D//4)))


def primes(bound):
    return [p for p in range(2, bound+1)
            if all(p % q for q in range(2, isqrt(p)+1))]


def kronecker(a, n):
    """Kronecker symbol (a/n), for integer a and positive integer n."""
    assert n > 0
    result = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        result *= 1 if a % 8 in (1, 7) else -1
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def Lminusone(D):
    assert fundamental(D)
    values = [kronecker(D, a) for a in range(1, D+1)]
    assert sum(values) == 0
    assert sum(a*c for a, c in enumerate(values, 1)) == 0
    return -Q(sum(a*a*c for a, c in enumerate(values, 1)), 2*D)


def zeta_minus_one(D):
    return -Lminusone(D)/12


class Ring:
    def __init__(self, D):
        assert fundamental(D)
        self.D = D
        self.t = D % 4
        self.c = (1-D)//4 if self.t else -D//4
        assert self.t*self.t-4*self.c == D

    def norm(self, z):
        a, b = z
        return a*a+self.t*a*b+self.c*b*b

    def mul(self, z, y):
        a, b = z
        d, e = y
        return (a*d-self.c*b*e, a*e+b*d+self.t*b*e)

    def small_primes(self):
        # Only prime ideals of residue degree one need a witness: an inert
        # rational prime generates its (principal) prime ideal itself.
        answer = []
        for p in primes(isqrt(self.D)//2+1):
            if 4*p*p > self.D:
                continue
            roots = [r for r in range(p)
                     if (r*r-self.t*r+self.c) % p == 0]
            answer.extend((p, r) for r in roots)
        return answer

    @staticmethod
    def ideal(P):
        p, r = P
        return [(p, 0), (-r, 1)]

    def product(self, I, J):
        return [self.mul(z, y) for z in I for y in J]

    @staticmethod
    def index(I):
        answer = 0
        for (a, b), (c, d) in combinations(I, 2):
            answer = gcd(answer, abs(a*d-b*c))
        assert answer > 0
        return answer

    def member(self, z, I):
        return self.index(I+[z]) == self.index(I)

    def generator(self, I):
        """Find, then independently check, a principal ideal generator.

        The search bound is only a way to find witnesses.  Completeness of
        the class-group certificate comes from Minkowski, not this search.
        Exhausting the search raises an error; it never means nonprincipal.
        """
        m = self.index(I)
        for b in range(10001):
            for target in (m, -m):
                disc = self.D*b*b+4*target
                if disc < 0 or isqrt(disc)**2 != disc:
                    continue
                for root in {isqrt(disc), -isqrt(disc)}:
                    if (root-self.t*b) % 2:
                        continue
                    z = ((root-self.t*b)//2, b)
                    if self.member(z, I):
                        assert abs(self.norm(z)) == m
                        return z
        raise AssertionError(f"No norm witness found for D={self.D}, I={I}")


def norm_minus_one_certificate(R):
    # A prime p=3 mod 4 dividing D forbids x^2-D*y^2=-4 modulo p.
    obstructions = [p for p in primes(R.D)
                    if p % 4 == 3 and R.D % p == 0]
    if obstructions:
        p = obstructions[0]
        squares = sorted({a*a % p for a in range(p)})
        assert (-4) % p not in squares
        return False, {"obstruction_prime": p, "squares_mod_p": squares}
    for b in range(1, 10001):
        d = R.D*b*b-4
        if d >= 0 and isqrt(d)**2 == d and (isqrt(d)-R.t*b) % 2 == 0:
            z = ((isqrt(d)-R.t*b)//2, b)
            assert R.norm(z) == -1
            return True, {"unit": z, "norm": -1}
    raise AssertionError(f"Unit-sign certificate not found for D={R.D}")


def class_certificate(D):
    R = Ring(D)
    Ps = R.small_primes()
    has_genus_extension = D % 5 == 0 and fundamental(D//5) and D//5 != 5
    relations = []
    if has_genus_extension:
        # Distinct fundamental discriminants 5, D/5, D define a biquadratic
        # field, and product= D^2 proves relative discriminant one.  All
        # three are positive, so the extension also splits at infinity.
        assert fundamental(5) and 5*(D//5) == D
        assert 5*(D//5)*D == D*D
        h = 2
        base = next(P for P in Ps if kronecker(5, P[0]) == -1)
        B = R.ideal(base)
        B2 = R.product(B, B)
        z = R.generator(B2)
        relations.append({"kind": "base_squared", "P": base,
                          "generator": z, "norm": R.norm(z),
                          "ideal_generators": B2})
    else:
        h, base, B = 1, None, None
    for P in Ps:
        I = R.ideal(P)
        assert R.index(I) == P[0]
        genus_value = (kronecker(D//5, P[0]) if P[0] == 5
                       else kronecker(5, P[0])) if h == 2 else 1
        if h == 2 and genus_value == -1:
            I = R.product(I, B)
            kind = "prime_times_base"
        else:
            kind = "prime_principal"
        z = R.generator(I)
        relations.append({"kind": kind, "P": P, "generator": z,
                          "norm": R.norm(z), "ideal_generators": I})
    negative, units = norm_minus_one_certificate(R)
    return {"D": D, "basis_polynomial": [R.c, -R.t, 1],
            "minkowski_bound_squared": str(Q(D, 4)),
            "prime_ideals": Ps, "base_ideal": base,
            "relations": relations, "class_number": h,
            "class_lower_bound": ({"genus_discriminants": [5, D//5, D],
                                    "biquadratic_discriminant": D*D}
                                   if h == 2 else {"trivial": True}),
            "negative_norm_unit": negative, "unit_certificate": units,
            "narrow_class_number": h if negative else 2*h}


def reduction_certificate():
    # The strict cubes avoid rounding both published decimal cutoffs.
    Ds = [D for D in range(2, 168) if fundamental(D) and
          ((D % 2 == 0 and D**3 < 600**2) or
           (D % 8 == 5 and D**3 < 2160**2))]
    assert 168**3 > 2160**2 and 72**3 > 600**2
    assert Ds == [5,8,12,13,21,24,28,29,37,40,44,53,56,60,61,
                  69,77,85,93,101,109,133,141,149,157,165]
    rows = [class_certificate(D) for D in Ds]
    survivors = []
    for row in rows:
        D = row["D"]
        if row["narrow_class_number"] == 1:
            row["decision"] = "narrow class number one"
            continue
        alphas = [4/zeta_minus_one(D)]
        if row["class_number"] == 2:
            alphas.append(4/(Lminusone(5)*Lminusone(D//5)))
        row["alpha"] = [str(a) for a in alphas]
        possible = [a for a in alphas if a.denominator == 1 and
                    (D % 2 or abs(a) >= 6)]
        if not possible:
            row["decision"] = "integrality / ramified lower bound"
        elif D == 40:
            # The indecomposability argument is already fully in the paper.
            R = Ring(40)
            P = R.ideal((3, 1))
            P2 = R.product(P, P)
            nu = (7, 2)
            assert R.member(nu, P2) and R.norm(nu) == R.index(P2) == 9
            assert 7 > 0 and 7*7 > 10*2*2
            row["D40_ideal_witness"] = {"P": [3,1], "nu": nu, "norm": 9}
            row["decision"] = "D=40 indecomposability argument in paper"
        else:
            survivors.append(D)
            row["decision"] = "retained for local computation"
    assert survivors == [12,21,24,28,69,77]
    expected = {
        12:(1,2,["24"]),21:(1,2,["12"]),24:(1,2,["8"]),
        28:(1,2,["6"]),40:(2,2,["24/7","10"]),44:(1,2,["24/7"]),
        56:(1,2,["12/5"]),60:(2,4,["2","5"]),69:(1,2,["2"]),
        77:(1,2,["2"]),85:(2,2,["4/3","5/2"]),93:(1,2,["4/3"]),
        133:(1,2,["12/17"]),141:(1,2,["2/3"]),165:(2,4,["6/11","5/6"])}
    actual = {r["D"]:(r["class_number"],r["narrow_class_number"],r["alpha"])
              for r in rows if "alpha" in r}
    assert actual == expected
    return {"arithmetic": "integer and rational; no external dependencies",
            "discriminants": Ds, "rows": rows, "survivors": survivors}


def divisor_data_integer(R, m):
    """Ideal divisor count and sigma_1 of the rational principal ideal (m)."""
    factors = []
    remaining = m
    for p in primes(m):
        v = 0
        while remaining % p == 0:
            remaining //= p
            v += 1
        if not v:
            continue
        roots = [r for r in range(p)
                 if (r*r-R.t*r+R.c) % p == 0]
        if R.D % p == 0:
            assert len(roots) == 1
            factors.append((p, 2*v))
        elif len(roots) == 2:
            factors.extend([(p,v),(p,v)])
        else:
            assert not roots
            factors.append((p*p,v))
    count, sigma = 1, 1
    for q, e in factors:
        count *= e+1
        sigma *= sum(q**j for j in range(e+1))
    return {"factor_norms_exponents": factors, "divisors": count, "sigma1": sigma}


def local_certificate():
    result = {"ramified": [], "inert": []}
    for D in [12,24,28]:
        alpha = 4/zeta_minus_one(D)
        t, power = 1, 6
        candidates = []
        while power <= alpha:
            # alpha=-3*eta(P)*2^t, ell=t+1.  If ell is even the
            # central character is ordinary, hence trivial (h=1).
            ell = t+1
            allowed = [1] if ell % 2 == 0 else [-1,1]
            if any(alpha == -power*eta for eta in allowed):
                candidates.append((ell, -1))
            t, power = t+1, 2*power
        assert not candidates
        result["ramified"].append({"D":D,"alpha":str(alpha),"candidates":candidates})
    for D, terminal_t in [(21,5),(69,2),(77,2)]:
        R = Ring(D)
        alpha = 4/zeta_minus_one(D)
        data = {m:divisor_data_integer(R,m) for m in [1,2,3]}
        assert data[2]["factor_norms_exponents"] == [(4,1)]
        # If 0<x,x'<4 and x=a+bw, then |b|sqrt(D)<4.
        # D>16 forces b=0, so the convolution at 4 has just 3 terms.
        assert D > 16
        decompositions = [(1,3),(2,2),(3,1)]
        convolution = sum(Q(data[x]["sigma1"]*data[y]["divisors"])*
                          Q(y,4)**terminal_t for x,y in decompositions)
        terminal = abs(alpha)*(convolution+4*Q(1,2)**terminal_t+
                               abs(alpha)*Q(1,4)**terminal_t)
        assert terminal < 15
        f3 = data[3]["sigma1"]
        # Expand the two formulas for G_(4), then cancel h_(2)^2:
        # G4 = h4+alpha*(h3+5*h2+f3)
        # G4-h4 = 2*alpha*h2+alpha^2-15*4^(ell-1).
        # Therefore 3*h2+h3 = alpha-f3-(15/alpha)*4^(ell-1).
        row = {"D":D,"alpha":str(alpha),"ideal_data":data,
               "decompositions_at_4":decompositions,"terminal_t":terminal_t,
               "normalized_upper_bound":str(terminal),"lower_bound":15,
               "f3":f3,"coefficient_relation":
               f"3*h2+h3 = {alpha-f3} - ({15/alpha})*4^(ell-1)"}
        if D in [69,77]:
            constant = alpha-f3-(15/alpha)*4
            ramanujan = 3*data[2]["divisors"]*2+data[3]["divisors"]*3
            assert constant == (-41 if D == 69 else -38)
            assert ramanujan == (21 if D == 69 else 18)
            assert abs(constant) > ramanujan
            row["weight_two_constant"] = str(constant)
            row["ramanujan_bound"] = ramanujan
            # The current manuscript excludes all ell>=2 directly: put
            # t=ell-2>=0. The required magnitude is >30*4^t, while
            # Ramanujan gives at most 12*2^t+9*3^t <=21*3^t.
            assert 30 > 21 and 4 > 3 and 3 >= 2
            assert alpha-f3 < 0
            row["all_source_weights_at_least_two_excluded"] = True
            row["infinite_weight_comparison"] = {
                "integer_parameter": "t=ell-2>=0",
                "required_magnitude_strictly_greater_than": "30*4^t",
                "Ramanujan_upper_bound": "12*2^t+9*3^t<=21*3^t",
                "base_and_ratio_checks": ["30>21", "4>3", "3>=2"],
            }
        else:
            assert R.norm((1,1)) == -3
            row["prime_above_3_generator"] = [1,1]
            row["generator_norm"] = -3
            row["resolution"] = "Pure coefficient proof in notes/d21_local_free.tex"
        result["inert"].append(row)
    return result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--local-only", action="store_true")
    args = parser.parse_args()
    result = {} if args.local_only else {"reduction":reduction_certificate()}
    result["local"] = local_certificate()
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(f"PASS: exact quadratic certificate written to {args.output}")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
