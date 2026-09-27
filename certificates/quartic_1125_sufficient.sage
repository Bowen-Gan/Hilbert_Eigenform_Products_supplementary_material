"""Free exact certificates sufficient for the D=1125 endpoint.

Run with Sage, or Python with passagemath-standard:
  python free_replacements/quartic_1125_sufficient.sage

Outputs S_2(1)=0 and dim S_5(epsilon)>=5, using the mathematical argument
in quartic_1125_sufficient_proof.md.  This does NOT assert the manuscript's
exact dimensions S_3=2 or S_5=6.  It proves the two facts needed to exclude
the remaining source weights, with no Magma call and no numerical zeta
reconstruction.
"""

from sage.all import *
from sage.version import version as SAGE_VERSION
from pathlib import Path
import argparse
import json

proof.all(True)


def elt(a):
    return [str(c) for c in a.list()]


def quat(q):
    return [elt(c) for c in q.coefficient_tuple()]


def main(output_path=None):
    x=polygen(QQ)
    F=NumberField(x**4-x**3-4*x**2+4*x+1,"a")
    a=F.gen()
    assert F.discriminant()==1125 and F.signature()==(4,0)
    assert F.class_group(proof=True).order()==1
    assert F.narrow_class_group().order()==2
    units=F.units(proof=True)
    assert all(u.norm()==1 for u in units)

    # Universal reduced-unit-group bound in degree four: a projective
    # element of order m gives Q(zeta_m+zeta_m^-1) subset F, hence phi(m)<=8.
    # The elementary inequality phi(m)^2>=m/2 reduces the list to m<=128.
    possible_orders=[m for m in range(1,129) if euler_phi(m)<=8]
    assert max(possible_orders)==30
    max_reduced_unit_order=max(2*max(possible_orders),12,24,60)
    assert max_reduced_unit_order==60

    # Mass upper bound using just the Euler factors at 2 and 3.
    # Elsewhere use zeta_F(2) <= zeta(2)^4, then cancel pi^8 exactly.
    ratio=QQ(1)
    local=[]
    for p in (2,3):
        factors=list(F.ideal(p).factor())
        degrees=[ZZ(P.residue_class_degree()) for P,e in factors]
        actual=prod(1/(1-QQ(p)**(-2*f)) for f in degrees)
        ratio*=actual*(1-QQ(p)**(-2))**4
        local.append({"p":p,"factors":[[int(e),int(P.residue_class_degree())]
                                           for P,e in factors]})
    assert local[0]["factors"]==[[1,4]]
    assert local[1]["factors"]==[[2,2]]
    mass_coefficient=ratio/(2**7*6**4)
    assert mass_coefficient==QQ(1)/826200
    # mass < 1125*sqrt(1125)/826200 < 1/20.
    mass_upper_squared=mass_coefficient**2*1125**3
    assert mass_upper_squared < QQ(1)/400
    assert 60**2*mass_upper_squared < 9
    # Norm surjectivity onto Cl^+(F) gives >=2 ideal classes; the mass
    # inequality and w<=60 give <3.  See proof note for these theorems.

    Q=QuaternionAlgebra(F,-1,-1)
    i,j,k=Q.gens()
    tau=a**3-3*a+1
    assert tau**2-tau-1==0
    z=(tau+(tau-1)*i+j)/2
    basis=[Q(1),i,z,i*z]
    basis_matrix=matrix(F,[q.coefficient_tuple() for q in basis])
    multiplication_table=[]
    for u in basis:
        row=[]
        for v in basis:
            coefficients=basis_matrix.transpose().solve_right(vector(F,(u*v).coefficient_tuple()))
            assert all(c.is_integral() for c in coefficients)
            row.append([elt(c) for c in coefficients])
        multiplication_table.append(row)
    gram=matrix(F,[[(u*v).reduced_trace() for v in basis] for u in basis])
    assert gram.det()==-1
    # Closure + unit trace discriminant proves this is a maximal order,
    # unramified at every finite prime.

    group={Q(1)}
    frontier=[Q(1)]
    while frontier:
        u=frontier.pop()
        for generator in (i,z):
            v=u*generator
            if v not in group:
                group.add(v)
                frontier.append(v)
            assert len(group)<=120
    assert len(group)==120
    assert all(g.reduced_norm()==1 for g in group)
    assert all(g in group for g in (-Q(1),-i,-z))
    # The universal projective-order bound above proves the full reduced
    # unit group is this 120-element group's quotient by {+1,-1}.

    elements=sorted(group,key=str)
    character_rows=[]
    character_sum=ZZ(0)
    for g in elements:
        t=g.reduced_trace()
        character=ZZ((t**3-2*t).norm())
        character_sum+=character
        character_rows.append({"element":quat(g),"trace":elt(t),
                               "weight_five_character":int(character)})
    invariant_dimension=QQ(character_sum)/120
    assert invariant_dimension==5

    result={"schema":"D1125-free-sufficient-endpoint-v1","sage_version":SAGE_VERSION,
        "polynomial":[1,4,-4,-1,1],"discriminant":1125,
        "ordinary_class_number":1,"narrow_class_number":2,
        "fundamental_units":[elt(u) for u in units],"fundamental_unit_norms":[1]*len(units),
        "possible_projective_element_orders":possible_orders,
        "universal_reduced_unit_group_bound":60,
        "local_euler_factors":local,"mass_upper_coefficient":str(mass_coefficient),
        "mass_upper_squared":str(mass_upper_squared),
        "mass_upper_less_than_one_twentieth":True,
        "quaternion_ideal_class_number":2,
        "tau":elt(tau),"maximal_order_basis":[quat(q) for q in basis],
        "multiplication_table":multiplication_table,
        "trace_gram":[[elt(c) for c in row] for row in gram.rows()],
        "trace_gram_determinant":"-1","norm_one_subgroup_order":120,
        "full_reduced_unit_group_order":60,
        "weight_five_character_sum":int(character_sum),
        "principal_ideal_class_weight_five_invariants":5,
        "weight_five_character_rows":character_rows,
        "conclusions":{"dim_S2_trivial_character":0,"dim_S5_totally_odd_character_at_least":5},
        "exact_S3_or_S5_dimension_claimed":False,
        "theorem_dependencies":["finite subgroups of SO(3)","Eichler mass formula",
            "surjectivity of the reduced-norm map to the narrow class group",
            "definite quaternionic modular-form decomposition",
            "Jacquet-Langlands correspondence"]}
    out=Path(output_path) if output_path else Path(__file__).with_suffix(".json")
    out.write_text(json.dumps(result,indent=2))
    print("D=1125: maximal quaternion order verified, binary icosahedral subgroup order 120.")
    print("Quaternion ideal classes =2; weight-two cusp space =0.")
    print("One ideal class contributes5 weight-five invariants; dim S5(epsilon)>=5.")
    print("Output:",out)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",help="write the JSON certificate to this path")
    args=parser.parse_args()
    main(args.output)
