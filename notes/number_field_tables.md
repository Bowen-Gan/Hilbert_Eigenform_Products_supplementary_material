# Complete finite number-field inputs

The corrected Bordeaux/PARI tables provide the complete external number-field
lists used by the weight-two finite reduction. The primary directory and format
specification are

- https://pari.math.u-bordeaux.fr/pub/pari/packages/nftables/
- https://pari.math.u-bordeaux.fr/pub/pari/packages/nftables/README.txt

Karim Belabas prepared the corrected PARI version in 2007; the individual files
in this distribution date from September 2008. The unmodified decompressed GP
texts are bundled in `data/number_field_tables`, with deterministic gzip wrappers
(`mtime=0`). `manifest.json` records source URLs, retrieval date, published table
ranges, full row counts, source-text and gzip SHA256 hashes, and every exact
source-row duplication. The package makes no byte-identity claim about the
upstream gzip header; the underlying downloaded GP text is unchanged.

| Degree | Source | Published complete range | Full source rows | Inclusive manuscript cutoff | Fields used |
|---|---|---:|---:|---:|---:|
| 3 | T33.gp.gz | D < 2000000 | 112444 | 4218 | 143 |
| 4 | T44.gp.gz | D < 1000000 | 13073 | 68133 | 552 |
| 5 | T55.gp.gz | D < 20000000 | 22740 | 209508 | 37 |
| 6 | T66.gp.gz | D < 10000000 | 398 | 2429064 | 40 |
| 7 | T77.gp.gz | D < 150000000 | 154 | 28162991 | 4 |

There are five identical duplicated quintic rows in the full upstream T55 text,
at discriminants 9262117, 13072837, 14731145, 17946025 and 18371721. These are
outside the manuscript cutoff. They are preserved and explicitly checked;
no duplicate occurs among the 37 quintic fields used here. Distinct defining
polynomials with the same discriminant are always retained as separate source
rows. Completeness and identification up to isomorphism are published table
inputs, rather than consequences of a count assertion in this package.

`number_field_table_inputs.py` verifies and parses every source row. GP's vector
of coefficients is highest-degree first; the output manifest uses the existing
constant-first JSON convention. `e2_candidate_enumeration.sage` recomputes maximal
orders, signatures, proof-enabled ordinary class numbers, dyadic factors and
necessary small-prime factors, and the certified finite Euler intervals. Its
ordinary/narrow class groups are Sage arithmetic results, not inferred only from
the table. The complete lower-degree input is 772 fields, and 21 require the
separate lattice or exceptional arithmetic certificates.

The previously omitted quartic is

`x^4 - 19*x^2 - 24*x + 16`, field discriminant `65808`.

Its power-basis discriminant is `4211712 = 64*65808`, so its index is eight.
The actual maximal-order factorization is `(2) = P1^2 * P2^2`, with both prime
norms equal to two. A repeated factor of the defining polynomial modulo two
would not establish this, because two divides the power-basis index. The new
Sage output records the maximal-order bases and factor ideals. The independent
Python verifier checks their ring arithmetic and also verifies ramification
without Dedekind's power-basis test: Frobenius on the four-dimensional algebra
`O_F/2O_F` has rank two, so this algebra is not reduced. Under the unchanged
screening order the field is already excluded by the class-number-one bound
`D < 18070.048`; its ramification supplies the additional exclusion stated in
the manuscript.

`verify_field_enumeration.py` uses standard-library rational arithmetic. It
checks source/output coverage, order discriminants, ring closure, power-basis
indices, norms and ideal closure, dyadic ramification, the exact Euler products,
and integer intersections using an independent rational enclosure of pi from
Machin's formula. It states its external inputs explicitly and does not claim
to independently enumerate number fields or recompute class groups.

The square-ramified quartic certificate reads all 17 source fields through 4446;
it finds discriminants 1600, 2000, 2624, 3600 and 4400, all with ordinary class
number one. The high-degree certificate uses all four septic source rows below
the degree-seven bound. Its three octic inputs and their completeness continue
to come from Voight's published Table 4, as stated in the manuscript and script.
