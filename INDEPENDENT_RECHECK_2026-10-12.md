# Independent re-verification — 12 October 2026

Every quantitative claim in the manuscript was re-derived from the paper's own text by
scripts written from scratch for this pass (`recheck_A.py`, `recheck_B.py`, `recheck_C.py`,
`recheck_D.py`, `recheck_E.py`, `attain_l3_recheck.py`, `focus1.py`), independently of the
four suites in the replication module (`audit_paper1.py`, `verify_unified.py`,
`attain_check.py`, `coop_create.py`). Roughly 300 checks; nothing contradicted the paper.

## What was re-derived

**Section 3 (cooperative).** Example 3.8: superadditivity, core membership, which constraints
bind, m = 4 and ℓ = 3, every 3×3 incidence minor ±1, the stationarity cubic
4s³−54s²+202s−168 = 2(2s³−27s²+101s−84), its irreducibility, the unique root in (3,7) to
7 digits, degree 3 for all four coordinates. Proposition 3.4 (3-player bargaining set):
the closed form agreed with a direct objection/counter-objection search on 486 random
imputations. Proposition 3.6: on a generic cell the τ-value's numerator is genuinely
quadratic in v over an affine denominator (found a cell with numerator degree 2, denominator
degree 1), confirming the corrected statement; on special cells it degenerates to affine,
which the "has a quadratic numerator" phrasing covers.

**Section 4 (market games).** Assignment Shapley (5/6, 7/3, 29/6); bankruptcy
v = (0,0,0,10,30,60,100) and Shapley (20,35,45); √(xy) exchange values and Shapley
(19/6−√2/2, 19/6−√2/2, 2/3+√2); the cube-root family — radicands {36,112,140,324} with
F₃-rank 4 giving 3⁴ = 81, and {2,5,24,28,104,116,279} on exactly the primes
{2,3,5,7,13,29,31} with rank 7 giving 3⁷ = 2187; Cobb–Douglas price 40/27; CES numerator
−19q²+7q+25 and p* = (7+√1949)²/1444.

**Section 5 (equilibrium).** Example 5.2: (0,1,1) is a strict pure equilibrium with payoffs
(2,3,1) against deviations 1, 0, 0; the three indifference polynomials match; the unique
interior solution and its minimal polynomials 9X²−16X+6, 6X²−16X+9, 18X²−16X+3.

**Section 6 (bargaining).** The facet formula verified symbolically. Bounded-chamber counts
computed by exact sign-vector enumeration plus LP feasibility and recession-cone tests for
eight arrangements: in every case b(𝒜) ≤ C(m−1,k), with equality in general position.
Example 6.7: ν = ((√13−2)/9, (5−√13)/12, (7+√13)/36), its three minimal polynomials, both
multipliers positive (4−√13 and √13−1), ρ = 1 interior to (18/19, 4/3).

**Theorem 6.11 (the headline).** The Fisher system for five agents with columns (1,i,i²) and
capacities (1,2,11/2) was re-eliminated by Gröbner basis; the paper's sextic f divides the
eliminant exactly and is irreducible; f has six real roots, exactly one positive, in (0,1);
y = Y(x) is pinned linearly with Y of degree 5 mod f; prices (0.852703, 0.713858, 0.494469)
and utilities (0.485194, 0.234836, 0.134327, 0.086061, 0.059582) to 6 decimals; all three
resources clear exactly and p·c = 5. Every price and every utility was reduced in ℚ[x]/(f)
and is irrational there, so by the maximality of the point stabiliser S₅ in S₆ each has
degree exactly six. **Gal(f) independently recomputed: order 720, nonsolvable.**
**Gal(f₅) independently recomputed: order 120, nonsolvable.**

**Theorem 6.8 (attainment).** Re-derived at ℓ = 2 for m = 3,4,5,6 by exact elimination, and
at ℓ = 3 for m = 4 (degree 3) and m = 5 (degree 6) by 180-digit solution plus PSLQ integer
relation detection — in each case the recovered minimal polynomial is irreducible over ℚ,
has degree exactly C(m−1,ℓ−1), and all its roots are real.

**Section 7 (two resources).** For three technologies (m = 3, 4, 6) all of the following were
verified symbolically: R₁ < H < A < R_m; g′ = WL/S²−1 > 0; g(0) = H and g(∞) = A;
E(t,g(t)) ≡ 0; deg_T E = m−1 with leading coefficient W(Z−A); E irreducible in ℚ(Z)[T];
gcd(N,D) = 1; the parametrisation's feasibility, clearing and stationarity; λ+μρ = W;
t = λ/μ; x_i = w_i/(S(t)(t+r_i)); C(t) = S/L−t; R₁ < C < ρ < A; dx_i/dρ and its sign rule;
λ′ > 0 and μ′ < 0; μ_t = −L and λ_t = S−tL = LC; C′ = 2(SM₃−L²)/L² > 0; C(0) = C₀ and
C(∞) = A; and the matching-capacity derivative w_i/(W(t+r_i)). Six-type technology:
H = 120/49, A = 7/2, N = D′, f₅ = −5E(T,16/5), f₅(6) = −31728, f₅(7) = 60408, five real roots
one apiece in (−6,−5), (−5,−4), (−3,−2), (−2,−1), (6,7); every coordinate and both prices of
degree 5; t recoverable from either price. Modular certificates recomputed: f mod 11 → (6),
mod 7 → (1,5), mod 17 → (1,2,3), squarefree and unit leading coefficient at each;
f₅ mod 19 → (5), mod 13 → (2,3). Congruence family: 247 = 13·19 and the stated capacities
give irreducible quintics with group of order 120; the density construction converges at the
stated rate. Appendix A: gcd(E,D) = T+R_k at matching capacities and 1 otherwise, the sign
pattern sgn E(−R_j,ρ) = sgn(ρ−R_j)(−1)^{j−1}, and m−1 distinct real roots with exactly one
positive.

**Section 8 (expansion).** Proposition 8.9's ∂u_i/∂c_k = (w_i/q_i²)(G⁻¹a_i)_k and
∂p/∂c_k = −G⁻¹e_k checked against 50-digit central differences for all three resources
(agreement to ~1e-20); the usage identities Σ_i a_{ji}∂u_i/∂c_k = δ_{jk} hold to 1e-61; each
expansion has both gainers and losers. Corollary 8.10: g₁₂/g₁₁ = C(t) symbolically.
Remark 8.11: the capacities at which all three prices equal one are exactly
(137/84, 271/168, 295/168), and each of a₁, a₂, a₃, a₂/a₁, a₃/a₁, a₃/a₂ fails to order the
gainers for at least one expanded resource while each aⱼ orders it for one.

**Table 1 and Section 9.** All 24 table entries reproduced exactly, with feasibility and
exact clearing at each; K = 112/65; the middle user's reversal (170/503 < 20/59 while
1132/503 > 124/59); C(t₀) < K < C(t₁); t* = 2+3√2, ρ* = 13/7+3√2/14 ≈ 2.16019,
x₂(ρ*) = 9/14−3√2/14, x₂′(2) = 1/12. Section 9.1: t⁰ = 5/13, x⁰ = (11/20, 9/40, 9/40),
the quadratic 3t²−8t−29 with root (4+√103)/3, H(ε) ≤ 9/5 and A(ε) ≥ 2 on [0,1].
Section 9.2: t₀ = 1 at ρ₀ = 617/223, x⁰_i = 140/(223(i+1)), p⁰ = 223/140;
ρ₁ = ρ₀+1/5513006 satisfies the congruence with both capacities in (120/49, 7/2);
E(T,ρ₁) is an irreducible quintic of Galois group order 120; max_i|x¹_i−x⁰_i| = 4.84e-8
(< 5e-8 as stated) and the largest price change 9.32e-7 (< 1e-6 as stated).

## Two proof-level repairs made in this pass

1. **Theorem 6.14, second application of IIA.** The proof concluded φ(U,0) = φ(Û,0) from
   Û ⊆ U, but IIA in the stated form needs φ(U,0) ∈ Û. That holds because Û₀ = U₀ — a point
   of Û with u ≥ 0 is dominated by a point of conv(U₀) = U₀ and so lies in U by
   comprehensiveness — and the proof now says so.
2. **Theorem 6.8, the box.** The birational map φ_k(t) = a_k + (b_k−a_k)t/(1+t) must be
   defined over ℚ, so the box needs rational endpoints; the proof now says so and notes that
   such a box exists because 𝒞 is nonempty and open.

## Where the recheck stopped short of the paper

The PSLQ route could not recover a degree-10 minimal polynomial at ℓ = 3, m = 6 within 180
digits — the coefficients exceed the detection bound. The paper claims attainment at ℓ = 3
only for m ∈ {4,5} (Appendix D), both of which were confirmed, so this is a limit of the
check and not of the result.

## Legacy suites

`audit_paper1.py` 48, `verify_unified.py` 46, `attain_check.py` 8, `coop_create.py` 16 —
**118 checks, 0 failures**, re-run against the retitled source.
