# Accretion/evaporation asymmetry of the adiabatic correction δφ_V

## Physical idea
At first order δφ ∝ ṁ ⇒ ṁ→−ṁ only flips the sign (trivial). The real asymmetry:
- **Accretion** (ṁ>0): **ingoing** Vaidya, clock = **advanced** time v=t+r∗.
- **Evaporation** (ṁ<0): **outgoing** Vaidya, clock = **retarded** time u=t−r∗
  (ds²=−f du²−2du dr).

The frozen reduction ∂_m F is IDENTICAL (it depends only on S); only the clock changes, and
v↔u flips the sign of the tortoise coordinate r∗.

## Decomposition (verified to 1e-15)
With v = E U₃ + r∗,  r∗ = (r−r₀) + 2m ln((r−2m)/(r₀−2m)):

    δφ/ṁ|accr = A∞ + B_hor        δφ/ṁ|evap = A∞ − B_hor
    A∞  = ∫ ∂_m F · E U₃ dr   (SYMMETRIC, polylog at INFINITY)
    B_hor = ∫ ∂_m F · r∗   dr   (ANTISYMMETRIC, tortoise)

**Net asymmetry** (at equal |ṁ|): δφ_accr − δφ_evap = 2|ṁ| B_hor,
carried ENTIRELY by the tortoise term; A∞ (the E U₃ sector) is common and CANCELS.

## Numbers (M=1, E=1.4, J=2.5)
- A∞ = 8.72756968   (polylog at ∞)
- B  = 6.11272429 = B_bulk(4.58262818, second kind) + B_hor-log(1.53009610, dilog D_k)
- accr = A+B = 14.84029397  (= clock v, direct match 3.6e-15)
- evap = A−B =  2.61484540  (= clock u, direct match 8.9e-16)
- asymmetry = 2B = 12.22544857
- (the three pieces 8.73/4.58/1.53 coincide with the split of the clock in vaidya_clock.py)

## Structure
B_hor by parts: 2m[A_m ln(r−2m)] − 2m ∫A_m/(r−2m) dr, and ∫A_m/(r−2m) = Σ c_k^M D_k
(HORIZON dilog). ⇒ the transcendental part of the asymmetry IS the dilog D_k already
derived. The accretion/evaporation asymmetry is the PHYSICAL face of the weight-2 horizon piece, as the
ergosphere penetration window (finite measure for one sign of J, zero for
the other) is the PHENOMENOLOGICAL face of the same horizon geometry.

Script: `vaidya_asymmetry.py`. Paper: eq:vaidya-asymmetry (main + PRD).
