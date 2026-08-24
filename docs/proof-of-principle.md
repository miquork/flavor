# Four-flavor proof-of-principle checkpoint

## Scope

This checkpoint answers two immediate questions:

1. Can the existing Z+jet and photon+jet tagged HDM inputs support a single
   simultaneous uds/g/c/b response fit?
2. Are the currently unimported dijet q/g-tag profiles usable as reco-level
   diagnostics before tag-purity unfolding?

The answer to both is yes, with important qualifications below.  None of these
outputs should yet be used as a CMSSW correction or FactorizedJetCorrector
payload.

## Joint Z+jet and photon+jet observable

For channel `s` and reconstructed tag `t`, the fitted observable is

```text
D_st = (HDM_data / HDM_MC)_st / (HDM_data / HDM_MC)_s,inclusive .
```

Using a tag divided by the inclusive balance in the same channel is essential.
At about 100 GeV the absolute inclusive Z/photon ratios differ by 0.80%, 0.84%,
and 1.83% in 2024, 2025, and 2026.  An absolute joint fit would incorrectly
absorb these channel-wide modes into flavor response.

The response share of truth component `f` is built from the MC count and its MC
response,

```text
A_stf = N_stf * R_stf / sum_h(N_sth * R_sth).
```

The unmatched component is kept in the response shares and fixed to unit
residual response.  Its gamma truth category is not considered validated yet;
the producer reads `GenJet_partonFlavour[iGenJet]` even when `iGenJet == -1`.
Because its share is tiny in the fitted window, it is retained only as a fixed
closure component for this checkpoint.

The fit linearizes the log double ratio.  Every matched flavor uses the same
generic form

```text
delta_f(pT) = a_f + b_f * ((pT / 100 GeV)^(-0.3) - 1).
```

The coefficients are independent between years.  Relative observables have
common intercept and slope null modes.  Both are fixed by requiring the
equal-channel 100 GeV inclusive reference mixture to have zero weighted
intercept and slope.  This is a convention, not a measured absolute response.
The nominal reference weights `(uds,g,c,b)` are approximately

| Year | uds | g | c | b |
|---|---:|---:|---:|---:|
| 2024 | 0.59231 | 0.20844 | 0.13995 | 0.05930 |
| 2025 | 0.59253 | 0.20830 | 0.13993 | 0.05923 |
| 2026 | 0.59242 | 0.20837 | 0.13997 | 0.05925 |

The central fit uses nominal MC shares.  `JetFlavor.C(year,true,0.2)` provides a
non-negative count-unfolded composition variant, but it must not be interpreted
as a corrected central value: without tag-SF, mistag, and background nuisances,
that procedure absorbs detector/category mismodeling into physics fractions.

## First fit result

The central window is 90--300 GeV.  Stored tag and inclusive errors are
propagated with the common inclusive denominator shared across tags in a native
bin.  A 0.2% diagonal floor is added solely to diagnose smooth-model closure,
and the displayed fit covariance is enlarged by `max(1, chi2/ndf)`.

At 100 GeV, relative to the declared reference mixture:

| Year | uds (%) | g (%) | c (%) | b (%) | chi2/ndf |
|---|---:|---:|---:|---:|---:|
| 2024 | -0.112 +/- 0.140 | -1.231 +/- 0.372 | +1.608 +/- 0.247 | +1.646 +/- 0.223 | 60.7/34 |
| 2025 | +0.091 +/- 0.140 | -1.207 +/- 0.374 | +0.598 +/- 0.250 | +1.920 +/- 0.231 | 60.1/34 |
| 2026 | +0.478 +/- 0.126 | -2.894 +/- 0.335 | +1.240 +/- 0.233 | +2.466 +/- 0.252 | 37.9/34 |

These are proof-of-principle values, not total uncertainties.  The dominant
2024/2025 tension comes from Z+jet; the per-channel chi2 values are preserved in
the result text files.  The curves are intentionally not extrapolated outside
90--300 GeV at this checkpoint.

## Dijet reco diagnostics

`DijetReco.C` reads:

```text
data: ../jecsys3/rootfiles/Prompt/Jet_v170/
      jmenano_data_cmb_<era>_JME_v170.root
MC:   ../jecsys3/rootfiles/Prompt/Jet_v171/
      jmenano_mc_out_Summer24MG_JMENANO_v171.root
```

from `GluonJets/tight` in data and `HLT_MC/GluonJets/tight` in MC.  It profiles
the bisector-average (`ab`) response for `qq`, `qg`, `gq`, and `gg` reco tag
pairs.  The current dijet definitions are a PNet QvG threshold of 0.45 after an
UParT heavy-flavor veto.  Scale factors are disabled in the producer.

The MPF0 reco contrast

```text
100 * [(gq - qg)_data - (gq - qg)_MC]
```

is already strongly year dependent:

| Year | 92.5 GeV | 195 GeV | 450 GeV | 900 GeV |
|---|---:|---:|---:|---:|
| 2024 | -0.586 +/- 0.189 | -0.620 +/- 0.065 | +0.048 +/- 0.034 | +0.050 +/- 0.042 |
| 2025 | +0.324 +/- 0.184 | +1.224 +/- 0.070 | +1.568 +/- 0.034 | +0.656 +/- 0.042 |
| 2026 | +0.786 +/- 0.421 | +1.198 +/- 0.148 | +2.129 +/- 0.070 | +1.799 +/- 0.086 |

MPF2 differs materially from MPF0, especially at low pT, while the neutral and
unclustered terms expose which recoil components drive the difference.  The
pair fractions also differ strongly between data and MC.  These observations
make the profiles useful constraints for a future likelihood, but they are not
yet quark/gluon response measurements.

## Tagger provenance audit

- Configured photon data (`Gam_w83`) and MC (`Gam_w73`) use legacy DeepFlavour:
  DeepFlavB > 0.7527; then average(CvB,CvL) > 0.3985; then DeepFlavQG >= 0.5
  for q and < 0.5 for g.  They are not the newer UParT/PNet definitions.
- The exact Z v115 producer was not identified publicly.  The normalized schema
  is known, and a non-production Z branch uses the same legacy sequence, but
  that is not sufficient provenance for a final result.
- Dijet v168/v170 uses the UParT heavy veto followed by PNet QvG 0.45.  Its
  `CvL <= 0.421` veto should be checked because 0.421 is otherwise documented as
  a CvB threshold.  The auxiliary efficiency diagnostic also reverses the q/g
  score convention, although the response-pair filling itself uses high score
  for q.
- Gamma data and MC come from different producer work points (w83 and w73), and
  the photon flavor code warns that Run 3 flavor histograms were not checked.

A regenerated input manifest must store producer repository/commit, NanoAOD or
JMENano campaign, jet collection, tagger branches and model, numerical working
points and boundaries, sequential priority, truth matching, and SF/mistag
payload identity/checksum.

## Required gates before a payload

1. Add a forward-folded joint likelihood for tag yields and responses.  Flavor
   abundances should be free per sample/bin; tag efficiencies and mistags need
   calibrated, correlated nuisance parameters.
2. Add Z+jet/top and photon+jet/EM-QCD components with constrained
   normalizations and shapes.
3. Bootstrap or jackknife event-level tag/inclusive and MPF/DB covariance.
   Scalar histogram errors cannot model the current overlap.
4. Propagate response templates, tag efficiencies, backgrounds, FSR/UE, pileup,
   HDM construction, and residual-JEC uncertainties with correlated toys or
   profiled nuisances.
5. Fix and validate gamma unmatched truth, align taggers across samples, and
   add an MC-as-data closure suite before unblinding updated data fits.
6. Add W->cs output for ud/s separation.  Until then the Run 2 numbers
   `R(c)=+0.6%`, `R(s)=-0.4%`, and `R(ud)=+0.2%` may be encoded only as a clearly
   labeled implementation prior/test, not as a Run 3 measurement.
7. Only after those gates, extend pT and eta coverage and fit MC truth-flavor JEC
   functions for the standalone JER-scaling application interface.
