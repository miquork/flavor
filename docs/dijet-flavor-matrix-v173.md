# Dijet v173 four-flavor matrix extraction

## Inputs and decoded schema

The first-round extraction uses:

```text
data/jmenano_data_cmb_2026BD_JME_v173.root
data/jmenano_data_out_2026BD_JME_v173.root
data/jmenano_mc_cmb_Summer24MG_JMENANO_v173_v2.root
data/jmenano_mc_out_Summer24MG_JMENANO_v173_v2.root
```

The matching public producer revision is NestorMancilla/dijet commit
`c6be7f9` ("v173 Flavour QvG"). The implemented pair code is

```text
combined_id = 10 * tag_id + probe_id
```

not `tag_id + 10*probe_id`. Both reconstructed and truth pairs use the same
ordering. The analysis follows the code and verifies the ordering through the
tag/probe-swapped response symmetry.

The reconstructed tag is a sequential classifier:

1. UParT AK4 B > 0.4648 gives b;
2. otherwise UParT AK4 CvL > 0.650 gives c;
3. otherwise ParticleNet QvG >= 0.45 gives uds;
4. otherwise a valid QvG score gives g;
5. invalid scores give undefined.

The analysis groups truth u/d/s into uds and merges truth undefined/unmatched
with g. Reconstructed undefined is also merged with g, although it is empty in
these files. The producer's auxiliary `doQvsG_Eff` block labels low QvG as the
quark pass and high QvG as the gluon pass, opposite to `getFlavorID`; those
auxiliary histograms must not be used as validation of the matrix labels until
that convention is fixed.

Both leading jets act as tag and probe in turn. `isdijet` requires the tag jet
to be in the barrel and the matrix block separately requires the probe jet to
be in the barrel, so both jets satisfy |eta| < 1.3.

## File and profile validation

`DijetFlavorMatrix.C` performs three independent checks before inference.

- For data it compares every populated y-category vector in every pT bin of
  all 16 `p2{m0,m2,mn,mu}{ab,ad,tc,pf}_flavormatrix` profiles with every HLT
  source directory. All 488 populated slices have at least one exact source.
  The `ab`/`ad` profiles follow ZeroBias and DiPFJetAve triggers; `tc`/`pf`
  follow ZeroBias and PFJet triggers. The 155--180 GeV `ab` slice is identical
  in `HLT_DiPFJetAve80` and `HLT_PFJet80`, so the content is correct but the
  source is intrinsically ambiguous from the file alone.
- Every populated MC pT slice in all 2D and 3D variants is an exact copy of the
  corresponding `HLT_MC` slice (1024/1024 tested profile slices).
- Summing the `TProfile3D` bin entries and weighted contents over the truth
  axis reproduces the matching `TProfile2D`. Across all variants the largest
  relative count difference is 3.5e-14 and the largest absolute profile-mean
  difference is 3.9e-14.

The combined files have inflated global `TProfile::GetEntries()` values. The
bin-local `GetBinEntries()` values, effective entries, and profile means are
internally consistent and are what the analysis uses. The public v173
`DijetHistosCombine.C` still does not contain the `GluonJets/FlavorMatrix`
copy logic that produced these files. The empirical merge is sound, but that
missing code/provenance should be committed and the global entries bookkeeping
fixed.

The detailed checks are written to:

```text
results/dijetFlavorMatrix/trigger_merge_audit.tsv
results/dijetFlavorMatrix/mc_combine_audit.tsv
results/dijetFlavorMatrix/profile3d_to_2d_closure.tsv
```

## Tagging inference

For each pT range, data provide a 16-bin reconstructed pair distribution
`d(R)`. MC provide the full 16-by-16 reconstructed/true pair distribution
`M(R,F)`. The truth-pair fractions are fixed to their MC values, as requested.
The inferred data joint distribution is the minimum-KL projection closest to
MC with the observed data reconstructed marginal and fixed MC truth marginal:

```text
W_data = argmin_W KL(W || M)
subject to sum_F W(R,F) = d(R)
           sum_R W(R,F) = pi_MC(F).
```

Iterative proportional fitting solves this system to a maximum marginal error
below 1e-10. The pair result is then marginalized over the tag and probe
orientations to obtain per-jet efficiencies `epsilon(reco|truth)` and purities
`P(truth|reco)`. This preserves pair correlations while returning the requested
single-jet quantities.

For the inclusive 60--2000 GeV diagnostic, the diagonal values are:

| Flavor | MC efficiency | Inferred data efficiency | Data/MC | MC purity | Inferred data purity | Data/MC |
|---|---:|---:|---:|---:|---:|---:|
| uds | 0.694 | 0.796 | 1.147 | 0.504 | 0.433 | 0.860 |
| g+unmatched | 0.774 | 0.635 | 0.820 | 0.863 | 0.895 | 1.038 |
| c | 0.480 | 0.554 | 1.154 | 0.542 | 0.485 | 0.895 |
| b | 0.761 | 0.851 | 1.119 | 0.811 | 0.627 | 0.773 |

These are conditional central values, not calibrated tagger scale factors. The
large off-diagonal changes show why external tagger constraints and uncertainty
propagation are required before physics use.

The `ab` and `ad` tagging results agree at the per-mille level. Using the
single-jet `tc`/`pf` coordinate instead of the bisector average changes the
inclusive diagonal efficiency ratios by about +2.3% (uds), -3.2% (g), +0.3%
(c), and +1.0% (b). `tc` and `pf` become identical after both orientations are
marginalized, as required by construction. Their pair-resolved distributions
remain the appropriate place to study orientation-specific migrations.

## Relative response fit

HDM is constructed from the merged component means, never event by event:

```text
R_HDM = (m0 - mn - mu) / (1 - mn/Rn - mu/Ru),
Rn = 1.00, Ru = 0.92.
```

In the two-jet bisector limit, an HDM balance `h` maps to a probe/tag response
ratio

```text
q = (h + 1) / (3 - h).
```

For each reconstructed/true pair cell the MC value is retained. The data model
multiplies its ratio by `k(probe truth)/k(tag truth)`, with the same `k_f` in
every reconstructed category of a given true flavor. The predicted data
profile in each reconstructed pair is the inferred-purity-weighted sum of the
transformed cell responses. A nonlinear least-squares fit uses all populated
reconstructed pairs and both orientations. Dijet data contain no absolute
reference, so `k_uds = 1` fixes the gauge and g/c/b are relative to uds.

Representative conditional results are:

| pT range (GeV) | g/uds | c/uds | b/uds | chi2/ndf |
|---|---:|---:|---:|---:|
| 85--125 | 0.9685 +/- 0.0020 | 0.9801 +/- 0.0035 | 0.9808 +/- 0.0034 | 14.8/13 |
| 250--350 | 0.97785 +/- 0.00045 | 0.98395 +/- 0.00086 | 0.98500 +/- 0.00088 | 7.7/13 |
| 500--600 | 0.98043 +/- 0.00039 | 0.98487 +/- 0.00087 | 0.98573 +/- 0.00087 | 8.5/13 |
| 800--1000 | 0.98375 +/- 0.00045 | 0.99648 +/- 0.00120 | 0.99085 +/- 0.00113 | 12.2/13 |
| 1200--1500 | 0.9903 +/- 0.0014 | 1.0082 +/- 0.0035 | 0.9975 +/- 0.0034 | 7.8/13 |
| 1500--1800 | 0.9879 +/- 0.0032 | 1.0211 +/- 0.0074 | 0.9952 +/- 0.0067 | 2.7/9 |

The uncertainties are conditional profile-statistical errors with a 0.02%
numerical/statistical floor per reconstructed-pair row. They exclude tagging,
truth-fraction, response-template, trigger, JEC, HDM-component covariance, FSR,
and model uncertainties. The nominal fit is therefore a proof of principle,
not a 0.1% precision claim. Results above 1.8 TeV currently fail the response
rank/row requirement even though the tagging IPF still converges.

## Outputs and missing improvements

The macro writes a ROOT result, full pT-binned tagging and response tables, and
PDF/PNG plots under the gitignored `results/dijetFlavorMatrix/` and
`plots/dijetFlavorMatrix/` directories. The response-cell plots show MC and
inferred data separately in 350--500 GeV; the main graph shows the fitted
single-flavor data/MC ratios.

Before these results can enter the joint Z/gamma/dijet fit or a payload:

1. store explicit TH2D/TH3D counts (sumw, sumw2, and preferably unweighted
   entries) instead of relying on profile internals;
2. commit the FlavorMatrix combination logic and embed producer/combine commit
   hashes and trigger provenance in the ROOT files;
3. correct the reversed q/g labels in the auxiliary efficiency block and
   document the tagger model/WP and official SF/mistag payloads;
4. propagate profile-component covariances and MC response-template statistics;
5. replace fixed truth fractions with constrained nuisance parameters and add
   calibrated tagging-efficiency/mistag constraints;
6. add event-level bootstrap/jackknife covariance across the two orientations,
   reconstructed pairs, and pT coordinates;
7. test MC-as-data closure and injected tag/response shifts;
8. study generator choice, FSR/ISR, unmatched treatment, and truth definition;
9. keep `ab` as the central response coordinate and use `ad/tc/pf` as migration
   and method systematics rather than applying the bisector response equation
   indiscriminately to their different projection axes.
