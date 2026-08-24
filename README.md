# Jet flavor response

Standalone ROOT/C++ tools for Run 3 CMS AK4 PUPPI residual jet-flavor response
studies.  The intended final scope is ud, s, g, c, and b response from 15 GeV to
4.5 TeV and eventually |eta| < 5.2.  The current checkpoint is a barrel-only
(|eta| < 1.3) four-flavor proof of principle using uds, g, c, and b.

## Current checkpoint

- `JetFlavor.C` jointly fits Z+jet and photon+jet tag-to-inclusive HDM balance
  double ratios.  It uses one generic smooth function for every flavor, with an
  independent fit for 2024, 2025, and 2026.
- `DijetReco.C` reads the primitive `GluonJets/tight` profiles from the sibling
  `jecsys3` checkout and makes reco-level tag-pair balance, fraction, and
  gq-minus-qg contrast plots.  Dijet is not yet included in the flavor fit.
- `tdrstyle_mod22.C` supplies the CMS plot style.

The fit is a diagnostic, not a JEC payload.  It uses nominal MC response shares,
an explicit relative-response anchor, stored scalar response errors, and a 0.2%
diagnostic floor over 90--300 GeV.  Tagging scale factors, mistags, backgrounds,
event-level covariance, and MC-template uncertainties are not yet profiled.

## Run

ROOT 6.34 has been used for this checkpoint.  From this directory:

```bash
root -l -b -q 'JetFlavor.C(2024,false,0.2)'
root -l -b -q 'JetFlavor.C(2025,false,0.2)'
root -l -b -q 'JetFlavor.C(2026,false,0.2)'

root -l -b -q 'DijetReco.C(2024)'
root -l -b -q 'DijetReco.C(2025)'
root -l -b -q 'DijetReco.C(2026)'
```

`JetFlavor.C` expects private inputs named
`data/jecdata{2024,2025,2026}FLAVOR.root`.  `DijetReco.C` defaults to
`../jecsys3`; pass a second argument if that checkout is elsewhere.

Generated `plots/`, `results/`, ROOT inputs, and private payloads are gitignored.
The ROOT outputs contain fit covariance, input rows/points, and provenance; the
text and CSV outputs provide human-readable summaries.

See [docs/proof-of-principle.md](docs/proof-of-principle.md) for the model,
results, provenance audit, limitations, and next implementation gates.
