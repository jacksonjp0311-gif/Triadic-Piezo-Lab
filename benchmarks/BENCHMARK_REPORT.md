# Benchmark Evidence Report — v8.0.0

Generated: 2026-07-14T21:45:38.553764+00:00

## Release conclusion

The benchmark lineage supports a bounded computational architecture: preserve rich local physics, compress only exposed coupling, certify relative to the observer and operating envelope, and retain temporal memory when cycle-to-cycle information persists. Negative results remain part of the evidence chain.

## Current evidence

| Family | Result | Interpretation |
|---|---:|---|
| Quasi-static pulse control | correlation `1.000000` | Scalar load scaling does not create a new spatial ranking |
| Triadic multigate | `29/29` passed; median reduction `88.54%` | Transfer, phase and extremal mismatch must all pass |
| Axisymmetric 2D | `15/15` passed; median exposed `134` of `1756` | Stable object is a frequency-closed dual manifold |
| Full 3D | `10/10` passed; median reduction `78.84%` | Volume, interface and global observables need separate authority |
| Observer manifold | `23` sentinels; final acceptance `100.00%` | Model validity is observer-conditioned |
| 4D temporal | Floquet max error `0.1070%`; hysteresis ratio `0.99555` | Memory and branch history can become state variables |

## Reproduction status

The portable 4D suite was rerun during package assembly. Its output is stored under `benchmarks/results/evolution6_4d/reproduced/`, with the console summary in `reproduction_log.txt`.

## Boundaries

- The 2D and 3D reduced models are certified relative to their full-order discretizations; that is not equivalent to continuum mesh convergence.
- Observer closure is bounded to the tested observer and operating envelopes.
- The optical particle layer in the HTML is qualitative.
- No result in this archive constitutes experimental device calibration or proof of a universal physical law.
