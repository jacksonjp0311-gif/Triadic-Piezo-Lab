# Piezoelectric Adaptive FEM — Two-Cycle Evolution Findings

## Experiment

Two solver evolution cycles were run against 1,025-node reference models across five deterministic scenarios and four degree-of-freedom budgets. A separate 24-case randomized holdout tested whether the final governance signals generalized.

### Cycle 1

Added wave-resolution pressure, mesh-size coverage, neighbor smoothing, and bulk marking to the original pulse-integrated residual method.

### Cycle 2

Added an entropy-derived locality authority that continuously allocates refinement between a localized residual channel and a global coverage channel:

```text
p_e = eta_e^2 / sum_j eta_j^2
A_L = clamp((1 - exp(-sum_e p_e log p_e))/N_e - 0.08)/0.48, 0, 1)
S_e = A_L S_local,e + (1-A_L) S_coverage,e
```

## Main deterministic results

| Scenario, 49 nodes | Uniform QoI error | Legacy pulse | Cycle 1 | Cycle 2 | Result |
|---|---:|---:|---:|---:|---|
| broadband_standard_8_280kHz | 0.4226% | 4.2067% | 2.4668% | 1.2525% | uniform required |
| dual_interface_5_120kHz | 0.5506% | 1.3323% | 1.8741% | 1.8310% | uniform required |
| localized_shifted_4_80kHz | 0.4193% | 5.1995% | 5.5154% | 0.3104% | adaptive win |
| localized_standard_3_60kHz | 0.0969% | 0.0446% | 1.1652% | 0.0472% | adaptive win |
| smooth_global_12_300kHz | 0.1317% | 0.1937% | 0.1937% | 0.1317% | near tie |

## Aggregate changes

- Cycle 2 reduced geometric-mean QoI error by **55.5%** relative to Cycle 1.
- Cycle 2 reduced geometric-mean QoI error by **47.2%** relative to the legacy pulse method.
- Mean mesh irregularity ratio fell by **54.8%**, from 10.50 to 4.75.
- Uniform refinement still had the best aggregate geometric-mean QoI error: 0.4702% versus 0.5616% for Cycle 2.

## Randomized holdout

- Cycle 2 beat uniform refinement in **13/24** cases.
- Median governed/uniform error ratio: **1.000**.
- Geometric-mean governed/uniform error ratio: **1.332**; rare large losses dominate this value.
- Locality alone was not a reliable gate. The best simple locality-times-unresolvedness threshold classified only **62.5%** of holdout cases correctly.

## Emergent conclusion

The benchmark does not support unconditional adaptive sparsening. It supports a three-gate architecture:

```text
G_dynamic   = evidence that the interrogation contains non-proportional dynamic information
G_locality  = concentration of pulse-integrated residual energy
G_necessity = evidence that the current representation is materially unresolved

CompressionAuthority = G_dynamic * G_locality * G_necessity
```

In the simulator, a final equal-DOF reference-validation gate is added: the adaptive map is displayed only when its measured error is at least 3% lower than the uniform candidate. Otherwise the tool renders the uniform safety map.

## Scientific status

The two-cycle method is a validated reduced-order computational experiment. It is not yet a general theorem or experimental piezoelectric result. The strongest current finding is that **residual concentration determines where refinement wants to go, but unresolvedness and quantity-of-interest relevance determine whether it should be allowed to go there**.
