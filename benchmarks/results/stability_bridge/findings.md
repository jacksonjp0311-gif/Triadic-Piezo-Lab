# Triadic Stability Bridge — Evolution Findings

## Status

This evolution replaces unstable geometric sparsening with a nested component-reduction architecture:

1. **Physical manifold** — retain the resolved finite-element geometry locally.
2. **Interface manifold** — hold the port coordinates invariant.
3. **Coherence gate** — admit local fixed-interface modes until all declared transfer certificates pass.

The established numerical machinery is Craig–Bampton-style component mode synthesis and congruence projection. The proposed contribution is the governed admission architecture and its three-channel certificate, not the underlying component-mode mathematics.

## Model and benchmark

- Fine model: 193 nodes / 192 free mechanical DOFs.
- Interface partition: 12 fixed global ports.
- Component basis: exact static constraint modes plus nested fixed-interface vibration modes.
- Tests: 5 deterministic regimes and 24 fixed randomized holdouts.
- Drive bands: 3 kHz to approximately 330 kHz.
- Local first fixed-interface frequencies: approximately 0.59–0.86 MHz.

The drive band therefore remains below the first local fixed-interface resonances in these tests. Static constraint modes carry most of the response; a small number of local modes correct dynamic phase and transfer error.

## Evolution cycle A — transfer-only bridge

Acceptance condition:

\[
\epsilon_H=
\frac{\|H_r-H_f\|_2}{\|H_f\|_2}
\leq 3\times 10^{-3}.
\]

Results across 29 scenarios:

- 29/29 passed the transfer target.
- Median admitted modes: 6.
- Median globally coupled DOFs: 18 = 12 ports + 6 modes.
- Median global reduction: 90.625% relative to 192 free DOFs.
- Geometric-mean transfer mismatch: 0.2209%.
- Positive reduced stiffness and mass and nonnegative proportional damping in every case.

However, transfer norm alone allowed phase RMS errors approaching one degree in some cases. Averaging hid narrow resonant disagreement.

## Evolution cycle B — three-channel stability certificate

Define

\[
\Delta_{\mathrm{stability}}
=
\max\left(
\frac{\epsilon_H}{0.003},
\frac{\epsilon_\phi}{0.25^\circ},
\frac{\epsilon_\infty}{0.02}
\right),
\]

where:

- \(\epsilon_H\) is global complex transfer mismatch;
- \(\epsilon_\phi\) is response-weighted phase RMS;
- \(\epsilon_\infty\) is regularized worst-frequency mismatch.

The reduced model is certified only when

\[
\boxed{\Delta_{\mathrm{stability}}\leq1.}
\]

Results across the same 29-case structure:

- 29/29 passed all three gates.
- Median admitted local modes: 10.
- Median globally coupled DOFs: 22 = 12 ports + 10 modes.
- Median global reduction: 88.542%.
- Geometric-mean transfer mismatch: 0.09885%.
- Median phase RMS: 0.08165 degrees.
- Maximum phase RMS: 0.24391 degrees, below the 0.25-degree gate.
- Maximum worst-frequency mismatch: 1.9933%, below the 2% gate.
- At the same mode count, a global lowest-frequency modal allocation passed only 17/29 cases; transfer-sensitive triadic allocation passed 29/29 and had the lower normalized stability deviation in 18/29 cases.
- Every projected model retained positive mass and stiffness and nonnegative proportional damping.

## Browser proof

The integrated browser model uses a 97-node fine actuator divided into 12 components. It computes:

- exact local static constraint modes;
- local generalized fixed-interface modes;
- a shared 12-port transformation;
- projected mass, damping, stiffness, and force matrices;
- complex frequency responses;
- greedy coherence-driven modal admission.

Standalone JavaScript kernel results:

### Localized 3–60 kHz regime

- Ports: 12.
- Admitted modes: 0.
- Reduced global DOFs: 12 from 96.
- Transfer mismatch: 0.18682%.
- Phase RMS: 0.01777 degrees.
- Worst-frequency mismatch: 0.21287%.
- Certificate passed.

### Broadband 8–280 kHz regime

- Ports: 12.
- Admitted modes: 5.
- Reduced global DOFs: 17 from 96.
- Transfer mismatch: 0.25130%.
- Phase RMS: 0.24043 degrees.
- Worst-frequency mismatch: 0.76483%.
- Certificate passed.

## Emergent information

### 1. Stability lives in nestedness

Moving or deleting wave-carrying nodes produced non-monotone phase errors. Fixed ports and nested local modal spaces avoid changing the global interconnection geometry as the model evolves.

\[
V_0\subset V_1\subset V_2\subset\cdots
\]

Each admitted mode enlarges the local representational space without invalidating prior interface coordinates.

### 2. Residual is a proposal, not authority

A residual can rank missing local behavior, but admission is governed by its effect on observable transfer coherence. The authority sequence becomes

\[
\text{local defect}
\rightarrow
\text{candidate mode}
\rightarrow
\text{port transfer test}
\rightarrow
\text{admit or reject}.
\]

### 3. Average agreement is not stability

A small global norm can conceal a narrow phase defect near resonance. Stability requires simultaneous control of distributed error, phase, and extremal disagreement.

### 4. The compressed object is the exposed state

The physical model remains rich locally. What becomes sparse is the globally exchanged state:

\[
\boxed{
\text{local physical richness}
+
\text{invariant sparse ports}
+
\text{certified modal corrections}
}
\]

### 5. Proposed Triadic Port-Coherence Principle

A falsifiable architectural hypothesis is:

> A resonant system can be compressed stably when its physical manifold is preserved locally, its interface manifold remains invariant, and all enrichment is admitted through a multi-channel transfer-coherence certificate.

This is not yet a theorem or a universal scientific law. It is supported by the present 1D frequency-domain benchmarks and now requires a 2D axisymmetric PZT–bond–mirror test.

## Next decisive experiment

Construct a 2D axisymmetric coupled electromechanical model and compare:

1. full fine coupled FEM;
2. exact frequency-dependent port condensation;
3. fixed-port Craig–Bampton reduction;
4. triadic adaptive component-mode enrichment;
5. uniform global modal allocation at equal reduced dimension.

Measure cavity-length transfer, mirror surface phase, interface stress, electric charge response, reduced solve time, local factorization cost, and port certificate robustness under randomized bond thickness, stiffness, damping, and electrode geometry.
