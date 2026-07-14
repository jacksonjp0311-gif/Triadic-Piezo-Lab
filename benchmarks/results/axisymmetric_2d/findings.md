# Axisymmetric Triadic Stability Bridge — 2D Evolution Findings

## Scope

This study extends the prior 1D piezoelectric stability bridge into a two-dimensional axisymmetric meridional continuum. The model contains a multilayer PZT stack, compliant bond, and overhanging flexible mirror. It uses four-node axisymmetric quadrilateral elements with radial and axial displacement, hoop strain, consistent mass, Rayleigh damping, and prescribed piezoelectric eigenstrain.

This is a computational reduced-order-model study. It is not experimental validation, a theorem, or proof of a new physical law.

## Full-order system

For harmonic frequency \(\omega\),

\[
\mathbf Z(\omega)\mathbf u(\omega)=\mathbf f_p,
\qquad
\mathbf Z(\omega)=\mathbf K-\omega^2\mathbf M+i\omega\mathbf C.
\]

The PZT region receives the imposed strain

\[
\boldsymbol\varepsilon_p=
\begin{bmatrix}
 d_{31}E_3 & d_{33}E_3 & d_{31}E_3 & 0
\end{bmatrix}^{\!T}.
\]

The fine benchmark mesh contains:

- 832 active Q4 elements;
- 1,756 free mechanical degrees of freedom;
- 280 nodal interface coordinates before port compression;
- 14 material/geometric subcomponents.

## Evolution 1 — nodal-port component synthesis

The first 2D implementation retained every shared interface coordinate and added local fixed-interface modes. On the 462-state development mesh:

- center-transfer-only cycle: 5/15 scenarios passed;
- transfer–geometry–energy cycle: 7/15 scenarios passed;
- median retained local modes: 42;
- median global reduction: 60.61%.

The interface itself had become the dimensional bottleneck.

## Evolution 2 — response-informed interface manifold

The nodal interface coordinates were replaced by a low-dimensional port basis derived from complex interface response snapshots. Local component modes remained nested and fixed-interface.

On the 462-state development mesh, full-band fitting passed all 15 scenarios with:

- median port modes: 8;
- median local modes: 45;
- median globally coupled states: 53;
- median global reduction: 88.53%.

However, an alternating-frequency holdout passed only 1/15 scenarios. Snapshot energy could fit sampled frequencies while missing resonances between them.

## Evolution 3 — adversarial resonance closure

Define the reduced space

\[
V_{S,q}=\operatorname{span}\!\left(\mathbf T_p\mathbf P_S,\,\boldsymbol\Phi_q\right),
\]

where:

- \(S\) is the sentinel-frequency set;
- \(\mathbf P_S\) is the response-informed port basis;
- \(\boldsymbol\Phi_q\) contains nested local fixed-interface modes.

The 2D certificate has three outer channels:

\[
\Delta_{2D}=\max\left(\Delta_T,\Delta_G,\Delta_E\right),
\]

where:

- \(\Delta_T\): complex tip transfer, weighted phase, and worst-frequency deviation;
- \(\Delta_G\): mirror-surface field and normalized radial shape disagreement;
- \(\Delta_E\): total strain-energy and material-energy-distribution disagreement.

A model is admitted only when

\[
\Delta_{2D}\le 1.
\]

If the certificate fails, the worst violating unsampled frequency is added:

\[
S_{k+1}=S_k\cup
\left\{
\arg\max_{\omega\in\Omega}\Delta(\omega;V_{S_k,q_k})
\right\}.
\]

The port basis and local-mode allocation are then rebuilt. This repeats until no frequency in the declared band violates the triad.

## Fine-manifold benchmark

Fifteen scenarios were tested: five deterministic regimes and ten fixed randomized geometries/material systems.

All 15 passed on the 1,756-state physical model.

| Metric | Result |
|---|---:|
| Passed scenarios | 15/15 |
| Median resonance sentinels | 8 |
| Maximum sentinels | 9 |
| Median port modes | 16 |
| Median local modes | 120 |
| Median globally coupled states | 134 |
| Minimum globally coupled states | 46 |
| Maximum globally coupled states | 276 |
| Median global reduction | 92.369% |
| Minimum global reduction | 84.282% |
| Maximum stability deviation | 0.98454 |
| Geometric-mean complex transfer error | 0.28963% |
| Geometric-mean mirror-surface error | 0.21341% |
| Geometric-mean material-energy error | 0.48123% |
| Geometric-mean phase RMS | 0.07687° |

The nominal broadband model retained:

- 1,756 physical states;
- 280 nodal interface coordinates before compression;
- 14 response-informed port modes;
- 105 local fixed-interface modes;
- 119 globally coupled states;
- 93.223% global reduction;
- stability deviation 0.82146.

## Physical mesh convergence

Reduction coherence and finite-element convergence are separate requirements.

Level 1 (462 DOFs) versus Level 2 (1,756 DOFs) showed large broadband discretization changes and was rejected as the final physical baseline.

Level 2 versus Level 3 (3,882 DOFs) produced:

| Regime | Tip transfer | Phase RMS | Surface field | Material energy |
|---|---:|---:|---:|---:|
| Nominal local | 0.979% | 0.0605° | 7.374% | 2.489% |
| Nominal broadband | 1.901% | 0.7243° | 14.501% | 4.792% |

Therefore the 1,756-state model is substantially improved but the broadband surface field is not yet fully mesh converged. The present results validate model reduction relative to that physical manifold; they do not establish continuum-exact mirror deformation.

## Emergent information

### 1. Interface invariance is a subspace property

In 2D, retaining every nodal interface coordinate is unnecessary. The invariant object is the transmissible interface response space, not a particular nodal coordinate list.

### 2. Stable compression requires two nested manifolds

\[
\boxed{
\text{response-informed interface manifold}
\;\oplus\;
\text{nested local physical manifold}
}
\]

The first communicates between components. The second preserves hidden material and wave physics.

### 3. Stability is closure, not fit

A basis that matches sampled frequencies can still miss an unsampled resonance. Band stability is achieved only after the sentinel set is closed against its own worst counterexample.

### 4. Physical refinement need not cause proportional global growth

Refining the physical mesh from 462 to 1,756 free DOFs increased the median certified global state count from 96 in the coarse adversarial study to 134 in the fine study. The hidden local manifold grew substantially, but globally exposed complexity remained small relative to the physical model.

## Proposed principle

> A component-based resonant continuum is stably compressible over a declared frequency band when its interface response subspace is closed against adversarial resonance, its local fixed-interface spaces are nested, and transfer, surface geometry, and energy localization pass simultaneously.

This is a falsifiable computational principle supported by the present benchmark family. The next required step is a mesh-converged axisymmetric model followed by a fully coupled electromechanical formulation with electric-potential degrees of freedom and experimental material parameters.
