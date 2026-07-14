# Triadic Stability Bridge — 3D Evolution Findings

## Status

A full three-dimensional reduced-order experiment was constructed for a piezoelectric stack, compliant bond, and flexible mirror. The model uses trilinear eight-node hexahedral elasticity elements, three mechanical displacement coordinates per free node, consistent mass, Rayleigh damping, and prescribed inverse-piezoelectric eigenstrain.

The experiment validates a model-reduction architecture. It is not yet a mesh-converged physical prediction or an experimental piezoelectric device model.

## Governing full-order system

For harmonic angular frequency \(\omega\),

\[
\left(\mathbf K-\omega^2\mathbf M+i\omega\mathbf C\right)\mathbf u(\omega)=\mathbf f_{\mathrm{piezo}}.
\]

The representative PZT eigenstrain is

\[
\boldsymbol\varepsilon_p=
\begin{bmatrix}
 d_{31}\mathcal E_3 & d_{31}\mathcal E_3 & d_{33}\mathcal E_3 & 0 & 0 & 0
\end{bmatrix}^{\!T}.
\]

## The three levels

The 3D model established that one certificate must operate across three nested levels.

### Level 1 — Hidden volume physics

This level monitors internal strain-energy transfer and representative volumetric displacement sensors:

\[
\Delta_V=
\max\left(
\Delta_{\mathrm{sensor}},
\Delta_{\mathrm{energy}},
\Delta_{\mathrm{energy\ distribution}}
\right).
\]

Failure at this level authorizes additional fixed-interface volume modes.

### Level 2 — Compressed interface manifold

This level monitors transmission through component boundaries:

\[
\Delta_I=
\max\left(
\Delta_{\mathrm{port}},
\Delta_{\mathrm{port\ shape}},
\Delta_{\mathrm{bond\ energy}}
\right).
\]

Failure at this level authorizes additional response-informed port modes.

### Level 3 — Global mirror observables

This level monitors the actuator-to-mirror response:

\[
\Delta_G=
\max\left(
\Delta_{\mathrm{center}},
\Delta_{\mathrm{phase}},
\Delta_{\infty},
\Delta_{\mathrm{surface}},
\Delta_{\mathrm{shape}},
\Delta_{\mathrm{tilt\,x}},
\Delta_{\mathrm{tilt\,y}}
\right).
\]

Failure at this level promotes the worst hidden frequency into the sentinel set. When that frequency is already represented, it authorizes additional local modes.

The complete certificate is

\[
\boxed{
\Delta_{3D}=\max(\Delta_V,\Delta_I,\Delta_G)\leq1.
}
\]

## Two-cycle result

### Cycle 1 — shared enrichment authority

All missing information was forced through a largely shared mode budget.

- Cases: 10
- Passed: 3
- Median exposed states: 176
- Median global reduction: 85.07%
- Maximum certificate deviation: 4.6495

Cycle 1 often retained fewer states, but it did not know which representational layer owned the error.

### Cycle 2 — level-owned enrichment authority

Each failing level was given a distinct admissible action:

\[
\Delta_V>1\Rightarrow\text{admit volume mode},
\]

\[
\Delta_I>1\Rightarrow\text{expand port manifold},
\]

\[
\Delta_G>1\Rightarrow\text{promote resonance sentinel}.
\]

If a promoted sentinel was already present, global failure could request local enrichment instead.

Final results across five deterministic and five fixed asymmetric holdouts:

- Passed: **10/10**
- Full-order mechanical states per case: **1,179**
- Median exposed states: **249.5**
- Median port modes: **29**
- Median local volume modes: **225**
- Median global reduction: **78.84%**
- Maximum final deviation: **0.98531**
- Global level passed: **10/10**
- Interface level passed: **10/10**
- Volume level passed: **10/10**

The nominal broadband model used:

- 1,179 physical DOFs
- 621 original nodal interface coordinates
- 14 compressed port modes
- 253 local volume modes
- 267 globally exposed states
- 77.35% global reduction
- final certificate deviation 0.37018

## Emergent information

### 1. A triad is not three measurements; it is three authorities

The first cycle measured several quantities but routed all failures through essentially one enrichment mechanism. That was not stable in 3D.

The second cycle became stable when each level was allowed to change only the representation it governs.

> Measurement identifies disagreement. Authority determines which manifold may evolve.

### 2. Three-dimensionality creates split families

The 2D axisymmetric model suppresses torsional, tilted, and independently split \(x/y\) mirror families. The 3D model required materially more local modes because those families are real independent directions in the response manifold.

The reduction decreased from approximately 92% in the fine 2D experiment to approximately 79% median in the present 3D experiment. That is not a regression. It is the additional informational cost of genuine three-dimensional behavior.

### 3. Stable compression is hierarchical

The architecture is now:

\[
\boxed{
\text{rich volume physics}
\rightarrow
\text{compressed interface information}
\rightarrow
\text{certified global coherence}.
}
\]

This maps naturally—but interpretively, not as an established physical identity—to:

- Energy: hidden strain and material-energy distribution.
- Information: port and interface transmission.
- Coherence: globally observable mirror response.

### 4. The certificate must close from the inside outward

A global mirror match does not prove the bond or PZT energy field is correct. An accurate interface response does not prove the final mirror shape is correct. A locally accurate volume does not guarantee coherent global phase.

All three are required simultaneously.

## Physical-manifold convergence boundary

A separate comparison was run between the 1,179-DOF model and a 7,494-DOF model at eight frequencies.

Across the full 5–260 kHz sample:

- Center-transfer relative difference: 25.72%
- Phase difference: 34.89°
- Common-node mirror-surface difference: 78.78%
- Material-energy difference: 45.69%

At 5–9 kHz, center displacement differed by less than approximately 0.8% and the lowest-frequency surface field differed by approximately 1%. The higher-frequency radial, torsional, and thickness modes are not yet physically mesh converged.

Therefore:

\[
\boxed{
\text{ROM coherence is validated relative to the current full-order model,}
}
\]

but

\[
\boxed{
\text{broadband continuum accuracy is not yet certified.}
}
\]

This physical-convergence check is the foundation beneath the three-level reduction certificate; it is not replaced by it.

## Next evolution

The next model should use the 7,494-DOF physical manifold or a locally enriched equivalent, then introduce electric-potential degrees of freedom:

\[
\begin{bmatrix}
\mathbf K_{uu} & \mathbf K_{u\phi}\\
\mathbf K_{\phi u} & -\mathbf K_{\phi\phi}
\end{bmatrix}
\begin{bmatrix}
\mathbf u\\\boldsymbol\phi
\end{bmatrix}
-
\omega^2
\begin{bmatrix}
\mathbf M & 0\\0&0
\end{bmatrix}
\begin{bmatrix}
\mathbf u\\\boldsymbol\phi
\end{bmatrix}
=
\begin{bmatrix}
\mathbf f\\\mathbf q
\end{bmatrix}.
\]

The three-level controller should then govern:

1. electromechanical volume modes;
2. mechanical/electrical interface ports;
3. optical and mirror observables.

The falsifiable target is to preserve all three channels on a mesh-converged 3D manifold while exposing substantially fewer globally coupled electromechanical states.
