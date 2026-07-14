# Evolution 6 — 4D Floquet–Optomechanical State Bridge

## Scope

This evolution adds a dimensionless reduced-order temporal layer to the existing 3D observer-conditioned laboratory. It is a validation and exploration model, not a calibrated prediction for a specific actuator.

The retained state is

\[
\mathbf y=(x,\dot x,\Re a,\Im a),
\]

where \(x\) is a normalized mechanical coordinate and \(a\) is the intracavity optical field.

The model combines:

\[
\ddot x+2\zeta\dot x+\left[1+h\cos(2\Omega t)\right]x+\beta x^3
=F\cos(\Omega t)+\mu |a|^2,
\]

\[
\dot a=\left[i(\Delta+g x)-\frac{\kappa}{2}\right]a+s.
\]

## Literature cross-check

The selected upgrades are grounded in established physics:

- Piezoceramics near resonance can exhibit Duffing-like jump phenomena, superharmonics, nonlinear damping, and quadratic/cubic constitutive terms.
- Floquet theory is the natural stability framework for periodically driven resonators near parametric instability.
- Cavity optomechanics includes dynamical backaction: the optical field has a finite relaxation time and can modify mechanical stiffness and damping.
- Periodically modulated cavity resonance can generate coherent Floquet sidebands, making a frequency/synthetic-dimension view physically relevant.
- Structure-preserving port-Hamiltonian reduction remains the appropriate future route for passive reduced models.

## Benchmark A — Floquet threshold

The damped Mathieu equation was integrated over one coefficient period and the monodromy spectral radius was evaluated.

For exact principal parametric resonance, the small-damping estimate is

\[
h_{\mathrm{crit}}\approx4\zeta.
\]

Across \(\zeta=0.005,0.01,0.02,0.04,0.07\), the maximum relative difference between the numerical threshold and \(4\zeta\) was:

\[
\boxed{0.1070\%}.
\]

The embedded 35 × 35 Floquet map covers frequency ratio 0.82–1.18 and modulation depth 0–0.34. Approximately 19.43% of that sampled map has positive Floquet growth.

## Benchmark B — nonlinear Duffing response

### Weakly nonlinear validation

First-harmonic balance was compared with direct RK4 integration in a single-valued weakly nonlinear regime.

Relative response error:

\[
\boxed{1.11\times10^{-6}}.
\]

### Strong nonlinear branch sweep

A phase-continuous upward/downward frequency sweep was compared with stable harmonic-balance branches.

- Upward sweep relative error: **0.2051%**
- Downward sweep relative error: **1.2395%**
- Harmonic-balance/direct hysteresis-area ratio: **0.99555**
- Maximum upward/downward branch separation: **0.87576** at frequency ratio approximately **1.1075**

This confirms that operating history becomes an independent state in the jump region.

## Benchmark C — cavity memory

A prescribed mechanical oscillation drove the complex cavity equation. Three observers were compared:

1. instantaneous adiabatic cavity map;
2. first derivative correction;
3. full dynamic optical quadrature state.

Under the declared gate

\[
\epsilon_{\mathrm{harmonic}}<5\%,\qquad
\epsilon_{\phi}<1^\circ,
\]

the derivative-corrected observer first passed at approximately

\[
\boxed{\kappa/\Omega=16}.
\]

Below this ratio, the optical amplitude and phase quadratures must remain dynamic states.

## Benchmark D — coupled 4D flow

The Duffing mechanical coordinate and optical quadratures were integrated together with radiation-pressure feedback.

The derivative-corrected reduced observer first met the coupled mechanical/optical gate at approximately

\[
\boxed{\kappa/\Omega=16}.
\]

At lower ratios it could reproduce mechanical motion more accurately than optical output, demonstrating that mechanical agreement does not certify the readout.

## Emerging temporal triad

The temporal certificate contains three independent channels:

\[
\boxed{
\Delta_T=
\max(\Delta_F,\Delta_B,\Delta_M)
}
\]

where:

- \(\Delta_F\): Floquet growth / parametric instability;
- \(\Delta_B\): nonlinear branch ambiguity and path dependence;
- \(\Delta_M\): optical observer-memory mismatch.

The emergent principle is:

> Time becomes a true model dimension when the observer retains information between mechanical cycles.

The fourth dimension is therefore not an animation coordinate. It is the retained optical state \((\Re a,\Im a)\), added when a memoryless observer fails its certificate.

## Laboratory changes

Evolution 6 adds:

- a 4D Floquet workbench opened with the **4D** header button or the `D` key;
- a precomputed Floquet stability tongue;
- direct and harmonic-balance Duffing branch plots;
- adiabatic versus derivative-corrected observer-memory curves;
- a live browser RK4 integrator for \((x,v,\Re a,\Im a)\);
- mechanical phase portrait, optical-quadrature orbit, and Poincaré points;
- stable, cusp-edge, and parametric-growth presets;
- a three-channel temporal certificate;
- 4D JSON export.

## Scientific boundary

The model is dimensionless and reduced-order. Its benchmark validates the numerical relationships implemented in the browser, not a specific PZT stack or optical cavity. A publishable next step requires parameter identification from a real actuator/cavity and a passive, causal reduced realization linked to the full 3D finite-element model.
