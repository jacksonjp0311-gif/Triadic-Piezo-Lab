# TRIAD // Classic Cinematic Scientific Evolution Laboratory

> **The vivid original 2035 resonant-mirror laboratory, restored as the canonical interface and packaged with the complete benchmark lineage.**

![Status](https://img.shields.io/badge/status-research%20preview-4dfcff)
![Version](https://img.shields.io/badge/version-8.0.0-ff38e4)
![Interface](https://img.shields.io/badge/interface-classic%20cinematic-ffd166)
![Distribution](https://img.shields.io/badge/distribution-one--click%20HTML-40ff9a)
![License](https://img.shields.io/badge/license-MIT-8f63ff)

**Creator and research direction:** James Paul Jackson  
**Release:** Evolution 8 — Classic Cinematic Evidence Package  
**Date:** July 2026

---

## Open the laboratory

### Fastest path

Double-click:

```text
OPEN_LAB.html
```

The laboratory uses pinned Three.js modules from jsDelivr. The initial load therefore requires an internet connection.

### Reliable one-click launch

**Windows:** double-click `Launch-Lab.cmd`  
**PowerShell:** run `Launch-Lab.ps1`  
**macOS / Linux:** run `./launch-lab.sh`

The launchers start a small local web server when Python is available, then open the laboratory automatically.

---

## What changed in v8

This release restores the original visual system as the canonical design:

- vivid cyan, magenta, violet, teal and gold laboratory lighting;
- bright glassmorphic side panels rather than the later near-black workbench;
- the animated PZT stack, electrodes, optical cavity, mirror and vacuum field;
- dense live metrology, plots and scientific controls;
- the integrated 2D axisymmetric bridge;
- the original cinematic boot sequence and emergent-analysis panel.

It adds two controls without changing that design language:

### Independent display matrix

Press the **sun button** in the header or press `B`.

You can independently adjust:

- **UI brightness** — dims panels, text and HUD without hiding the specimen;
- **UI color intensity** — lowers neon saturation;
- **specimen exposure** — controls Three.js tone-mapping exposure;
- **holographic bloom** — controls glow without changing physical data.

`Shift+B` restores the original classic profile.

### In-app evidence matrix

Press the **Σ button** or `E` to open a concise benchmark lineage containing accepted results, bounded claims, rejected hypotheses and scientific limitations.

---

## Laboratory capabilities

### Interactive resonant mirror chamber

- multilayer PZT actuator and electrodes;
- flexible circular mirror and cavity field;
- real-time damped modal response;
- circular Bessel-like mirror modes;
- frequency response and oscilloscope plots;
- qualitative dynamical-Casimir-inspired particle visualization;
- component inspection and floating holographic annotations.

### Adaptive and reduced-order experiments

- residual-governed axial FEM;
- two-cycle adaptive refinement and equal-DOF safety gate;
- exact dynamic condensation;
- fixed-port component-mode stability bridge;
- 2D axisymmetric PZT–bond–mirror continuum;
- transfer, geometry and energy certification.

### Packaged evidence lineage

The archive also carries the later benchmark results that extended the architecture into:

- three-level 3D volume/interface/global control;
- observer-conditioned electrical and optical witnesses;
- nonlinear Duffing branches;
- Floquet parametric stability;
- optical memory as a fourth state dimension.

Those later datasets are packaged as evidence without replacing the restored classic interface.

---

## Benchmark highlights

| Evolution | Reproduced or archived result | Boundary discovered |
|---|---:|---|
| Adaptive FEM | Localized target reached with 20 elements where the tested uniform sequence required 48 | Quasi-static scalar pulses add no new spatial ranking; broadband waves may require global resolution |
| Triadic multigate | 29/29 passed; median 88.54% global reduction | Transfer average alone can hide phase and extremal error |
| Axisymmetric 2D | 15/15 passed; median 134 exposed states from 1,756 | Reduced-model coherence is distinct from physical mesh convergence |
| Full 3D | 10/10 passed; median 249.5 exposed states from 1,179 | Torsion and split mode families require level-owned enrichment |
| Observer manifold | Base ROM passed 5.09%; witness closure passed 100% of 864 tested points | Validity belongs to the model–observer–operating-envelope pair |
| 4D temporal | Floquet threshold max error 0.107%; hysteresis-area ratio 0.99555 | Time becomes state when growth, branch history or optical memory persists |

Negative controls and failed branches are preserved under `benchmarks/results/negative_controls/`.

---

## Reproduce and verify

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

Verify archived benchmark JSON and package structure:

```bash
python benchmarks/run_benchmarks.py --quick
python tests/verify_package.py
```

Rerun the portable 4D Floquet/Duffing/optical-memory benchmark:

```bash
python benchmarks/run_benchmarks.py --reproduce-4d
```

The largest 2D, 3D and observer suites are shipped with their original source and recorded outputs. They are computational research runs rather than mandatory startup tasks.

---

## Controls

| Key | Action |
|---|---|
| `B` | Open/close display matrix |
| `Shift+B` | Restore original cinematic brightness |
| `E` | Open/close benchmark evidence |
| `2` | Open/close 2D axisymmetric laboratory |
| `Space` | Start/stop oscillation |
| `R` | Reset laboratory state |
| `P` | Random preset |
| `F` | Fullscreen |
| `Esc` | Close active overlays |

---

## Package layout

```text
Triadic-Piezo-Scientific-Evolution-Lab-v8.0.0-Classic-Cinematic/
├── OPEN_LAB.html
├── Launch-Lab.cmd
├── Launch-Lab.ps1
├── launch-lab.sh
├── README.md
├── QUICKSTART.txt
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
├── manifest.json
├── SHA256SUMS.txt
├── requirements.txt
├── benchmarks/
│   ├── run_benchmarks.py
│   ├── scripts/
│   ├── reference_sources/
│   └── results/
│       ├── benchmark_summary.json
│       ├── adaptive_fem/
│       ├── stability_bridge/
│       ├── axisymmetric_2d/
│       ├── triadic_3d/
│       ├── observer_manifold/
│       ├── evolution6_4d/
│       └── negative_controls/
├── docs/
├── research/
├── source/
│   └── ORIGINAL_CLASSIC_STYLE_REFERENCE.html
└── tests/
```

---

## Design invariants

1. **One-click remains sacred.** The distributable application is one HTML file.
2. **The classic visual language is canonical.** New features must fit its cinematic glass-and-neon system.
3. **UI and specimen luminance remain independent.** Accessibility does not require flattening the physical scene.
4. **Compression targets exposed coupling—not the manifold carrying coherence.**
5. **Every level owns its correction.** Residuals propose; certificates authorize.
6. **Observers participate in validity.**
7. **Failures remain evidence.**
8. **No equation enters only for appearance.**

---

## Scientific scope

TRIAD is a numerical research simulator and scientific visualization package. It is not experimental calibration of a manufactured device, proof of a universal physical law or a production engineering qualification tool. High-frequency 2D and 3D continuum results retain documented mesh-convergence limits. The dynamical-Casimir particle layer is qualitative.

Read `docs/SCIENTIFIC_BOUNDARIES.md` before publishing quantitative claims.

---

## Citation

```text
Jackson, James Paul. TRIAD Classic Cinematic Scientific Evolution
Laboratory, v8.0.0 (2026). Interactive piezoelectric, reduced-order,
observer-conditioned and temporal-stability research package.
```

---

## License

MIT License. See `LICENSE`.
