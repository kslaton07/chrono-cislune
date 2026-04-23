[EDITED] Chrono::FSI-SPH Modified Cam Clay (MCC) Parameter Guide {#manual_fsi_mcc_rheology}
====================================================================================

\tableofcontents

Scope
-----

This page is a practical guide for configuring the CRM `MCC` rheology in Chrono::FSI-SPH.

It focuses on:

- physical intuition for each MCC parameter
- exact API calls used in Chrono
- code-level behavior of the return-mapping implementation
- two calibration workflows:
  - classical geomechanics test-based identification
  - empirical Chrono calibration when only interaction tests are available

Primary references:

- Itasca Cam-Clay documentation (parameter interpretation and classical calibration):  
  https://docs.itascacg.com/itasca900/common/models/camclay/doc/modelcamclay.html
- Chrono MCC implementation files:  
  `src/chrono_fsi/sph/physics/SphFluidDynamics.cu`,  
  `src/chrono_fsi/sph/physics/SphForceWCSPH.cu`,  
  `src/chrono_fsi/sph/physics/SphDataManager.cu`

Motivation Behind MCC - THE MEMORY OF TERRAIN (PRESSURE)
--------------------------------------------------------
In general, terrain has "memory", past interactions alter its composition. Most importantly, the application of pressure can permanently alter the volume of the terrain (i.e. how much it has been compacted).

Below is a "Compression Plane" or "v-ln(p)" plot. It is used to track the behavior described above. It depicts normal consolidation (i.e. soil that hasn't been compressed previously).

<img src="http://www.projectchrono.org/assets/manual/Chrono_FSI-MCC_NCL.png" width="800">

FIGURE TITLE: Normal consolidation and swelling behavior in `v-ln(p)` space.

As you increase pressure (x-axis), volume decreases (y-axis), shown by the line labeled NCL. If at any point you remove the pressure, the terrain partially restores its volume. This is shown by the "swelling lines."

In actuality, there is a swelling line for every point on the NCL — all parallel, sharing an identical slope of _kappa_. However, only one matters at any given time: the one at _p_c_, the maximum pressure the terrain has ever experienced. _p_c_ determines which swelling line the terrain currently follows. Keep it in mind, as it reappears in the next section!


The slopes of the NCL and swelling line (_lambda_ and _kappa_, respectively), determine the overall deformation behavior of the terrain.
If the terrain were purely elastic, _kappa_ = _lambda_, meaning the volume would fully restore upon removing pressure. However, in the figure above, since _kappa_ < _lambda_, the terrain undergoes plastic (permanent) deformation. Overall, the larger the gap between _lambda_ and _kappa_, the greater the irreversible compaction.

Motivation Behind MCC - THE MEMORY OF TERRAIN (PRESSURE AND SHEAR) 
-----------------------------------------------------------------
In the previous section, we only considered normal pressure (like squeezing the terrain) and how that impacts its composition. In reality, there are two types of stress that can cause terrain to permanently deform. Here we introduce shear stress, the second phenomenon that can alter soil composition.

For example, imagine a stick being pushed into sand. The sand directly below the stick compresses due to normal pressure, however the sand along the sides of the stick is forced to "slide past it." This sliding force is shear stress.

Below is a "p-q Plot" or "Stress Space Diagram", which captures how both mean pressure _p_ (squeezing terrain equally from all sides) and shear stress _q_ (sliding or distorting the soil) affect terrain composition.

<img src="http://www.projectchrono.org/assets/manual/Chrono_FSI-MCC_yield.png" width="800">

FIGURE TITLE: MCC elliptical yield surface in `p-q` space with CSL slope `M` and cap size `p_c`.

The yield curve characterizes the type of deformation the terrain experiences based on its current _p_ and _q_ combination. Note that this curve is an ellipse, which is a modeling choice of MCC to simply capture how combined pressure and shear drive terrain deformation. If the current _p_-_q_ combination is inside the ellipse, the terrain deforms elastically. If it is outside, the terrain yields plastically. The ellipse ends at _p_c_ on the right — the same "history pressure" from before.

The CSL is where the soil fails due to shear stress alone. _M_ is the slope of this line, defined such that the CSL always passes through the top of the ellipse
(together, _M_ and _p_c_ determine the shape of the ellipse). This intersection between the yeild curve and the CSL captures the stress state where the combination of pressure and shear is just enough to cause plastic deformation.

Connecting back to the previous section: when a _p_-_q_ combination causes the soil to yield, _p_c_ increases and the ellipse grows wider to reflect the new stress history. The shape of the ellipse also varies between terrain types, as _M_ differs from soil to soil.

Chrono Implementation: All MCC Terrain Parameters and Intuition
-----------------------------------------------------------------
Below is a summary of all MCC simulation parameters and how they influence terrain behavior.

- **lambda**: Slope of the NCL. Controls how much the terrain permanently compacts under virgin compression. Larger lambda = softer, more compressible terrain.

- **kappa**: Slope of the swelling lines. Controls how much the terrain elastically recovers after pressure is removed. Larger kappa = more elastic springback. Must be smaller than lambda.

- **M**: Slope of the CSL in p-q space. Controls shear strength — how much sliding/distorting force the terrain can resist before failing. Larger M = stronger resistance to shear.

- **p_c** (per particle): The preconsolidation pressure — the maximum pressure a particle has ever experienced. Sets the size of the yield ellipse. Larger p_c = larger elastic zone, stiffer early response.

NOTE: This may seem inconsistent with the previous section, where p_c described the entire terrain. In reality, the NCL is a material-wide property, but each particle tracks its own p_c based on its individual stress history. Think of a sandbox — the whole sandbox shares the same NCL, but only the portion you pressed has an elevated p_c.

- **v_lambda**: The specific volume of the terrain at a reference pressure of 1000 Pa on the NCL. Essentially sets the vertical position of the NCL on the compression plane. Determined from lab compression data.

- **Young's Modulus**: How stiff the material is under stretching/compression. Under elastic deformation in MCC, this helps set the rate of the elastic response. 

- **Poisson's Ratio**: When you compress something in one direction, it tends to expand sideways, Poisson's ratio describes how much it expands/buldges. It works together with young's modulus to capture elastic behavior (something buldging out but restoring its intial composition). 

- **density**: Bulk density of the terrain. Affects inertia (heavier terrain accelerates slower) and helps determine the initial pressure at each particle (i.e influencing the initial _p_c_ and setting a particles starting position on the NCL). 

Chrono Implementation: General Workflow for Setting Up MCC Terrain
------------------------------------------------------------------

1. Set `rheology_model = RheologyCRM::MCC`.
2. Set all MCC constitutive parameters explicitly: `mcc_M`, `mcc_kappa`, `mcc_lambda`, `mcc_v_lambda`.
3. Initialize per-particle preconsolidation pressure `p_c0` (via `pc` in `AddSPHParticle`, or callback-based scaling).
4. Reuse SPH solver settings from `manual_fsi_sph_parameter_selection` (especially artificial viscosity, shifting, free-surface handling, and time-step strategy).

Chrono Implementation: Specfic API Map for Setting Up MCC Terrain Parameters
----------------------------------------------------------------------------

Constitutive parameters are set through `ElasticMaterialProperties`:

~~~{.cpp}
ChFsiFluidSystemSPH::ElasticMaterialProperties mat_props;
mat_props.rheology_model = RheologyCRM::MCC;

mat_props.density = 1700.0;
mat_props.Young_modulus = 1.0e6;
mat_props.Poisson_ratio = 0.3;

mat_props.mcc_M = 1.1;        // CSL slope
mat_props.mcc_kappa = 0.02;   // swelling/recompression slope
mat_props.mcc_lambda = 0.10;  // NCL slope
mat_props.mcc_v_lambda = 2.0; // specific volume at p1 = 1000 Pa

sysSPH.SetElasticSPH(mat_props);
~~~

Important: `p_c0` is not in `ElasticMaterialProperties`. It is per particle and passed through `pc` at particle creation:

~~~{.cpp}
// Explicit particle creation path
double p0 = rho0 * g * depth;
double pc0 = pre_pressure_scale * p0;
sysSPH.AddSPHParticle(pos, rho0, p0, mu0, v0, tau_diag, tau_offdiag, pc0);
~~~

If you use `ChFsiProblemSPH` with a properties callback:

~~~{.cpp}
class Props : public ChFsiProblemSPH::ParticlePropertiesCallback {
  public:
    void set(const ChFsiFluidSystemSPH& sysSPH, const ChVector3d& pos) override {
        p0 = sysSPH.GetDensity() * std::abs(sysSPH.GetGravitationalAcceleration().z()) * (z0 - pos.z());
        rho0 = sysSPH.GetDensity();
        mu0 = sysSPH.GetViscosity();
        v0 = VNULL;
        pre_pressure_scale0 = 2.0;  // ChFsiProblemSPH sets pc = p0 * pre_pressure_scale0
    }
    double z0;
};
~~~

Code-level note:

- In `ChFsiProblemSPH.cpp`, `pc` is built as `p0 * pre_pressure_scale0`.
- In direct `AddSPHParticle`, you provide `pc` explicitly.

Chrono Implementation: MCC Math-to-Code Map
-------------------------------------------

Chrono's MCC branch is implemented in `TauEulerStep(...)` (`SphFluidDynamics.cu`) and uses stress-rate terms from
`CrmRHS(...)` (`SphForceWCSPH.cu`).

| Concept | Implementation detail |
|---|---|
| Trial stress | `sig_N_tr = sig_n + dt * deriv_sig` in `TauEulerStep(...)` |
| Invariants | `p = -(1/3)tr(sigma)`, `q = sqrt(3 J2)` |
| Yield function | `f = q^2 + M^2 p (p - p_c)` |
| Elastic moduli | `K = v p / kappa`, `G = 3K(1-2nu)/(2(1+nu))` (then clamped) |
| Plastic correction | quadratic solve for `delta_lambda`, return to yield surface |
| Hardening | `p_c <- p_c * (1 + delta_eps_vp * v / (lambda - kappa))`, with floor |
| Tension handling | if trial/corrected `p < 0`, stress is zeroed |
| Free-surface handling | stress and pressure zeroed for free-surface flagged particles |

Current Chrono safeguards (important for interpretation):

- `K` and `G` are clamped to `[0.1, 1.0]` of reference values derived from `E` and `nu`.
- `p_c` has a hard floor (`>= 100 Pa` in current code).
- Free-surface particles skip hardening updates and are stress-zeroed.
- EOS is not used in CRM; pressure is obtained from stress (`p = -(1/3)tr(sigma)`).

Chrono Implementation: Core MCC Parameters Variable Names and Intuition (Again) 
-------------------------------------------------------------------------------

Set these through `ElasticMaterialProperties`, plus per-particle `pc`.

| Parameter | Physical meaning | Code default | Practical recommendation | API |
|---|---|---|---|---|
| `mcc_M` | Critical state line slope (`q = M p`), controls shear strength at critical state. | `0` | Set from triaxial data or friction angle mapping. Typical terramechanics starts are about `0.8-1.4`. | `mat_props.mcc_M` -> `SetElasticSPH` |
| `mcc_kappa` | Swelling/recompression slope in `v-ln(p)`; controls elastic volumetric stiffness (`K ~ v p / kappa`). | `0` | Must be positive and usually smaller than `lambda`; increasing `kappa` softens elastic reload response. | `mat_props.mcc_kappa` -> `SetElasticSPH` |
| `mcc_lambda` | Normal consolidation slope in `v-ln(p)`; controls virgin compressibility and hardening rate. | `0` | Must exceed `kappa`; larger `lambda` gives more compressible virgin response. | `mat_props.mcc_lambda` -> `SetElasticSPH` |
| `mcc_v_lambda` | Specific volume intercept at reference pressure `p1 = 1000 Pa`. | `2.0` | Set from isotropic compression line intercept in `v-ln(p)` space. | `mat_props.mcc_v_lambda` -> `SetElasticSPH` |
| `pc` (per particle) | Initial preconsolidation pressure `p_c0`; sets initial yield cap size / OCR effect. | `1e3` in `AddSPHParticle` signature | Set from OCR/test data or use depth-dependent `p0 * pre_pressure_scale`. This is highly influential in penetration/sinkage response. | `AddSPHParticle(..., pc)` or callback `pre_pressure_scale0` |
| `density` | Bulk density scale for inertia and pressure/stress coupling. | `1000` | Use measured bulk density for intended initial state. | `mat_props.density` -> `SetElasticSPH` |
| `Young_modulus`, `Poisson_ratio` | Reference elastic moduli used to compute `K_bulk` and `G_shear` and clamping bounds. | `1e6`, `0.3` | Keep physically plausible; very large `E` can force smaller stable steps in explicit integration. | `mat_props.Young_modulus`, `mat_props.Poisson_ratio` -> `SetElasticSPH` |

Critical consistency checks:

- `mcc_kappa > 0`
- `mcc_lambda > mcc_kappa` (hardening denominator uses `lambda-kappa`)
- `mcc_M > 0`
- `pc > 0` and initial `p0 > 0` where NCL-based initialization is used

Naming note:

- In current MCC equations in Chrono, `kappa` is used as the swelling/recompression slope and `lambda` as the NCL slope.
  This follows standard MCC notation used in the equations and implementation.

Traditional Parameter Selection (Physical Tests)
------------------------------------------------

This is the standard CSSM/MCC workflow used in geomechanics (see Itasca reference above).

1. Isotropic compression + unload/reload tests:
   - Fit `v-ln(p)` lines.
   - `lambda` = slope of normal consolidation line (NCL).
   - `kappa` = slope of swelling/recompression lines.
   - `v_lambda` = intercept at reference pressure `p1` (Chrono uses `p1 = 1000 Pa` in initialization code).
2. Triaxial tests to determine critical state slope:
   - `M = q/p` at critical state.
   - If only friction angle `phi` is available (triaxial compression convention), a common mapping is:
     `M = 6 sin(phi) / (3 - sin(phi))`.
   - Equivalent Itasca form:
     `N_phi = (1 + sin(phi)) / (1 - sin(phi))`, `M = 3 (N_phi - 1) / (N_phi + 2)`.
3. Initial cap size from stress history:
   - Set `p_c0` from isotropic yield stress or OCR:
     `p_c0 = OCR * p0`.

Useful conversion if lab data is in void ratio:

- `v = 1 + e`
- if using compression/recompression indices from `e-log10(p)` data, typical conversions are:
  `lambda = Cc / 2.303`, `kappa = Cr / 2.303`

Empirical Chrono Calibration Workflow (When Lab Data Is Missing)
-----------------------------------------------------------------

1. Set `M` first from a target friction angle (or a prior material estimate).
2. Set `kappa`, `lambda`, and `p_c0` to match penetration/sinkage compliance and hardening trends.
3. Use complementary tests (for example Cone Penetrometer Test + Normal Bevameter Plate Sinkage Test) to reduce parameter ambiguity.
4. Keep SPH numerics fixed while tuning constitutive parameters (see `manual_fsi_sph_parameter_selection`).

Interpretation while tuning:

- Increase `lambda`: softer virgin compression, larger settlement.
- Increase `kappa`: softer elastic unload/reload response.
- Decrease `(lambda-kappa)`: faster `p_c` growth under plastic compaction (faster hardening).
- Increase `p_c0`: larger initial elastic domain (stiffer early response).
- Increase `M`: higher shear resistance for given mean pressure.

Current Defaults and Demo Ranges
--------------------------------

Current constructor defaults for MCC-specific fields are:

- `mcc_M = 0`
- `mcc_kappa = 0`
- `mcc_lambda = 0`
- `mcc_v_lambda = 2.0`

These defaults are placeholders and are not physically usable for MCC without explicit user input.

Values seen in Chrono demos span a broad range depending on scenario scale:

- `mcc_kappa`: about `0.01` up to `0.2`
- `mcc_lambda`: about `0.02` up to `1.0`
- `M`: often computed from a chosen friction coefficient/angle
- `pre_pressure_scale` for `p_c0`: from about `2` (large-scale terrain) to much higher in penetration-focused setups

Hard-coded Specific-Volume Initializer (Current Code)
-----------------------------------------------------

Chrono currently initializes specific volume (`Sv`) in `FsiDataManager::AddSphParticle` using:

- `Sv = v_lambda - lambda * ln(pc / p1) + kappa * ln(pc / confining_stress)`
- with hard-coded `p1 = 1000 Pa`
- and `confining_stress = pres` (the particle's initial pressure)

Location:

- `src/chrono_fsi/sph/physics/SphDataManager.cu`

Implications:

- This initializer is always executed at particle insertion.
- It assumes positive `pc` and positive `pres` (to avoid invalid logarithms).
- If you initialize particles with near-zero pressure or leave `pc` at nonphysical values, MCC state initialization can be inconsistent.

Practical guidance:

1. In MCC runs, initialize particles with physically meaningful `p0` and `pc0`.
2. Keep `pc0` consistent with your intended OCR/state profile (depth-dependent if needed).
3. Set `mcc_v_lambda` explicitly when you have `v-ln(p)` calibration data.

Common Failure Modes and First Fixes
------------------------------------

- Very soft early response:
  increase `p_c0` and/or reduce `lambda`.
- Overly stiff response with little plasticity:
  reduce `p_c0`, increase `lambda`, or reduce `M`.
- Unstable volumetric behavior:
  verify `lambda > kappa > 0`, and reduce step size (or use variable time step).
- Free-surface particles carrying unphysical stress:
  check `free_surface_threshold` and keep CRM free-surface handling enabled.
