# aaPhysicalSky review

Read from source. Only the struct layout in section 3 was compiled. Nothing was rendered or compared with reference values.

Line numbers in `physicalSkyClass.h` refer to the working tree with the uncommitted `addSun` change (`Vector3`, `#include "vector.h"`). HEAD is one line lower from line 30 on.

## 1. Units

### Correct
- Sun table in `spectralData.h` is W·cm⁻²·sr⁻¹·µm⁻¹ radiance. The chain ×1e4, ÷1e3, ×Ω (`physicalSunClass.h:146-148`) gives W/m²/nm irradiance.
- Sun ×683.002 gives lux. `sunDiscCol` = lux/Ω is about 1.6e9 cd/m².
- Perez zenith Y: kcd/m² ×1000 is correct.
- Perez and zenith coefficient tables match the paper.
- β_m is near 1.2e-5 m⁻¹ at 550 nm. β_p is plausible: the `pow(0.01, V-3)` factor gives about 8e-5 m⁻¹ at T = 2. Without it, β_p is 100× too big. The AP maths is in SI metres.
- Ground colour ×`candelaConversion` ÷`multiplier` cancels correctly.

### Wrong or mixed
1. `ap_CalculateS0` (`physicalSkyClass.h:688-690`):
   - `delTheta` is 4.5° and `delPhi` is 9°. The comments say 9° and 18°.
   - 20 × 4.5° covers 90° of zenith angle. This is correct for the sky dome.
   - 20 × 9° covers only 180° of azimuth. Half the sky is missing.
   - `delThetaPhi` is in degrees², about 3283× too big.
2. Sun and sky inscatter weights are not photometric-consistent. `getSkySpectralRadiance` scales the spectrum so that ∫S·ȳ is cd/m². So the sun term (W/m²/nm) needs ×683.002, but it gets ×100. The sky term needs ×1 with radians and the full 2π, but it gets ×1e-4. Effective factors: sky 3283 × 1e-4 × ½ ≈ 0.16, sun 100/683 ≈ 0.15. The ratio is nearly right. The total is about 6× too weak.
3. `toXYZ` (`spectrum.h:373-379`) uses the 10° CMF (columns 4-6). Preetham x,y and 683 lm/W are 2° observer. Use columns 1-3.
4. No chromatic adaptation. `XYZtoACEScg` uses the D60-referenced AP1 matrix. The sky is D65. Its comment still says "sRGB primaries and D65". Match the host renderer's XYZ→working-space convention, not this shader's. `clampRGB` hard-zeroes negatives instead of gamut-mapping.
5. `h0`, `s = 32000`, `Rl` and `max_visibility_distance` are treated as metres. No conversion from scene units.
6. `h0 = -camera_height` makes air denser with altitude, since e^(−αh0) = e^(+α·height). Harmless at 1 m, wrong at 1 km.
7. `ap_InitS0` builds the inscatter table with h0 = +1 (`ap_ExtinctionFactor(1.f, …)`). `shader_evaluate` evaluates it with h0 = −camera_height. The two do not agree.
8. Below-horizon inscatter integrates underground to s = 32 km with signed cosθ. Extinction uses `fabs`. e^(α₁·32001) is about 4e11 for Mie. `Helper1`'s `t` probably underflows to 0, not NaN. Check ground rays by running.
9. `addSun` uses `half_dot = dot*0.5` in `rescale`. Limb darkening and blur are about constant across the disc (blur ≈ 0.785, limb ≈ 0.95). The disc is about 25% dimmer than the light's irradiance, with no limb darkening. `sun_disk_scale` is read and stored but never used.

## 2. Model bugs
- `perez()` line 400: the denominator uses `cos2gamma`. Preetham needs cos²θs. Correct only at the zenith.
- CIE overcast (`:429`): Y = Yz(1 + cosθ/3). The CIE law is Yz(1 + 2cosθ)/3. `Y_occ` is computed and not used.
- Rayleigh phase (`:707-709`, `:727-729`): `psi` is already a dot product (cosψ). The code then takes `cos(psi)`, so the phase is 1 + 0.9324·cos²(cosψ). This inverts the shape: the maximum is at 90° instead of the minimum. `ap_GetNeta(psi)` is correct because it does `acos` itself.
- Sky radiance is attenuated twice. Preetham L is already the radiance seen at the ground.
  - `ap_CalculateS0` multiplies L by `ap_ExtinctionFactor(1, 32000)` before it integrates.
  - `shader_evaluate` applies Fex(32 km) to sky misses and adds Lin.
- `arnoldSkyShader.cpp:258`: `_xyY.y = overcastWhitePoint_D65_x` sets y to x.
- `addGround:486`: `AiV3Cross(v_cross, v_camDir)` reads an uninitialised `v_cross` (`AtVector() = default` in Arnold 7.4) and assigns the result to `v_h_blur_lower`. The ground decision is undefined behaviour. Replace with a smoothstep on `direction.z`.
- `horizon_height` is applied to Arnold `.y` before the `y_is_up` swap, so it is wrong when `y_is_up` is false. The ground adds a second offset (`v_h.z = hh·-0.1`).
- `getVolumeEffects`: the "unoccluded" colour is a magenta debug value (1, 0, 0.2). `sg->Ci` is never used.
- No real aerial perspective. The closed forms are used with s = 32000 for every pixel. A surface hit should use s = `Rl`: L = L_surface·Fex(Rl) + Lin(Rl). `atmos_depth` is ignored.
- Default output is not Preetham. Turn these off before any test against reference values:
  - `node_update` (`arnoldSkyShader.cpp:164-165`) rewrites `multB`. At defaults, Perez B is scaled by 0.55-1.0 with sun height.
  - Perez C is scaled by `multD`. `mult_c` has no effect.
  - `attenuateSky` lerps extinction towards 1 for sun zenith 60-90°.
  - `dusk` fades the sky for sun zenith 90-96°.
- Smaller:
  - `maxTurbidity` is really the minimum.
  - No upper turbidity clamp (Preetham is valid to about 10-12).
  - `clamp(multD)` mutates the member. When `mult_d` is outside [0, 2], every update rebuilds the sky.
  - `far_clip` is a static shared by all nodes.
  - `AiUniverseGetCamera` can return null.
  - `isVisible` compares `sg->Rt == AI_RAY_ALL_DIFFUSE` and `AI_RAY_ALL_SPECULAR`. These are masks, so the diffuse and glossy toggles never match.

## 3. SSE2
- Confirmed bug (`spectralData.h:7-17`, `spectrum.h:141-143`). Each `PreethamData` member is `ALIGN(16)`. Compiled: sizeof 128, S0 at offset 80, S1 at 96. Each member is a double plus 8 bytes of padding. `_mm_load_pd(&spectral_data[i].S0)` loads {S0[i], padding}. Every odd band of `spectrum(x,y)` is zero, so the sky spectrum is a 20 nm comb.
  - Sky colour: mild shift, because `getSkySpectralRadiance` renormalises by Y.
  - AP: large effect. The sky part of the odd bands of S0 is zero, so those bands get only the sun term.
  - Fix: struct of arrays (`double S0[38]`) or scalar loops.
- `operator-(spectrum)`, `operator-(Real)` and `operator/(spectrum)` compute `in - this` and `in / this`. Latent, unused. `operator-(void)` is not const.
- `spectrum(Real[])` uses an aligned load, so an unaligned caller crashes. The current caller is a local `double[38]`, which x86-64 aligns to 16.
- Little benefit. The loops are 19 iterations of 2 doubles. GCC -O3 auto-vectorises plain loops wider with AVX2. The hand SSE2 locks you to 128-bit and to x86. The `USE_AVX` path is dead.
- Wasted work: default construction zeroes memory, and each `operator+` builds a zeroed temporary and overwrites it. `exponentSpectrum`, `toXYZ` and `lerp` are scalar. Prefer plain loops over `alignas(64) float[40]` with `#pragma omp simd`.
- The real cost is elsewhere:
  - `ap_CalculateS0` recomputes `getSkySpectralRadiance` × `ap_ExtinctionFactor` for the same 400 directions for each of 441 table entries. These do not depend on the view direction. Hoisting gives about 400× faster `reinit`.
  - `shader_evaluate` spends about 300 `EXP` calls per sample. Sky and AP depend only on direction, so bake them.
- `Helper1` subtracts near-equal integrals with 1/K³ terms. This is probably why double is needed. Validate `t2exp` against `std::exp` at low K.

## 4. Is aerial perspective applicable to a path tracer?
Yes, as a post-hit modifier, and only if scenes are large.
- The closed forms depend on direction, start height and distance. Bake Fex and Lin into small tables at setup. Apply L·Fex(s) + Lin(s) to a hit with s = segment length. This is deterministic, needs no random draws and no heap on the sample path.
- It is not a participating medium. The density is exponential, so a homogeneous medium gives wrong colour at long range.
- Inscatter is unoccluded. Terrain does not shadow it (no light shafts).
- It changes what the camera sees. It does not light the scene.
- Range depends on turbidity. At T = 2, Fex(500 m) ≈ 0.95 at 550 nm, nearly invisible. At the default T = 5, β_p is about 4× larger, so Fex(500 m) ≈ 0.85 and Lin is visible.
- Convert scene distance to metres with the scene's meters-per-unit.
- Apply to hits only. Sky miss rays already include it.
- h0 sets the scope:
  - A table at camera height covers camera segments only.
  - Secondary segments (a lake that reflects a far mountain) need h0 as a table dimension, or no AP.
  - Decide this before you build the table.
- AOVs: Fex scales every lobe, and Lin is additive. To keep beauty = Σ lobe planes, apply Fex to the lobe planes and give Lin its own plane. Data planes stay unaffected.

## 5. Integrating into a renderer
1. Make the core free of host-SDK types. Plain float3/double types, no `omp`, no SSE2. Do the physics once at setup, with nothing allocated on the sample path.
2. Sky as a baked environment map. Evaluate into a working-space float image at prepare time (about 256×128, finer near the horizon). Feed it to the renderer's environment light so importance sampling and MIS apply.
3. Sun as a distant light, 0.533° diameter, irradiance from `sunLightCol`. Keep the disc out of the baked map to avoid double counting and a spike in the environment sampler.
4. AP in the integrator as a post-hit modifier, behind a switch.
5. Parameters: `turbidity`, sun direction, ground albedo, `multiplier`, `camera_height`, `overcast`. Drop `candela_conversion`, `saturation`, `tonemap`, the Perez multipliers and the dusk hacks.
6. Fix before porting, in this order:
   1. The `spectrum` layout bug.
   2. The Rayleigh cos(cosψ) bug.
   3. The double sky attenuation.
   4. The Perez denominator and the overcast law.
   5. The 2° observer.
   6. The `ap_CalculateS0` grid (full 2π) and units.
   7. Remove the B/C rewrite, `attenuateSky` and `dusk`.
   8. The ground undefined behaviour.
7. Model choice: Preetham is unreliable below T ≈ 2 and with the sun near the horizon (Zotti, Wilkie, Purgathofer 2007). Hosek-Wilkie is the usual replacement for the sky dome. The AP closed forms are Preetham-specific and can stay next to it.
8. Tests, with all hacks off: baked sky luminance at a few angles against Preetham's published values, and sun illuminance near 1e5 lux at the zenith for T = 2.
