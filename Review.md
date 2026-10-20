# aaPhysicalSky review

Read from source. Only the struct layout in section 3 was compiled. Nothing was rendered or compared with reference values.

## 1. Units

### Correct
- Sun table in `spectralData.h` is W·cm⁻²·sr⁻¹·µm⁻¹ radiance. The chain ×1e4, ÷1e3, ×Ω (`physicalSunClass.h:146-148`) gives W/m²/nm irradiance.
- Sun ×683.002 gives lux. `sunDiscCol` = lux/Ω is about 1.6e9 cd/m².
- Perez zenith Y: kcd/m² ×1000 is correct.
- β_m is near 1.3e-5 m⁻¹ at 550 nm. β_p is plausible. The AP maths is in SI metres.
- Ground colour ×`candelaConversion` ÷`multiplier` cancels correctly.

### Wrong or mixed
1. `ap_CalculateS0` (`physicalSkyClass.h:688-690`):
   - `delPhi = RadsToDegs(π/nPhi)` is 9°, not 18°. 20 steps cover only 180° of azimuth.
   - `delThetaPhi` is in degrees², about 3283× too big.
   - The comments say 9° and 18°, but `delTheta` is 4.5°.
2. Sun and sky inscatter weights are not photometric-consistent. Sun (W/m²/nm) is ×100. Sky (cd/m²) is ×1e-4. The sun term should be ×683.002. The sky term should be ×1 with radians and the full 2π. Both factors are about 0.15× of correct, so the ratio is nearly right and the total is about 6× too weak. Not checked numerically.
3. `toXYZ` (`spectrum.h:373-379`) uses the 10° CMF (columns 4-6). Preetham x,y and 683 lm/W are 2° observer. Use columns 1-3.
4. No chromatic adaptation. `XYZtoACEScg` uses the D60-referenced AP1 matrix. The sky is D65. `clampRGB` hard-zeroes negatives instead of gamut-mapping.
5. `h0`, `s = 32000`, `Rl` and `max_visibility_distance` are treated as metres. No conversion from scene units.
6. `h0 = -camera_height` makes air denser with altitude, since e^(−αh0) = e^(+α·height). Harmless at 1 m, wrong at 1 km.
7. Below-horizon inscatter integrates underground to s = 32 km with signed cosθ. Extinction uses `fabs`. e^(α₁·32000) is about 4e11 for Mie. Check for NaN or huge values on ground rays.
8. `addSun` uses `half_dot = dot*0.5` in `rescale`. Limb darkening and blur are about constant (blur ≈ 0.78, limb ≈ 0.95). The disc is about 26% dimmer than the light's irradiance, with no limb darkening. `sun_disk_scale` is never used.

## 2. Model bugs
- `perez()` line 400: the denominator uses `cos2gamma`. Preetham needs cos²θs. Correct only at the zenith.
- CIE overcast (`:429`): Y = Yz(1 + cosθ/3). The CIE law is Yz(1 + 2cosθ)/3.
- `arnoldSkyShader.cpp:258`: `_xyY.y = overcastWhitePoint_D65_x` sets y to x.
- `addGround:486`: `AiV3Cross(v_cross, v_camDir)` reads an uninitialised `v_cross` and assigns the result to `v_h_blur_lower`. The ground decision is undefined behaviour. Replace with a smoothstep on `direction.z`.
- `horizon_height` is applied to Arnold `.y` before the `y_is_up` swap, so it is wrong when `y_is_up` is false. The ground adds a second offset (`v_h.z = hh·-0.1`).
- `getVolumeEffects`: the "unoccluded" colour is a magenta debug value (1, 0, 0.2). `sg->Ci` is never used.
- No real aerial perspective. The closed forms are used with s = 32000 for every pixel. A surface hit should use s = `Rl`: L = L_surface·Fex(Rl) + Lin(Rl). `atmos_depth` is ignored.
- Smaller:
  - `maxTurbidity` is really the minimum.
  - No upper turbidity clamp (Preetham is valid to about 10-12).
  - `clamp(multD)` mutates the member and forces a rebuild on every update.
  - `far_clip` is a static shared by all nodes.
  - `AiUniverseGetCamera` can return null.

## 3. SSE2
- Confirmed bug (`spectrum.h:75-77`). Each `PreethamData` member is `ALIGN(16)`. Compiled: sizeof 128, S0 at offset 80, S1 at 96, so the stride is 16 bytes. `_mm_load_pd(&spectral_data[i].S0)` loads {S0[i], padding}. Every odd band of `spectrum(x,y)` is zero, so the sky spectrum is a 20 nm comb and chromaticity is skewed. Fix: struct of arrays (`double S0[38]`) or scalar loops.
- `operator-(spectrum)`, `operator-(Real)` and `operator/(spectrum)` compute `in - this` and `in / this`. Latent, unused. `operator-(void)` is not const.
- `spectrum(Real[])` uses an aligned load, so an unaligned caller crashes.
- Little benefit. The loops are 19 iterations of 2 doubles. GCC -O3 auto-vectorises plain loops wider with AVX2. The hand SSE2 locks you to 128-bit and to x86. The `USE_AVX` path is dead.
- Wasted work: default construction zeroes memory, and each `operator+` builds a zeroed temporary and overwrites it. `exponentSpectrum`, `toXYZ` and `lerp` are scalar. Prefer plain loops over `alignas(64) float[40]` with `#pragma omp simd`.
- The real cost is elsewhere:
  - `ap_CalculateS0` recomputes `getSkySpectralRadiance` × `ap_ExtinctionFactor` for the same 400 directions for each of 441 table entries. Hoisting gives about 400× faster `reinit`.
  - `shader_evaluate` spends 400+ `EXP` calls per sample. Sky and AP depend only on direction, so bake them.
- `Helper1` subtracts near-equal integrals with 1/K³ terms. This is probably why double is needed. Validate `t2exp` against `std::exp` at low K.

## 4. Is aerial perspective applicable to a path tracer?
Yes, as a post-hit modifier, and only if scenes are large.
- The closed forms depend on direction, camera height and distance. Bake Fex and Lin into small tables at setup. Apply L·Fex(s) + Lin(s) to a hit with s = ray length. This is deterministic, needs no random draws and no heap on the sample path.
- It is not a participating medium. The density is exponential, so a homogeneous medium gives wrong colour at long range.
- It changes what the camera sees. It does not light the scene.
- It pays off on mountains, ocean horizons and very large scenes. Under about 1 km it is nearly invisible (Fex at 500 m ≈ 0.95 at 550 nm).
- Convert scene distance to metres with the scene's meters-per-unit.
- Apply to hits only. Sky miss rays already include it.
- Use camera height for h0.
- Apply to beauty only. Leave AOVs and data planes unaffected.

## 5. Integrating into a renderer
1. Make the core free of host-SDK types. Plain float3/double types, no `omp`, no SSE2. Do the physics once at setup, with nothing allocated on the sample path.
2. Sky as a baked environment map. Evaluate into a working-space float image at prepare time (about 256×128, finer near the horizon). Feed it to the renderer's environment light so importance sampling and MIS apply.
3. Sun as a distant light, 0.533° diameter, irradiance from `sunLightCol`. Keep the disc out of the baked map to avoid double counting and a spike in the environment sampler.
4. AP in the integrator as a post-hit modifier, behind a switch.
5. Parameters: `turbidity`, sun direction, ground albedo, `multiplier`, `camera_height`, `overcast`. Drop `candela_conversion`, `saturation`, `tonemap` and the dusk hacks.
6. Fix before porting: the Perez denominator, overcast law, observer, S0 constants, the `spectrum` layout bug and the ground UB.
7. Tests: baked sky luminance at a few angles against Preetham's published values, and sun illuminance near 1e5 lux at the zenith for T = 2.
