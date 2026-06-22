---
tags: [antenna_design, design_notes]
aliases: [Series aperture-coupled patch antenna array design]
---

# Aperture-coupled patch antenna design

## Theory of operation and electrical model

The aperture-coupled topology isolates the radiating element from the feed network using separate dielectric layers divided by a common ground plane. Energy transfers magnetically through a ground-plane aperture (slot).

- **Equivalent circuit:** The slot acts as a transformer (patch’s coupled impedance) and a parallel LC circuit (slot resonance and fringing fields) in series with the microstrip feed.
- **Smith chart locus:** The impedance locus typically forms a circle. Tuning aims to size and shift this circle to pass through the $50\,\mathrm{\Omega}$ centre ($Z=1$).

## Substrate selection: theory vs practice

A key advantage of this architecture is the independent optimization of feed and antenna substrates. However, theoretical ideals often clash with the physical realities of dense, high-frequency arrays.

- **Antenna substrate:** Requires a low-loss, electrically thick material ($h \geq \lambda/10$) with low relative permittivity ($\epsilon_\mathrm{r} < 3$). This encourages outward radiation at the patch edges, widens impedance bandwidth, and reduces surface-wave propagation (crucial for MIMO channel isolation).
- **Feed substrate:** Theory suggests an electrically thin substrate ($h < \lambda/50$) with high permittivity ($\epsilon_\mathrm{r} > 5$) to confine fields into a guided wave, maximize aperture coupling, and prevent spurious feed radiation.
- **The trace-width reality check:** At lower microwave bands (e.g. 9.7 GHz), high-$\epsilon_\mathrm{r}$ feed substrates (e.g. RO3006, $\epsilon_\mathrm{r}=6.5$) yield practical $50\,\mathrm{\Omega}$ trace widths (e.g. 0.38 mm). When scaling to 79 GHz automotive bands, however, $\epsilon_\mathrm{r} \geq 5$ forces microscopic trace widths, compounding conductor losses and etching tolerances. Consequently, millimetre-wave arrays practically require an $\epsilon_\mathrm{r}$ between 3 and 4 for the feed layer.
- **The surface-wave reality check:** At X-band (e.g. 9.7 GHz on $787\,\mathrm{\mu{}m}$ RO3006), the substrate electrical thickness ($\approx 0.025\lambda_0$) inherently excites the dominant $\text{TM}_0$ surface wave. In large MIMO arrays (e.g. $8 \times 4$), this propagates and couples into adjacent elements, causing ‘scan blindness’ – reflection spikes at specific scan angles. When scaling to 79 GHz, the $h/\lambda_0$ ratio spikes further. Failing to drastically thin the substrate allows higher-order modes (e.g. $\text{TE}_1$) to proliferate, leading to severe power loss and spatial decorrelation. Transitioning to a substrate-integrated waveguide (SIW) cavity-backed topology is the ultimate solution for preserving virtual aperture integrity in millimetre-wave MIMO.

## Parametric tuning guide

When tuning in a full-wave solver, it is essential to map each physical ‘knob’ to its electrical response.

- **Patch length:** The primary determinant of operating frequency (inversely proportional).
- **Patch width:** Primarily affects resonant resistance (wider patch = lower resistance).
- **Slot length:** The primary control for coupling level. Increasing it expands the diameter of the impedance locus on the Smith chart. The slot must remain strictly sub-resonant (ideally $\approx \lambda_\mathrm{g}/10$) to suppress excessive back-radiation. If $\lambda_\mathrm{g}/2$ is needed to achieve a match, the feed substrate is likely too thick or its permittivity too low.
- **Slot width:** Has a minor effect on coupling. Maintain a standard length-to-width ratio of $\approx 10:1$. Modifying the shape to a ‘dog-bone’ or ‘bow-tie’ can increase magnetic polarizability without expanding the overall footprint.
- **Feed stub length:** The open-ended microstrip extending past the slot tunes out excess slot reactance. It is theoretically slightly less than $\lambda_\mathrm{g}/4$.[^1] Shortening the stub rotates the impedance locus downwards (capacitive direction) on the Smith chart.
- **Alignment offsets:** Maximum coupling requires the patch to be perfectly centred over the slot, and the feedline to be centred and perpendicular to the slot. Skew degrades coupling, cross-polarization isolation, and broadside gain.

## Personal insights and tuning strategy

- **Over-coupling for bandwidth:** Instead of forcing a high-Q resonance exactly at $Z=1$, it is beneficial to ‘overshoot’ the coupling (by lengthening the slot). This causes the impedance locus to intersect the $R=1$ circle at two points symmetrically around the centre, significantly broadening the operational bandwidth.
- **Reactance compensation:** Changing the slot length inherently alters its reactance. Always counteract this by re-tuning the stub length.
- **Resonance pulling:** The slot length pulls the patch resonance down. Slight patch length variations will be necessary to re-centre the frequency. _Caution:_ Even a 0.1 mm change to the patch length at X-band drastically alters the resonance – tune it gently.
- **Stub length:** While classical 1D theory dictates a $\lambda_\mathrm{g}/4$ stub, the actual length can be merely 35–40% of $\lambda_\mathrm{g}/4$. This severe discrepancy stems from localized dielectric loading (fields fringing into the high-$\epsilon_\mathrm{r}$ patch substrate drastically lower the local $\lambda_\mathrm{g}$ over the aperture), fringing capacitance at the stub’s open end (acting as a virtual length extension), and the distributed nature of the coupling slot in 3D space blurring the exact reference plane.

---

# Series-fed arrays

Transitioning from a single element to a series-fed array fundamentally alters the electromagnetic problem from a localized load to a cascaded microwave network.

## Array termination and topology

Intermediate elements in a series array do not possess individual tuning stubs; the feedline simply crosses the slot and continues. The termination after the final element dictates the fundamental behaviour of the array.

- **Resonant arrays:** The line is terminated with an open-ended tuning stub after the final element, establishing a standing wave along the feedline. This is a highly radiation-efficient solution, yet fundamentally limited (typically 1–2%) in its operational bandwidth; as frequency shifts, the standing wave nodes drift away from the slots, degrading the amplitude distribution.
- **Travelling-wave arrays:** The line is terminated with a matched $50\,\mathrm{\Omega}$ load, absorbing residual power and preventing standing waves. This offers much wider impedance bandwidth, but the main beam inherently squints (tilts off broadside) over frequency – a phenomenon highly detrimental to FMCW radar spatial resolution.
	- [I] The beam squint is – in theory – completely eliminated in centre-fed travelling-wave arrays due to symmetric squints in each direction effectively cancelling out.
- **Polarization purity advantage:** A major advantage of utilizing an aperture-coupled topology for a resonant series array is the absolute suppression of feed network radiation. Confining the cascaded feed behind a solid ground plane yields the exceptional cross-polarization discrimination (XPD) essential for polarimetric processing.

### The paradox of broadside travelling-wave arrays

Series-fed travelling-wave arrays terminated in a (virtual) matched load present a fundamental design paradox when configured to radiate broadside:

- **Broadside phasing condition:** for in-phase broadside radiation, the electrical phase delay between adjacent elements must be a multiple of $360^\circ$. This dictates that the physical inter-element feedline length must be equal to an integer multiple of the guided wavelength ($L = n\lambda_\mathrm{g}$).
- **Constructive reflection buildup:** under this exact condition, any localized reflections generated by the individual elements travel back towards the input port. Because the round-trip distance between successive elements is exactly a multiple of $2n\lambda_\mathrm{g}$, the reflections from all elements add up constructively at the input.
- **The paradox:** The very condition that aligns the main beam to broadside maximizes the input reflection, generating a massive standing wave and destroying the travelling-wave condition. Conversely, squinting the beam (out-of-phase feeding) introduces a phase mismatch that causes inter-element reflections to interfere destructively at the input, yielding excellent return loss but pointing the main beam off-broadside.

To resolve this conflict, the following architecture paths were established.

**1. Travelling-wave mode with intentional squint**
- Keep the travelling-wave architecture to preserve bandwidth and stability.
- Intentionally apply a progressive phase shift to squint the beam slightly off-broadside in elevation.
	- [i] A slight squint in elevation should be easy to eliminate during mechanical mounting.
- This progressive phase shift prevents constructive reflection buildup and restores the input match.

**2. Resonant array made with input matching network**
- If a strict broadside beam is required, switch to a standing-wave (resonant) topology.
- Remove the progressive phase shift across elements.
- Integrate a dedicated feed matching network (e.g. a localized microstrip transformer) at the array input to match the resulting real (but non-$50\,\Omega$) input impedance back to the system impedance.

---

## Amplitude tapering and coupling

If identical slots are placed uniformly along a series feed, the exponential decay of incident power raises side-lobe levels and skews the beam. Achieving a specific amplitude taper in a travelling-wave array requires a distinct energy coupling distribution at each element to compensate for this progressive leakage.

Unlike resonant arrays that maintain uniform energy levels via standing waves, a travelling-wave array requires the first element to be strongly under-coupled and the terminal element to be strongly over-coupled. This gradient is especially critical when omitting a traditional physical load in favour of a virtual load to absorb residual energy.

### Radiated-to-available power ratio (RAPR) and synthesis strategy

This required gradient is quantified by the radiated-to-available power ratio (RAPR), defined as the ratio of the power radiated by an element to the power available to it. Following the design procedure outlined by Kang et al. ,[^2] the array synthesis involves calculating the desired RAPR profile to achieve the target sidelobe level (SLL) and mapping it to physical element dimensions. To cover the required dynamic range of RAPR (e.g., from <1% near the feed point to >70% at the terminal end), the array must employ different element topologies (‘types’) to independently manage both high and low coupling states without introducing reflections.

### Results using simple rectangular slots only

There is no requirement to simulate further parameter combinations for these specific variables, as the sweep fully captures the peak resonance. If the tapering profile subsequently demands an element with $>24\%$ coupling, sweeping the same variables will yield minimal improvement. Instead, fundamental parameters must be modified (e.g. widening the slot width, increasing the substrate thickness) or a larger fraction of power must be dissipated in the matched load to artificially lower the required coupling percentages.
- **Low-coupling elements:** Weak coupling levels extend down to 0.45%. This is excellent for heavily tapered profiles (such as high sidelobe-suppression Chebyshev or Taylor tapers) where the edge elements ($1$ and $N$) must extract only a minute fraction of the power.
- **High-coupling elements:** The peak coupling terminates at 24.4%. For a travelling-wave array, this acts as the physical upper bound.

From an electromagnetism and antenna array design perspective, whether a maximum coupling of 24.4% is sufficient depends heavily on three factors: the size of the array ($N$), the feeding topology (travelling-wave versus resonant), and the target radiation efficiency.

> [!warning] For a small array ($N = 7$), a limit of 24% is insufficient, necessitating a slot modification.

### Physical principles and efficiency degradation

In a series-fed travelling-wave array, power propagates down the transmission line, and each successive element extracts a percentage of the power remaining on the line. Because the array contains only seven elements, fewer elements are available to radiate the total input power. If a high radiation efficiency is required (e.g. 90%, where only 10% is dissipated in the terminal load), the individual elements must exhibit a high coupling coefficient.

This constraint is exacerbated significantly by amplitude tapering (such as Taylor or Chebyshev distributions):
- In a tapered array, the central elements (e.g. elements 3, 4, and 5) are responsible for radiating the vast majority of the power.
- As the guided wave reaches the centre of the array, the middle elements must mathematically couple very strongly to achieve the required radiated amplitude.
- If the physical slot is constrained to a maximum power extraction of 24.4%, the mathematical taper ratio can only be satisfied by flooding the entire feedline with excess power. The centre element achieves its 24.4% maximum, whilst the substantial leftover power propagates down the remainder of the line and is dissipated as heat in the terminal matched load.

> [!example] Mathematical illustration
> If a heavily tapered profile is applied to an $N = 7$ array whilst restricting the maximum coupling of any single element to $\le 24.4\%$, approximately 40% to 50% of the total input power is mathematically forced into the terminating load. Consequently, the antenna radiation efficiency drops to approximately 50%.
> - **Resonant arrays:** If designing a centre-fed standing-wave (resonant) array, 24% is entirely sufficient. Power is not extracted in a single pass; instead, it reflects and establishes a standing wave.
> - **Large arrays ($N > 20$):** In long travelling-wave arrays, the power is distributed across many elements, meaning no single element is required to couple strongly. For $N = 30$, a maximum coupling of 24% is fully adequate.

---

## Phase matching and physical spacing

To generate a standard automotive fan beam (broad in elevation, highly directive in azimuth) that points broadside without grating lobes, the physical and electrical spacing must be carefully decoupled.

- **Physical spacing:** Patch phase centres must remain physically spaced at $\approx 0.5\lambda_0$ (and strictly $< 0.7\lambda_0$) to avoid grating lobes.
- **Electrical phase:** The fields at each patch must be strictly in phase to radiate broadside, requiring a total electrical phase shift of precisely $360^\circ$ between adjacent elements.
- **U-bend compensation:** Because the guided wavelength ($\lambda_\mathrm{g}$) on the feed substrate is physically shorter than the free-space $\lambda_0$, fitting a full $360^\circ$ electrical delay into a smaller $0.5\lambda_0$ physical gap requires physically lengthening the microstrip line. This can be achieved by routing the line through U-bends or meanders between the slots.

### Travelling-wave phase compensation

In travelling-wave series arrays, tuning elements to different coupling levels inherently shifts their transmission and radiation phases. Based on the synthesis logic by Kang et al.,[^2] realizing a highly accurate beam requires tracking the explicit **radiation phase ($\angle E_{\mathrm{rad}}$)** of each element.

To maintain the correct progressive phase ($\Phi_{\mathrm{prog}}$) across the array, the precise electrical length ($\Delta\Phi_{\mathrm{feed}}$) required for the connecting feed line between element $i$ and element $i+1$ is calculated as:
$$
\Delta\Phi_{\mathrm{feed}} = \angle S_{21,i} + \angle E_{\mathrm{rad},i+1} - \angle E_{\mathrm{rad},i} - \Phi_{\mathrm{prog}}
$$
This elegantly compensates for individual radiation phase variations (which can otherwise squint the beam or degrade the sidelobe level) by appropriately altering the physical meander length of the feed line segment separating each unique element pair.

### Array scaling (element count and physical spacing)

When scaling the array to maximize aperture directivity, two parameters dictate performance bounds: element count ($N$) and physical inter-element spacing ($d$).

- **Element count ($N$):** Increasing the patch count inherently smooths the required amplitude taper, significantly reducing the maximum necessary coupling coefficient ($C_i$) at the central elements. This allows the array to radiate a higher percentage of the total power before reaching the terminal load, increasing overall radiation efficiency. However, in series-fed topologies, power must physically traverse the lossy microstrip feedline. As $N$ increases, the $+10\log(N)$ directivity gain is eventually outpaced by cumulative dielectric and conductor insertion losses. For X-band arrays on standard high-frequency substrates (e.g. RO3006), realized gain typically saturates between 12 and 16 elements.
- **Physical spacing ($d$):** Expanding the spacing to approach the theoretical $\lambda_0$ limit increases directivity per element but is highly detrimental for series-fed travelling-wave arrays. First, travelling-wave arrays inherently squint over frequency. Setting $d \approx \lambda_0$ guarantees that even a minor frequency-induced squint will immediately drag a grating lobe into visible space. Second, maintaining the required broadside electrical phase across a wider physical gap demands a much longer interconnecting feedline (e.g. stretching from $1\lambda_g$ to $2\lambda_g$). This longer line possesses a steeper phase slope, drastically worsening the beam’s frequency dispersion (squint rate). Maintaining a spacing of $0.5\lambda_0$ to $0.7\lambda_0$ safely mitigates both issues.

---

## Characterizing aperture coupling for series-fed arrays

Synthesizing a specific amplitude taper (such as a Taylor or Chebyshev distribution) requires precise quantification of each element’s coupling strength. This is achieved by placing a single slot-coupled patch in a two-port unit-cell environment. By sweeping the slot geometry and feedline parameters, the power coupled to the patch ($|C|^2$) is extracted from the scattering parameters of the lossless network:

$$
|C|^2 = 1 - |S_{11}|^2 - |S_{21}|^2
$$

> [!example]+ Technical considerations for extraction
> - **De-embedding:** reference planes in the full-wave solver must be de-embedded to the centre of the slot to remove phase rotation and dielectric loss from the access lines.
> - **Sub-resonant behaviour:** the slot must remain strictly sub-resonant. High $|S_{11}|$ values introduce impedance mismatches and unwanted standing waves, disrupting the travelling-wave amplitude distribution.
> - **Loading effects:** mutual coupling between adjacent patches can alter individual input impedances. High-fidelity designs require a second-pass optimization to account for these loading effects.

---

## High-coupling elements & impedance matching

A fundamental challenge during array synthesis is that the amplitude taper dictates the radiation profile, but the required coupling elements disrupt the feed line. Each coupling structure acts as a complex series impedance, $Z_\mathrm{slot} = R_\mathrm{slot} + jX_\mathrm{slot}$. Even if an element radiates perfectly, the uncompensated reactance and the discrete introduction of resistance cause small reflections. Across multiple elements, these reflections accumulate, breaking the travelling-wave condition and generating a standing wave.

To achieve the high radiated-to-available power ratio (RAPR > 70%) required for the central elements of a heavily tapered array, the slot must extract a massive amount of power.

> [!failure]+ Infeasibility of the feedline gap in aperture-coupled design
> In direct series-fed microstrip arrays (e.g., Kang et al., 2020), a physical gap in the feedline over the patch acts as a series capacitor, forcing strong coupling (the ‘through-element’). However, applying this gap directly over an _aperture-coupled_ slot completely fails. Because aperture coupling relies on the continuous magnetic field of the travelling wave looping through the slot, a gap severs the longitudinal current path, collapsing the magnetic field. This results in massive reflections ($S_{11} \approx 0\,\mathrm{dB}$) and essentially zero coupling.

### Two-part local thinning and impedance transformation

To achieve high RAPR in an aperture-coupled design without introducing severe reflections, a two-part feedline modification is implemented. We treat the entire coupling region (the slot and its localized trace) as a single ‘black box’ series impedance block, because a transverse slot interrupting the longitudinal ground return current fundamentally manifests as a series-loaded impedance to the travelling wave.

1. **Maximized magnetic coupling.** A modified hourglass or bow-tie slot increases magnetic polarizability. Simultaneously, the microstrip line is drastically thinned (e.g. down to 0.1 mm) strictly in the region immediately traversing the slot. This high-impedance section forces a very high current density, blasting magnetic flux into the aperture. By extracting the complex impedance of this _entire_ composite structure (slot + thinned trace) from simulation, any localized impedance jumps within the coupler are fully accounted for within a single macroscopic series impedance, $Z_\mathrm{slot}$.
2. **Quarter-wave transformer.** If an element extracts 70% of the power, the equivalent series resistance $R_\mathrm{s}$ of the composite coupler is roughly $2.33 Z_0$ ($\approx 116\,\mathrm\Omega$). The total input impedance is $Z_\mathrm{in} = R_\mathrm{s} + Z_0 \approx 166\,\mathrm\Omega$. Connecting this directly to the incoming $50\,\mathrm\Omega$ line creates severe reflection. To solve this, an integrated transmission line matching network (a medium-thin trace acting as a transformer) is placed immediately upstream to _precede_ the slot. _Crucially, there is no transformer after the slot._ The line must immediately return to the nominal $Z_0$ width, as the continuing line physically _is_ the $50\,\mathrm\Omega$ load in the equivalent circuit.

### Slot impedance extraction and transformer analytics

To grasp the complex slot impedance we are compensating for, we extract it from a two-port unit cell. Assuming the composite coupler is a series impedance, the total load impedance is $Z_L = Z_\mathrm{slot} + Z_0$.[^3] Substituting this into the reflection coefficient definition yields $S_{11} = Z_\mathrm{slot} / (Z_\mathrm{slot} + 2Z_0)$. Rearranging for $Z_\mathrm{slot}$ yields
$$
Z_\mathrm{slot} = 2Z_0 \frac{S_{11}}{1 - S_{11}}
$$

To analytically estimate the matching section, we treat the preceding microstrip transformer as a transmission line of impedance $Z_1$ and electrical length $\theta$. Equating its input impedance – when terminated in the disrupted load $Z_L = Z_0 + R_\mathrm{slot} + jX_\mathrm{slot}$ – to the main line impedance $Z_0$, we get
$$
Z_0 = Z_1 \frac{(Z_0 + R_\mathrm{slot} + jX_\mathrm{slot}) + jZ_1 \tan(\theta)}{Z_1 + j(Z_0 + R_\mathrm{slot} + jX_\mathrm{slot}) \tan(\theta)}
$$

Separating this into real and imaginary components and solving the real part for $\tan(\theta)$ provides the required electrical length (dictating the physical transformer length):
$$
\tan(\theta) = -\frac{Z_1 R_\mathrm{slot}}{Z_0 X_\mathrm{slot}}
$$

Substituting this expression into the imaginary part equation and rearranging for $Z_1$ yields the required transformer impedance (dictating the required transformer trace width)
$$
Z_1 = \sqrt{Z_0^2 + \frac{Z_0}{R_\mathrm{slot}}(R_\mathrm{slot}^2 + X_\mathrm{slot}^2)}
$$

> [!note]+ Note on the transformer equation
> While a standard quarter-wave transformer matching a terminal load $R_L$ uses the geometric mean ($Z_1 = \sqrt{Z_0 R_L}$), our architecture places the slot in _series_ with the continuing feedline. The load is therefore an additive series combination ($Z_L = Z_0 + Z_\mathrm{slot}$). The resulting derivation perfectly accounts for this cascaded series load, yielding the elegantly expanded square root of addition shown above, which mathematically collapses back to the classical geometric mean if the line were simply terminated ($Z_0 \rightarrow 0$).

> [!example]+ Parametric mapping strategy
> For automated array synthesis, we must cleanly map the parameter space so the Python optimizer can interpolate the exact dimensions needed to hit a target ACC while maintaining $\text{Im}(Z) \approx 0$. Because exploring a 5D parameter space (slot dimensions + transformer dimensions) is computationally prohibitive, the search space is divided into three logical families:
>
> 1. **High-coupling regime (up to -2 dB ACC):** Utilizes the full impedance transformer and heavily thinned coupling section. The trace thinning over the slot is fine-tuned to dial the extreme coupling, while a coarse sweep over the transformer length ensures optimal $S_{11}$ tracking.
> 2. **Vanilla hourglass (up to -4 dB ACC):** For mid-range elements, the transformer and feedline thinning are completely disabled (returned to nominal $Z_0$). The coupling is dialed purely by sweeping the hourglass flare width (`W2`).
> 3. **Simple slot (-6 dB ACC):** For the weakest elements near the feed point, the hourglass shape is neutralized into a standard rectangular slot (`W1 = W2`), and the coupling is dialed by sweeping the slot width.
> 
> In all sequences, the `slotLength` is finely swept. This is the most critical design rule: fine resolution on `slotLength` guarantees the optimizer will always capture the exact frequency crossing where $\text{Im}(Z_\mathrm{slot}) \approx 0$, regardless of how the other parameters distort the local reactance.

## Results using modified slot shapes and impedance transformers

To overcome the coupling limits of the standard rectangular slot, an $N = 8$ series-fed travelling-wave array was synthesized using a ‘modified-hourglass’ slot profile – which significantly increases magnetic polarizability – combined with localized microstrip feedline thinning (down to $0.1\,\mathrm{mm}$) and integrated quarter-wave impedance transformers. Full-wave electromagnetic simulation sweeps of the progressive phase shift ($\Phi_{\mathrm{prog}}$) were conducted at $9.7\,\mathrm{GHz}$ to characterize the input match, bandwidth, and radiation pattern properties.

### Return loss and bandwidth collapse

A key challenge in travelling-wave series arrays is the compensation of the series-loaded slot impedances. The progressive phase sweep highlights a severe return loss sensitivity:

| Progressive phase ($\Phi_{\mathrm{prog}}$) | $S_{11}$ at $9.7\,\mathrm{GHz}$ (dB) | $f_{\mathrm{min}}$ (GHz) | $f_{\mathrm{max}}$ (GHz) | $-10\,\mathrm{dB}$ return loss bandwidth (%) |
| :----------------------------------------- | :---------------------------------- | :---------------------- | :---------------------- | :------------------------------------------ |
| $-15.0^\circ$                              | $-12.40$                            | $9.682$                 | $9.762$                 | $0.82$                                      |
| $-10.0^\circ$                              | $-7.77$                             | –                       | –                       | $0.00$                                      |
| $-5.0^\circ$                               | $-5.51$                             | –                       | –                       | $0.00$                                      |
| $0.0^\circ$                                | $-4.36$                             | –                       | –                       | $0.00$                                      |
| $5.0^\circ$                                | $-3.80$                             | –                       | –                       | $0.00$                                      |
| $10.0^\circ$                               | $-3.93$                             | –                       | –                       | $0.00$                                      |
| $15.0^\circ$                               | $-3.01$                             | –                       | –                       | $0.00$                                      |

At the designed progressive phase shift of $-15.0^\circ$, the array achieves a solid input match of $-12.40\,\mathrm{dB}$ at $9.7\,\mathrm{GHz}$. However, the operational impedance bandwidth is extremely narrow, measuring only $0.82\%$ ($80\,\mathrm{MHz}$ absolute bandwidth). For all other progressive phase angles, the input match collapses entirely ($S_{11} > -10\,\mathrm{dB}$).

This behaviour is very likely due to the cascaded quarter-wave matching transformers. Matching high-impedance series slots to the main feedline using multiple cascaded high-ratio $\lambda_g/4$ sections introduces high frequency sensitivity: slight deviations from the centre frequency rotate the impedance locus away from the match point, causing reflections to accumulate constructively at the input.

### Radiation performance

The progressive phase shift controls the beam pointing (squint) and radiation pattern metrics:

| Progressive phase ($\Phi_{\mathrm{prog}}$) | Peak gain (dBi) | Beam squint ($^\circ$) | Sidelobe level (dB) | Front-to-back ratio (dB) |
| :----------------------------------------- | :-------------- | :--------------------- | :------------------ | :----------------------- |
| $-15.0^\circ$                              | $12.77$         | $0.0$                  | $-12.65$            | $19.45$                  |
| $-10.0^\circ$                              | $12.17$         | $0.0$                  | $-12.05$            | $21.59$                  |
| $-5.0^\circ$                               | $11.39$         | $-5.0$                 | $-11.75$            | $22.89$                  |
| $0.0^\circ$                                | $10.95$         | $-5.0$                 | $-11.36$            | $26.30$                  |
| $5.0^\circ$                                | $10.63$         | $-5.0$                 | $-10.98$            | $27.75$                  |
| $10.0^\circ$                               | $10.76$         | $-5.0$                 | $-10.20$            | $28.46$                  |
| $15.0^\circ$                               | $10.16$         | $-5.0$                 | $-10.14$            | $27.57$                  |

At $\Phi_{\mathrm{prog}} = -15.0^\circ$, the phase compensation correctly aligns the slot phase centres, steering the main beam precisely to broadside ($0.0^\circ$ squint) with a peak gain of $12.77\,\mathrm{dBi}$ and an azimuth sidelobe level (SLL) of $-12.65\,\mathrm{dB}$.

Interestingly, steering the beam to broadside simultaneously minimizes input reflection ($S_{11} = -12.40\,\mathrm{dB}$). This is because the feed phase compensation cancels out the transmission phase shifts introduced by the slot-coupled junctions, aligning both the radiation phases and the reflection phase cancellation at the input. Under other phase conditions (e.g. $\Phi_{\mathrm{prog}} = 0.0^\circ$), the beam squints to $-5.0^\circ$ and the realized gain degrades to $10.95\,\mathrm{dBi}$, primarily due to input mismatch loss rather than pattern degradation. Throughout the sweep, **polarization purity remains exceptionally high**, specifically over 30 dB across the entire azimuth fan ($\varphi\in[-90^\circ,90^\circ]$, $\theta\in[-15^\circ,15^\circ]$), validating the shielding benefit of the aperture-coupled ground plane.

### Architectural implications and centre-fed transition

The severe bandwidth limitation ($0.82\%$) demonstrates that while the local thinning and quarter-wave transformer approach successfully achieves high slot coupling, it destroys the bandwidth benefits of the travelling-wave array. To resolve this, two potential design paths are identified:

1. **Pivot to a centre-fed series-corporate hybrid topology (Recommended):**
	- The input port is split symmetrically at the centre of the array (e.g. using a broadband corporate T-junction on an underlying routing layer) to feed two identical series branches propagating outwards.
	- Because the binomial/Chebyshev amplitude distribution peaks in the centre and decays towards the edges, the peak power requirement coincides with the physical feed location. Feeding from the centre means the highest radiation is demanded where the feedline power is also at its peak. The required coupling coefficients ($C_i = P_{\mathrm{rad},i}/P_{\mathrm{inc},i}$) are thus drastically reduced and flattened.
	- Lower coupling requirements remove the need for high-ratio matching transformers or thinned lines. Keeping the feedline at its nominal width ($Z_0$) throughout preserves the wide impedance bandwidth of the array.
	- Progressive phase delays on the two symmetric branches are equal and opposite, causing frequency-induced beam squint over a chirp sweep to cancel out.
2. **Implement transverse shunt-stub matching networks:**
	- Replace the longitudinal quarter-wave transformers with localized L- or Pi-section matching networks constructed from transverse open shunt stubs.
	- While stubs can provide a slightly more broadband impedance match and save longitudinal feedline space, they add layout complexity and increase the risk of spurious coupling with adjacent patch cavities.

> [!tip]+ Investigate terminal element first
> Before diving into the options of centre-fed travelling-wave arrays, I should investigate the question of how to improve the RAPR capabilities of a single terminal element. This element’s performance will be more easily verified in the current side-fed architecture, where the residual power poses greater problems; therefore, the effect should be more visible. Ultimately, combining an improved ‘virtual load’ element with centre-fed architecture should pave the way to a very good antenna design for my system prototype.

---

[^1]: This applies to standard designs operating with a strongly sub-resonant slot. In our modified topologies – such as dual-polarized elements – the slot must be larger (closer to resonance) to utilize double-resonance tuning. This disrupts the idealized ‘separation of concerns’ between the feed and patch layers.
[^2]: Y. Kang, E. Noh, and K. Kim, ‘Design of Traveling-Wave Series-Fed Microstrip Array With a Low Sidelobe Level’, IEEE Antennas and Wireless Propagation Letters, vol. 19, no. 8, pp. 1395–1399, Aug. 2020, doi: 10.1109/LAWP.2020.2989916.
[^3]: For $Z_0$, we utilize the simulated line impedance ($Z_\mathrm{line}$) rather than an ideal 50 Ω to avoid systematic phase errors.
