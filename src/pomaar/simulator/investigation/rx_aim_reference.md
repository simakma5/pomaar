# Receiver aiming changes the span only (reference, 2026-10-21)

**Claim.** In the multi-receiver scenes (parametric beams 60 x 60 deg, horizontal displacement in
the x-y principal plane, bistatic angles 0 ... 16 deg, sphere, trihedral and dihedral at 4 m),
turning the receivers towards the target (rxAim 1) instead of leaving them on boresight
(rxAim 0) multiplies each receiver's whole 2 x 2 scattering matrix by one complex scalar g(beta).
Every polarimetric descriptor is invariant under S -> g S, so aiming changes the span by |g|^2 and
nothing else.

**Argument.**

1. Invariance (exact, by construction). Every descriptor of `signatures` other than the span is a
   channel ratio, a phase difference, a span-normalized fraction (Pauli) or a function of the
   span-normalized coherency matrix (H, A, mean alpha), so g cancels.
2. Scalar fit (`rx_aim_scalar_fit.csv/png` in `data/processed/bistaticity/parametric_beams/rx_aim/`,
   from `aim_check.py`): least squares
   g = <S_bore, S_aim> / ||S_bore||^2 per receiver, angle and frequency, residual
   ||S_aim - g S_bore|| / ||S_aim||:

   | Design | residual, lobe median | lobe max | max anywhere (nulls) |
   | --- | --- | --- | --- |
   | Sphere | -96 dB | -79 dB | -79 dB |
   | Trihedral roll 0 / 90 | -94 / -90 dB | -49 / -48 dB | -37 / -39 dB |
   | Dihedral roll 0 / 22.5 / 45 / 67.5 / 90 | -88 / -103 / -99 / -84 / -79 dB | -62 / -65 / -70 / -55 / -50 dB | -47 / -32 / -39 / -39 / -39 dB |

   ("lobe": span within 20 dB of monostatic, at 79 GHz.) The residual peaks sit at the nulls of
   the targets' bistatic patterns, where the antenna gain slope across the target (about 0.1 dB
   per degree at 16 deg against the 0.6 deg the 40 mm target subtends at 4 m) re-weights
   cancelling contributions. That weighting is still the same for every polarization.
3. g is the receive beam's gain. A receiver on boresight sees the target at beta off its axis
   (the right angle sits at the Tx), and |g|^2 = 3 (beta / 30 deg)^2 dB (the -3 dB point of a
   60 deg beam) to within 0.1 dB for all eight designs (0.83 dB at 16 deg), with |arg g| < 0.8 deg.
   Path lengths are equal (aiming only turns each receiver about its own phase centre).

The descriptor-level view is `summary/summary_aim_difference_<family>.png` (same results folder): differences
of 0.01 ... 0.1 dB and below 0.6 deg inside the lobe. Larger values occur only in nulls and in
Pauli fractions in dB that sit at the floor (e.g. the sphere's k2 at -80 dB).

**Scope: why this is sufficient here and not in general.** It holds because the parametric
beams have identical H and V patterns and no cross-pol, and because the displacement lies in a
principal plane, where the antenna's polarization basis is not rotated. For real antennas the
receive side becomes S_meas = P_Rx(theta) S P_Tx(0). P_Rx(theta) is the 2 x 2 polarization
matrix of the receive pattern at the off-axis angle theta = beta, and aiming changes it:

* the E- and H-plane gains differ (in the azimuth plane the H element sees its E-plane pattern,
  the V element its H-plane pattern), which scales |VV / HH|;
* the pattern's own cross-pol at theta leaks co-pol into cross-pol;
* off the principal planes, the Ludwig-3 basis turns relative to the scattering basis (second
  order, about az x el).

These are one-way antenna-pattern properties. They can be read from the far-field patterns
directly, with no SBR+ run, and they apply equally to a monostatic target off boresight.
