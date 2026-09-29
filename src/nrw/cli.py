"""Command-line interface.

    python -m nrw.cli --input data/raw/sample01.s2p --thickness-mm 2.0 \\
        --method iterative_eps_mu --output results/sample01

Writes results.csv (all quantities, flags and uncertainty columns), summary.json
(settings, warnings, physics and KK statements) and figures into --output.
"""
from __future__ import annotations

import argparse
import sys

from .api import METHODS, extract_material_parameters
from .branch import METHODS as BRANCH_METHODS
from .io import read_measurement
from .uncertainty import UncertaintySpec


def build_parser():
    p = argparse.ArgumentParser(prog="nrw", description="Extract complex eps_r, mu_r from T/R S-parameters.")
    p.add_argument("--input", required=True, help=".s2p (Touchstone v1) or .csv file")
    p.add_argument("--thickness-mm", type=float, required=True, help="sample length (mm)")
    p.add_argument("--method", default="auto", help=f"one of {METHODS} (hyphens allowed)")
    p.add_argument("--branch-method", default="anchored", help=f"one of {sorted(BRANCH_METHODS)}")
    p.add_argument("--geometry", default="coax", choices=["coax", "waveguide"])
    p.add_argument("--waveguide-a-mm", type=float, default=None, help="broad-wall width a (mm) for TE10")
    p.add_argument("--offset1-mm", type=float, default=0.0, help="port-1 reference plane to sample face (mm)")
    p.add_argument("--offset2-mm", type=float, default=0.0, help="port-2 reference plane to sample face (mm)")
    p.add_argument("--holder-length-mm", type=float, default=None,
                   help="sample-holder length (mm); enables the reference-plane-invariant BJ formulation")
    p.add_argument("--nonmagnetic", action="store_true", help="assume mu_r = 1 (method auto -> baker_jarvis)")
    p.add_argument("--noise-sigma", type=float, default=1e-3, help="rms complex noise per S-parameter")
    p.add_argument("--fmin-ghz", type=float, default=None)
    p.add_argument("--fmax-ghz", type=float, default=None)
    p.add_argument("--uncertainty", action="store_true", help="run Monte Carlo uncertainty propagation")
    p.add_argument("--mc-draws", type=int, default=200)
    p.add_argument("--s-mag-rel", type=float, default=0.005, help="MC: relative |S| std (both S11, S21)")
    p.add_argument("--s-phase-deg", type=float, default=0.5, help="MC: phase std in degrees")
    p.add_argument("--thickness-std-um", type=float, default=5.0)
    p.add_argument("--ref-plane-std-um", type=float, default=5.0)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--touchstone-backend", default="builtin", choices=["builtin", "skrf"])
    p.add_argument("--output", required=True, help="output directory")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    data = read_measurement(args.input, backend=args.touchstone_backend)
    f = data.frequency
    sel = slice(None)
    if args.fmin_ghz is not None or args.fmax_ghz is not None:
        import numpy as np

        lo = -np.inf if args.fmin_ghz is None else args.fmin_ghz * 1e9
        hi = np.inf if args.fmax_ghz is None else args.fmax_ghz * 1e9
        sel = (f >= lo) & (f <= hi)
    geometry = {"kind": args.geometry}
    if args.geometry == "waveguide":
        if args.waveguide_a_mm is None:
            print("error: --waveguide-a-mm is required for --geometry waveguide", file=sys.stderr)
            return 2
        geometry["a_m"] = args.waveguide_a_mm * 1e-3
    unc = None
    if args.uncertainty:
        unc = UncertaintySpec(s11_mag_rel=args.s_mag_rel, s11_phase_deg=args.s_phase_deg,
                              s21_mag_rel=args.s_mag_rel, s21_phase_deg=args.s_phase_deg,
                              thickness_m=args.thickness_std_um * 1e-6, ref_plane_m=args.ref_plane_std_um * 1e-6,
                              n_draws=args.mc_draws, seed=args.seed)
    res = extract_material_parameters(
        frequency=f[sel], s11=data.s11[sel], s21=data.s21[sel], s12=data.s12[sel], s22=data.s22[sel],
        thickness=args.thickness_mm * 1e-3, geometry=geometry, method=args.method.replace("-", "_"),
        branch_method=args.branch_method.replace("-", "_"),
        reference_plane_offsets=(args.offset1_mm * 1e-3, args.offset2_mm * 1e-3),
        holder_length=None if args.holder_length_mm is None else args.holder_length_mm * 1e-3,
        assume_nonmagnetic=args.nonmagnetic, noise_sigma=args.noise_sigma, uncertainty=unc)
    res.metadata["input_file"] = str(args.input)
    res.metadata["touchstone_z0_ohm"] = data.z0
    out = res.save(args.output, plots=not args.no_plots)
    print(f"wrote {out}/results.csv, summary.json" + ("" if args.no_plots else " and figures"))
    for w in res.warnings:
        print(f"WARNING: {w}")
    for q, r in res.causality.items():
        print(f"KK[{q}]: {r.statement}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
