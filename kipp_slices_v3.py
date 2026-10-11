# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""
Render MESA output as a series of "star slice" images, one per model.

A classic Kippenhahn diagram puts time on the x-axis and mass on the y-axis.
This script instead draws each model as a disk, like cutting the star in half:

  - left half  = net energy generation (red = fusion wins, blue = neutrino cooling wins)
  - right half = mixing (convection, overshoot, semiconvection, thermohaline, ...)
  - dashed circles = He / CO / Fe core boundaries

Two data sources:

  --source profiles (default)  reads LOGS/profile*.data. One frame per profile, but
                               every zone is available, including physical radius.
  --source history             reads LOGS/history.data. One frame per model (much
                               smoother animation), but mass coordinate only, since
                               history.data stores boundaries as q = m / M_star.

Two radial coordinates:

  --coord mass (default)  disk radius = mass coordinate m (Msun). The scale is fixed
                          to the largest star mass, so you can watch winds strip it.
  --coord radius          disk radius = physical radius r (Rsun), profiles only.
                          The scale follows each model's surface, since the star
                          grows from a few Rsun to ~1000 Rsun.

Usage:
    uv run kipp_slices_v3.py                                # profiles, mass -> slices_profiles_mass/
    uv run kipp_slices_v3.py --coord radius                 # profiles, physical radius
    uv run kipp_slices_v3.py --source history --every 10    # every 10th history model
    uv run kipp_slices_v3.py --start 1000 --end 1200
    uv run kipp_slices_v3.py --logs LOGS_1_pre_zams LOGS_2_to_zams   # pick stages by hand
    ffmpeg -framerate 30 -pattern_type glob -i 'slices_profiles_mass/*.png' -pix_fmt yuv420p kipp.mp4
"""

import argparse
import glob
import os
import re
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Wedge

# --- Colors -----------------------------------------------------------------

# Mixing colors, copied from star/private/pgstar_colors.f90.
MIX_COLORS = {
    1: ("#87CEFA", "Convection"),
    2: ("#F5F5DC", "Overshoot"),
    3: ("#708090", "Semiconvection"),
    4: ("#CC99FF", "Thermohaline"),
    5: ("#20B2AA", "Rotational"),
    6: ("#CD5C5C", "Rayleigh-Taylor"),
    7: ("#FF7F50", "Minimum"),
    8: ("#D2B48C", "Anonymous"),
    9: ("#0066FF", "Leftover convection"),
}

# pgstar's 51-step burning colormap: blue (index 1) -> near white -> red (index 51).
BURN_CMAP = [
    "#0000FF", "#0909FF", "#1212FF", "#1B1BFF", "#2424FF", "#2F2FFF", "#3838FF",
    "#4141FF", "#4B4BFF", "#5454FF", "#5D5DFF", "#6666FF", "#6F6FFF", "#7878FF",
    "#8181FF", "#8C8CFF", "#9595FF", "#9E9EFF", "#A7A7FF", "#B0B0FF", "#B9B9FF",
    "#C3C3FF", "#CCCCFF", "#DADAFF", "#E9E9FF", "#FFF8F8", "#FFE9E9", "#FFDADA",
    "#FFCCCC", "#FFC3C3", "#FFB9B9", "#FFB0B0", "#FFA7A7", "#FF9E9E", "#FF9595",
    "#FF8A8A", "#FF8181", "#FF7878", "#FF6F6F", "#FF6666", "#FF5D5D", "#FF5454",
    "#FF4B4B", "#FF4141", "#FF3838", "#FF2D2D", "#FF2424", "#FF1B1B", "#FF1212",
    "#FF0909", "#FF0000",
]

SPACE_COLOR = "#000000"  # figure background
STAR_COLOR = "#1C1C1C"   # star body where nothing is burning or mixing
INK = "#E6E6E6"          # primary text
MUTED_INK = "#9A9A9A"    # secondary text

CORE_LINES = [
    ("he_core_mass", "He core", "--"),
    ("co_core_mass", "CO core", "-."),
    ("fe_core_mass", "Fe core", ":"),
]

SCALAR_KEYS = ["star_mass", "star_age", "center_h1", "center_he4", "center_c12",
               "center_o16", "log_center_T", "log_L"] + [k for k, _, _ in CORE_LINES]


def burn_color(burn_type, bmax, bmin):
    """Same mapping as draw_burn_for_step in star/private/pgstar_kipp.f90."""
    n = len(BURN_CMAP)
    mid = n // 2
    if burn_type > 0:
        frac = 1.0 - max(0.0, min(1.0, burn_type / bmax))
        index = n - int(0.6 * frac * (n - mid))
    else:
        frac = 1.0 - max(0.0, min(1.0, burn_type / bmin))
        index = 1 + int(0.6 * frac * mid)
    return BURN_CMAP[index - 1]


def burn_label(burn_type):
    """burn_type = n means 10^(n-1) <= |eps_nuc - eps_neu| < 10^n erg/g/s."""
    lo, hi = abs(burn_type) - 1, abs(burn_type)
    kind = "fusion" if burn_type > 0 else "ν cooling"
    return f"{kind} $10^{{{lo}}}$–$10^{{{hi}}}$"


# --- Reading MESA files -----------------------------------------------------

def read_mesa_table(path):
    """Read a MESA history/profile file into (header dict, column dict)."""
    with open(path) as f:
        f.readline()
        header_names = f.readline().split()
        header_vals = f.readline().split()
        f.readline()
        f.readline()
        names = f.readline().split()
    header = {}
    for k, v in zip(header_names, header_vals):
        try:
            header[k] = float(v)
        except ValueError:
            header[k] = v.strip('"')
    data = np.loadtxt(path, skiprows=6, ndmin=2)
    return header, {name: data[:, i] for i, name in enumerate(names)}


def load_history(path):
    """One model per history row. Regions are (type, m_inner, m_outer) in Msun."""
    _, cols = read_mesa_table(path)

    # After a restart, MESA appends rows whose model numbers overlap earlier ones.
    # Keep only the last occurrence of each model number.
    last = {}
    for i, m in enumerate(cols["model_number"].astype(int)):
        last[m] = i
    rows = [last[m] for m in sorted(last)]

    def slots(prefix):
        n = 0
        while f"{prefix}_type_{n + 1}" in cols:
            n += 1
        if n == 0:
            raise SystemExit(f"No {prefix}_type_* columns in {path}. Add "
                             f"{'mixing' if prefix == 'mix' else 'burning'}_regions "
                             "to history_columns.list.")
        return [(cols[f"{prefix}_type_{i}"], cols[f"{prefix}_qtop_{i}"])
                for i in range(1, n + 1)]

    mix_slots, burn_slots = slots("mix"), slots("burn")

    def regions(slot_cols, i, M, unused):
        out, qbot = [], 0.0
        for types, qtops in slot_cols:
            t, qtop = types[i], min(qtops[i], 1.0)
            if t == unused:
                break
            out.append((int(t), qbot * M, qtop * M))
            qbot = qtop
            if qbot >= 1.0:
                break
        return out

    models = []
    for i in rows:
        M = float(cols["star_mass"][i])
        m = {k: float(cols[k][i]) if k in cols else 0.0 for k in SCALAR_KEYS}
        m.update(model_number=int(cols["model_number"][i]),
                 mix=regions(mix_slots, i, M, unused=-1),
                 burn=regions(burn_slots, i, M, unused=-9999))
        models.append(m)
    return models


def zone_runs(types, inner, outer):
    """Merge neighboring zones with the same type into (type, inner, outer) runs."""
    out = []
    start = 0
    for k in range(1, len(types) + 1):
        if k == len(types) or types[k] != types[start]:
            out.append((int(types[start]), float(inner[start]), float(outer[k - 1])))
            start = k
    return out


def load_profiles(logdir):
    """One model per profile file. Regions are (type, m_inner, m_outer) in Msun."""
    paths = glob.glob(os.path.join(logdir, "profile*.data"))
    paths.sort(key=lambda p: int(re.search(r"profile(\d+)\.data$", p).group(1)))
    if not paths:
        raise SystemExit(f"No profile*.data files in {logdir}")

    by_model = {}  # later profile numbers win if a model number repeats
    for path in paths:
        header, cols = read_mesa_table(path)
        for need in ("mass", "radius", "mix_type", "eps_nuc", "non_nuc_neu"):
            if need not in cols:
                raise SystemExit(f"{os.path.basename(path)} has no '{need}' column. "
                                 "Add it to profile_columns.list.")

        # Zone 1 is the surface. Flip to center -> surface. 'mass' and 'radius'
        # are at each zone's outer face, so the inner face is the next zone down.
        mass = cols["mass"][::-1]
        radius = cols["radius"][::-1]
        m_center = 0.0  # no inner boundary cut out of this model
        m_inner = np.concatenate([[m_center], mass[:-1]])

        # burn_type, computed the same way as set_burn_types in star/private/history.f90.
        net = (cols["eps_nuc"] - cols["non_nuc_neu"])[::-1]
        with np.errstate(divide="ignore"):
            mag = np.where(net != 0, 1 + np.log10(np.abs(net)), 0.0)
        burn_type = (np.sign(net) * np.floor(np.maximum(0.0, mag))).astype(int)
        mix_type = cols["mix_type"][::-1].astype(int)

        m = {k: float(header[k]) if k in header else 0.0 for k in SCALAR_KEYS}
        if "logT" in cols:
            m["log_center_T"] = float(cols["logT"][-1])
        if "photosphere_L" in header:
            m["log_L"] = float(np.log10(header["photosphere_L"]))
        m.update(model_number=int(header["model_number"]),
                 mix=zone_runs(mix_type, m_inner, mass),
                 burn=zone_runs(burn_type, m_inner, mass),
                 mass_grid=np.concatenate([[m_center], mass]),
                 radius_grid=np.concatenate([[0.0], radius]))
        by_model[m["model_number"]] = m
    return [by_model[k] for k in sorted(by_model)]


# --- Drawing ----------------------------------------------------------------

def draw_frame(frame):
    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=frame["dpi"])
    fig.patch.set_facecolor(SPACE_COLOR)
    ax.set_facecolor(SPACE_COLOR)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0.22, right=0.78, top=0.88, bottom=0.08)

    R = frame["rmax"] * 1.05
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)

    S = frame["surface"]
    ax.add_patch(Circle((0, 0), S, color=STAR_COLOR, zorder=0))

    # Left half: burning. Skip burn_type 0 (below 1 erg/g/s), like pgstar does.
    for t, lo, hi in frame["burn"]:
        if t == 0 or hi <= lo:
            continue
        ax.add_patch(Wedge((0, 0), hi, 90, 270, width=hi - lo,
                           color=burn_color(t, frame["bmax"], frame["bmin"]),
                           linewidth=0, zorder=1))

    # Right half: mixing. mix_type 0 is radiative, so leave it as the star body.
    for t, lo, hi in frame["mix"]:
        if t <= 0 or hi <= lo:
            continue
        ax.add_patch(Wedge((0, 0), hi, -90, 90, width=hi - lo,
                           color=MIX_COLORS.get(t, ("#FFFFFF", ""))[0],
                           linewidth=0, zorder=1))

    # Surface outline and the divider between the two halves.
    ax.add_patch(Circle((0, 0), S, fill=False, edgecolor=INK, linewidth=1.2, zorder=3))
    ax.plot([0, 0], [-S, S], color=INK, linewidth=0.8, zorder=3)

    for label, style, r, mass in frame["cores"]:
        ax.add_patch(Circle((0, 0), r, fill=False, edgecolor=INK,
                            linestyle=style, linewidth=1.0, zorder=4))
        ax.text(0, r, f" {label} {mass:.2f} M$_\\odot$", color=INK, fontsize=8,
                va="bottom", ha="left", zorder=5)

    ax.text(-S * 0.5, S * 1.04, "Energy generation", color=INK,
            ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text(S * 0.5, S * 1.04, "Mixing", color=INK,
            ha="center", va="bottom", fontsize=11, weight="bold")

    # Header text.
    age = frame["star_age"]
    age_str = f"{age / 1e6:.4f} Myr" if age >= 1e5 else f"{age:,.0f} yr"
    fig.text(0.5, 0.95, f"Model {frame['model_number']}   ·   t = {age_str}   ·   "
             f"M = {frame['star_mass']:.3f} M$_\\odot$", color=INK, ha="center", fontsize=13)
    fig.text(0.5, 0.915,
             f"center:  X(H) = {frame['center_h1']:.3f}   Y(He) = {frame['center_he4']:.3f}   "
             f"X(C) = {frame['center_c12']:.3f}   X(O) = {frame['center_o16']:.3f}   "
             f"log T$_c$ = {frame['log_center_T']:.2f}   log L = {frame['log_L']:.2f}",
             color=MUTED_INK, ha="center", fontsize=9)
    fig.text(0.5, 0.03, frame["footer"], color=MUTED_INK, ha="center", fontsize=8)

    # Legends: burning on the left, mixing on the right.
    burn_handles = [Patch(color=burn_color(t, frame["bmax"], frame["bmin"]),
                          label=burn_label(t)) for t in frame["burn_levels"]]
    leg = fig.legend(handles=burn_handles, loc="center left", bbox_to_anchor=(0.01, 0.5),
                     title="erg g$^{-1}$ s$^{-1}$", frameon=False, fontsize=8,
                     title_fontsize=9, labelcolor=INK)
    leg.get_title().set_color(INK)

    mix_handles = [Patch(color=c, label=name) for c, name in
                   (MIX_COLORS[t] for t in frame["mix_levels"])]
    mix_handles.append(Patch(facecolor=STAR_COLOR, edgecolor=MUTED_INK, label="Radiative"))
    mix_handles += [Line2D([], [], color=INK, linestyle=s, label=l) for _, l, s in CORE_LINES]
    fig.legend(handles=mix_handles, loc="center right", bbox_to_anchor=(0.99, 0.5),
               frameon=False, fontsize=8, labelcolor=INK)

    out = os.path.join(frame["outdir"], f"model_{frame['model_number']:06d}.png")
    fig.savefig(out, facecolor=SPACE_COLOR)
    plt.close(fig)
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", choices=["profiles", "history"], default="profiles")
    p.add_argument("--coord", choices=["mass", "radius"], default="mass",
                   help="radial coordinate of the disk (radius needs --source profiles)")
    p.add_argument("--logs", nargs="+", default=None,
                   help="MESA LOGS directories, in stage order (default: all LOGS_* folders)")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: slices_<source>_<coord>/)")
    p.add_argument("--every", type=int, default=1, help="render every Nth model")
    p.add_argument("--start", type=int, default=None, help="first model_number")
    p.add_argument("--end", type=int, default=None, help="last model_number")
    p.add_argument("--jobs", type=int, default=os.cpu_count(), help="parallel processes")
    p.add_argument("--dpi", type=int, default=120)
    args = p.parse_args()

    if args.coord == "radius" and args.source == "history":
        p.error("--coord radius needs --source profiles (history.data has no radius info)")

    if args.outdir is None:
        args.outdir = os.path.join(here, f"slices_{args.source}_{args.coord}")

    # Default: every per-stage LOGS_* folder in order (LOGS_1_..., LOGS_2_..., ...),
    # falling back to plain LOGS/ if there aren't any.
    if not args.logs:
        args.logs = sorted(glob.glob(os.path.join(here, "LOGS_*"))) or [os.path.join(here, "LOGS")]

    # Stitch stages together. If two stages share a model number, the later one wins.
    by_model = {}
    for logdir in args.logs:
        if args.source == "history":
            path = os.path.join(logdir, "history.data")
            if not os.path.exists(path):
                print(f"  skipping {logdir} (no history.data)")
                continue
            stage = load_history(path)
        else:
            stage = load_profiles(logdir)
        print(f"  {os.path.basename(logdir)}: {len(stage)} models "
              f"({stage[0]['model_number']}-{stage[-1]['model_number']})")
        by_model.update((m["model_number"], m) for m in stage)
    if not by_model:
        raise SystemExit("No models found.")
    models = [by_model[k] for k in sorted(by_model)]

    # Burning colors are scaled to the whole run (pgstar scales to the plotted
    # window), so a given color means the same thing in every frame.
    all_burn = np.array([t for m in models for t, _, _ in m["burn"]])
    bmax = max(0.1, all_burn.max())
    bmin = min(-0.1, all_burn[all_burn < -1].min()) if (all_burn < -1).any() else -0.1
    burn_levels = sorted({int(t) for t in all_burn if t != 0}, reverse=True)
    mix_levels = sorted({t for m in models for t, _, _ in m["mix"] if t in MIX_COLORS})

    models = [m for m in models
              if (args.start is None or m["model_number"] >= args.start)
              and (args.end is None or m["model_number"] <= args.end)][::args.every]

    max_mass = max(m["star_mass"] for m in models)
    frames = []
    for m in models:
        if args.coord == "radius":
            to_r = lambda x, m=m: float(np.interp(x, m["mass_grid"], m["radius_grid"]))
            surface = m["radius_grid"][-1]
            rmax = surface
            footer = "radius = physical radius r (R$_\\odot$), scaled to fit each frame"
        else:
            to_r = lambda x: x
            surface = m["star_mass"]
            rmax = max_mass
            footer = "radius ∝ mass coordinate m (M$_\\odot$), not physical radius"

        cores = [(label, style, to_r(m[key]), m[key]) for key, label, style in CORE_LINES
                 if m.get(key, 0.0) > 0.01 * m["star_mass"]]
        frame = {k: m[k] for k in SCALAR_KEYS + ["model_number"]}
        frame.update(
            outdir=args.outdir, dpi=args.dpi, surface=float(surface), rmax=float(rmax),
            footer=footer, bmax=bmax, bmin=bmin,
            burn_levels=burn_levels, mix_levels=mix_levels, cores=cores,
            mix=[(t, to_r(lo), to_r(hi)) for t, lo, hi in m["mix"]],
            burn=[(t, to_r(lo), to_r(hi)) for t, lo, hi in m["burn"]])
        frames.append(frame)

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Rendering {len(frames)} models from {args.source} ({args.coord} coordinate) "
          f"to {args.outdir}  (burn_type range {bmin:g} to {bmax:g})")
    with Pool(args.jobs) as pool:
        for n, _ in enumerate(pool.imap_unordered(draw_frame, frames, chunksize=8), 1):
            if n % 100 == 0 or n == len(frames):
                print(f"  {n}/{len(frames)}")


if __name__ == "__main__":
    main()
