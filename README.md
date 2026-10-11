<img width="2560" height="1600" alt="test_12msolzsol" src="https://github.com/user-attachments/assets/60f62b73-7f21-4b11-a7cc-9a41437eb4be" />
# Neo-Kippenhahn: "star slice" Kippenhahn diagrams for MESA

A classic **Kippenhahn diagram** puts time (model number) on the x-axis and mass on the y-axis. It packs a lot in, but it's hard to read if you're new to stellar evolution. It's a bit like slicing a hot dog the long way.

This program draws the same information as **a cross-section of the star, one image per model**. Flip through the images, or turn them into a movie, and you can watch the star's insides change over time.

In short, instead of looking at stellar evolution exclusively through conventional diagrams, viewers can follow changes in nuclear energy generation, mixing regions, and core boundaries through a sequence of stellar slices.

Each image shows:

| Part of the image | What it shows |
|---|---|
| **Left half** | Energy generation: red where fusion releases energy, blue where neutrinos carry energy away |
| **Right half** | Mixing: convection, overshoot, semiconvection, thermohaline mixing, ... |
| **Dashed circles** | Boundaries of the helium, carbon-oxygen, and iron cores |
| **Header** | Model number, age, mass, central composition, central temperature, luminosity |

---

## 1. What you need

- A finished (or partly finished) **MESA run** with a `LOGS` folder.
- **Python 3.10 or newer**, plus the `numpy` and `matplotlib` packages.

The easiest setup uses [uv](https://docs.astral.sh/uv/). It installs the packages for you the first time you run the script:

```bash
uv run kipp_slices_v3.py
```

If you'd rather use a regular virtual environment:

```bash
python3 -m venv .venv
.venv/bin/pip install numpy matplotlib
.venv/bin/python kipp_slices_v3.py
```

> **Using VS Code's Run button?** Press `Cmd+Shift+P` → **Python: Select Interpreter** → pick `./.venv/bin/python`. Otherwise you'll get `ModuleNotFoundError: No module named 'matplotlib'`.

---

## 2. Set up MESA so it saves the right data

You only have to do this once, *before* running MESA.

### History data (needed for `--source history`)
In `history_columns.list`, make sure these two lines are turned on, with no `!` in front:

```
burning_regions 80
mixing_regions 40
```

To save **every** model instead of every 5th (MESA's default), add this to the `&controls` section of your inlist:

```fortran
history_interval = 1
```

### Profile data (needed for `--source profiles`)
Your `profile_columns.list` needs these columns: `mass`, `radius`, `logT`, `mix_type`, `eps_nuc`, `non_nuc_neu`.

MESA only writes a profile every `profile_interval` steps (100 in this test suite), and **you get one image per profile file**. Lowering `profile_interval` gives smoother animations, but the files add up fast: each one can be several MB.

### ⚠️ Multi-stage runs: stop stages from erasing each other
Many MESA projects, including this test suite, run several inlists in a row (pre-main sequence → ZAMS → helium burning → ...). **Each new stage overwrites `LOGS/history.data`**, so by the end only the last stage is left.

The fix is to give each stage its own folder. Add a line to the `&controls` section of each stage's inlist:

```fortran
&controls
   log_directory = 'LOGS_3_he_burn'
   ...
```

Number the folders in stage order (`LOGS_1_...`, `LOGS_2_...`, ...). The script finds every `LOGS_*` folder automatically and stitches them into one timeline.

---

## 3. Running it

Run these commands from the folder that holds your `LOGS` folders:

```bash
# Every model from history.data: smoothest animation (recommended)
uv run kipp_slices_v3.py --source history

# One image per profile file, mass coordinate (default)
uv run kipp_slices_v3.py

# One image per profile file, real physical radius
uv run kipp_slices_v3.py --coord radius

# Quick preview: every 100th model
uv run kipp_slices_v3.py --source history --every 100

# Only a range of models
uv run kipp_slices_v3.py --source history --start 1000 --end 1200

# Pick specific LOGS folders yourself
uv run kipp_slices_v3.py --logs LOGS_1_pre_zams LOGS_2_to_zams
```

Images are saved as `model_000001.png`, `model_000002.png`, ... in a folder named after the mode, such as `slices_history_mass/`.

### All options

| Option | Default | What it does |
|---|---|---|
| `--source` | `profiles` | Read from `profiles` or `history` |
| `--coord` | `mass` | Disk radius = `mass` coordinate, or physical `radius` (profiles only) |
| `--logs` | all `LOGS_*` folders, or `LOGS` | Which LOGS folders to read, in stage order |
| `--outdir` | `slices_<source>_<coord>/` | Where to save the images |
| `--every` | `1` | Draw every Nth model |
| `--start`, `--end` | whole run | First and last model number to draw |
| `--jobs` | all CPU cores | How many images to draw at once |
| `--dpi` | `120` | Image resolution (higher = bigger files) |

### Turning the images into a movie
With [ffmpeg](https://ffmpeg.org/) installed:

```bash
ffmpeg -framerate 30 -pattern_type glob -i 'slices_history_mass/*.png' -pix_fmt yuv420p kipp.mp4
```

---

## 4. How to read the images

### Mixing (right half)
The colors match MESA's built-in `pgstar` Kippenhahn plot:

| Color | Mixing type | What's going on |
|---|---|---|
| `#87CEFA` light blue | Convection | Hot gas rises, cool gas sinks, like boiling water |
| `#F5F5DC` beige | Overshoot | Convective blobs coasting a bit past the edge of a convective zone |
| `#708090` slate gray | Semiconvection | Slow, partial mixing where a composition gradient resists convection |
| `#CC99FF` lilac | Thermohaline | Slow mixing driven by heavier gas sitting on top of lighter gas |
| dark gray | Radiative | No mixing: energy moves outward as light, not by moving gas |

### Energy generation (left half)
Each colored ring shows the **net** energy rate: energy from fusion minus energy lost to neutrinos, in erg per gram per second.

- **Red = fusion wins.** Darker red means more energy. "10³–10⁴" means between 1,000 and 10,000 erg/g/s.
- **Blue = neutrino cooling wins.** This shows up in late stages, after core helium burning, once the core is hot and dense enough to make lots of neutrinos.
- **Dark gray** means less than 1 erg/g/s either way.

The color scale is fixed for the whole run, so the same shade of red means the same energy rate in every image. (MESA's `pgstar` rescales colors for each plot, so its colors can shift.)

### Mass coordinate vs. physical radius
- **`--coord mass`** (default): distance from the center = how much mass is inside that point, in solar masses. The scale is fixed to the star's starting mass, so you can **watch the star shrink** as winds blow mass away. This view shows the core most clearly.
- **`--coord radius`**: distance from the center = real size, in solar radii. This is the "true" picture, but for a red supergiant the core turns into a tiny dot, because the envelope is about 1000× bigger than the core. Seeing that for yourself is worth it once.

---

## 5. Things to know

- **History vs. profiles:** `history.data` has a row for *every* step (smooth animations) but only stores layer boundaries, so it can't do the physical-radius view. Profiles store every layer in full detail, radius included, but usually only every 100 steps.
- **History merges small layers.** MESA keeps at most 40 mixing layers and 80 burning layers per model in `history.data`, and thin layers near the surface can get merged together. Profiles don't have this limit.
- **Repeated model numbers** (from restarts, or two stages that overlap) are handled for you: the newest one wins.

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'matplotlib'` | You're using a Python without the packages. Run with `uv run ...` or `.venv/bin/python ...` (see section 1). |
| `No mix_type_* columns in ...` | Turn on `mixing_regions` and `burning_regions` in `history_columns.list`, then rerun MESA. |
| `... has no 'mix_type' column` | Add the missing column to `profile_columns.list`, then rerun MESA. |
| The images only cover part of the evolution | Each stage overwrote `history.data`. Give each stage its own `log_directory` (section 2). |
| Only one image every 100 models | That's how often MESA wrote profiles. Use `--source history`, or lower `profile_interval`. |
| No blue (cooling) anywhere | Normal before core helium runs out. Neutrino cooling only takes over in later stages. |

