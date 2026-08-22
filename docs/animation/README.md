# KLASSIFY explainer animation

A ~2 minute [manim](https://www.manim.community/) animation walking through the
algorithm, one section per pipeline step.

| Section        | Shows                                                          |
| -------------- | -------------------------------------------------------------- |
| `S0Problem`    | an F1 read whose two halves come from two homologous chromosomes |
| `S1Build`      | `build` — keeping only k-mers unique to a single chromosome      |
| `S2Classify`   | `classify` — the 300 / 50% / 10% chimeric-read filter            |
| `S3Control`    | `regions` — subtracting the parent reads to drop artefacts       |
| `S4Breakpoint` | `breakpoint` — the prefix/suffix scan that locates the switch    |
| `S5Pair`       | `cluster-pairs` — pairing the two landing sites                  |
| `S6Recap`      | the pipeline end to end                                          |

## Rendering

```bash
pip install manim          # also needs ffmpeg, and LaTeX for the one formula
./build.sh                 # 1080p60 -> klassify_algorithm.mp4 + .gif
./build.sh -ql             # 480p15, much faster while iterating
```

`build.sh` renders each section, concatenates them into `klassify_algorithm.mp4`,
and derives `klassify_algorithm.gif` from it. Only the GIF is checked in — the
`.mp4` and manim's `media/` directory are ignored, since `build.sh` regenerates
them. The GIF is the one the top-level `README.md` embeds, so it is deliberately
kept small: 900 px wide, 10 fps, 64 colours, no dithering.

A single section:

```bash
manim -qh klassify_algorithm.py S4Breakpoint
```

## Staying honest

Every threshold on screen is hard-coded here to match the current Rust default.
Nothing parses `src/*.rs`, so the two have to be changed together:

| On screen         | Source                                              |
| ----------------- | --------------------------------------------------- |
| `k = 24`          | `KMER_SIZE` in `src/build.rs`                       |
| 300 / 50% / 10%   | `KMER_THRESHOLD`, `SCORE_THRESHOLD`, `MINOR_SCORE_THRESHOLD` in `src/classify.rs` |
| 10 kb bins        | `BINSIZE` in `src/utils.rs`                          |
| 5-100x support    | `MIN_READ_SUPPORT`, `MAX_READ_SUPPORT` in `src/regions.rs` |
| 30 k-mers a side  | `DEFAULT_KMER_THRESHOLD` in `src/breakpoint.rs`      |

`choose_breakpoint()` in `klassify_algorithm.py` is a direct port of the
function of the same name in `src/breakpoint.rs`, down to the `ab_max >= ba_max`
tie-break. The score curve in `S4Breakpoint` is that function's real output on
the synthetic hit track in `demo_hits()` — the peak position, the score, and the
per-side k-mer counts are computed, not drawn by hand. If the Rust changes, the
curve changes with it.
