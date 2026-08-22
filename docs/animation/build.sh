#!/usr/bin/env bash
# Render the KLASSIFY explainer, stitch the sections into one film, and derive
# the GIF that README.md embeds.
#
#   ./build.sh          # 1080p60
#   ./build.sh -ql      # 480p15, for quick iteration
#
# Requires manim (pip install manim) and ffmpeg.
set -euo pipefail

QUALITY="${1:--qh}"
case "$QUALITY" in
  -ql) DIR=480p15 ;;
  -qm) DIR=720p30 ;;
  -qh) DIR=1080p60 ;;
  -qk) DIR=2160p60 ;;
  *) echo "unknown quality flag: $QUALITY" >&2; exit 1 ;;
esac

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

SECTIONS=(S0Problem S1Build S2Classify S3Control S4Breakpoint S5Pair S6Recap)

manim "$QUALITY" -a --media_dir media klassify_algorithm.py

VIDEODIR="media/videos/klassify_algorithm/$DIR"
LIST="$(mktemp)"
trap 'rm -f "$LIST"' EXIT
for s in "${SECTIONS[@]}"; do
  printf "file '%s/%s.mp4'\n" "$PWD/$VIDEODIR" "$s" >> "$LIST"
done

ffmpeg -v error -y -f concat -safe 0 -i "$LIST" -c copy klassify_algorithm.mp4
echo "wrote $PWD/klassify_algorithm.mp4"

# The GIF is what README.md embeds, so keep it small: 900px wide (about the
# width GitHub renders a README image at), 10 fps, and a 64-colour palette with
# no dithering, which this flat-coloured design survives without visible loss.
TMPDIR_PAL="$(mktemp -d)"
PALETTE="$TMPDIR_PAL/palette.png"
trap 'rm -f "$LIST"; rm -rf "$TMPDIR_PAL"' EXIT
ffmpeg -v error -y -i klassify_algorithm.mp4 \
  -vf "fps=10,scale=900:-1:flags=lanczos,palettegen=stats_mode=diff:max_colors=64" \
  "$PALETTE"
ffmpeg -v error -y -i klassify_algorithm.mp4 -i "$PALETTE" \
  -lavfi "fps=10,scale=900:-1:flags=lanczos[x];[x][1:v]paletteuse=dither=none:diff_mode=rectangle" \
  klassify_algorithm.gif
echo "wrote $PWD/klassify_algorithm.gif ($(du -h klassify_algorithm.gif | cut -f1))"
