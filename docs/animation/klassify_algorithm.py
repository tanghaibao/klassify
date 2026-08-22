"""Manim animation explaining how KLASSIFY finds recombination breakpoints.

Render every section (they concatenate into one film, see build.sh)::

    manim -qh klassify_algorithm.py

Render a single section, e.g. the breakpoint refinement::

    manim -qh klassify_algorithm.py S4Breakpoint

Constants shown on screen are the real defaults from the Rust source:
`build.rs` (k = 24), `classify.rs` (300 k-mers / 50% / 10%), `utils.rs`
(10 kb bins, 20 kb chaining), `regions.rs` (support 5-100) and
`breakpoint.rs` (30 k-mers per side).
"""

from __future__ import annotations

from manim import *

config.background_color = "#0d1117"

# --- palette ---------------------------------------------------------------
FG = "#e6edf3"
MUTED = "#8b949e"
DIM = "#30363d"
A_COL = "#58a6ff"  # SoChr01A
B_COL = "#ff7b39"  # SoChr01B
C_COL = "#3fb950"  # SoChr01C
D_COL = "#bc8cff"  # SoChr01D
NOHIT = "#484f58"  # k-mer not unique to any single chromosome
KEEP = "#3fb950"
DROP = "#f85149"
HILITE = "#f0c674"

SAFE_W = 13.4  # keep everything inside the 14.2-unit frame

# --- the algorithm, ported verbatim from src/breakpoint.rs -----------------


def choose_breakpoint(labels, ra, rb, kmer_threshold=30):
    """Port of `choose_breakpoint` in src/breakpoint.rs.

    Returns (idx, left_label, right_label, count_left, count_right, ab, ba)
    where `ab[i] = pref_a[i] + suf_b[i]` and `ba[i] = pref_b[i] + suf_a[i]`.
    """
    n = len(labels)
    eq_a = [lab == ra for lab in labels]
    eq_b = [lab == rb for lab in labels]

    pref_a, pref_b = [0] * n, [0] * n
    pref_a[0], pref_b[0] = int(eq_a[0]), int(eq_b[0])
    for i in range(1, n):
        pref_a[i] = pref_a[i - 1] + int(eq_a[i])
        pref_b[i] = pref_b[i - 1] + int(eq_b[i])

    suf_a, suf_b = [0] * n, [0] * n
    for i in range(n - 1, 0, -1):
        suf_a[i - 1] = suf_a[i] + int(eq_a[i])
        suf_b[i - 1] = suf_b[i] + int(eq_b[i])

    ab = [pref_a[i] + suf_b[i] for i in range(n)]
    ba = [pref_b[i] + suf_a[i] for i in range(n)]
    ab_max, ab_idx = max(ab), ab.index(max(ab))
    ba_max, ba_idx = max(ba), ba.index(max(ba))

    if ab_max >= ba_max:  # same tie-break as the Rust `ab_max >= ba_max`
        idx, left, right = ab_idx, ra, rb
        c_left, c_right = pref_a[ab_idx], suf_b[ab_idx]
    else:
        idx, left, right = ba_idx, rb, ra
        c_left, c_right = pref_b[ba_idx], suf_a[ba_idx]
    return idx, left, right, c_left, c_right, ab, ba


def demo_hits():
    """A synthetic but realistic k-mer hit track for one chimeric read.

    `A` / `B` are the two homologous chromosomes named in the read label,
    `C` is a k-mer unique to some *other* chromosome (counted by neither side).
    """
    labels = ["A"] * 38 + ["B"] * 42
    for i in (11, 25):
        labels[i] = "B"  # noise on the A side
    for i in (5, 33, 70):
        labels[i] = "C"  # unique to a third chromosome
    for i in (44, 61):
        labels[i] = "A"  # noise on the B side
    return labels


LABELS = demo_hits()
IDX, LEFT_LAB, RIGHT_LAB, C_LEFT, C_RIGHT, AB, BA = choose_breakpoint(LABELS, "A", "B")
COLOR_OF = {"A": A_COL, "B": B_COL, "C": C_COL}


# --- shared helpers --------------------------------------------------------


def caption(text: str, size: int = 22, color: str = MUTED):
    return Text(text, font_size=size, color=color).to_edge(DOWN, buff=0.5)


def chrom_bar(label: str, color: str, width: float = 4.6, height: float = 0.32):
    bar = RoundedRectangle(
        corner_radius=0.06,
        width=width,
        height=height,
        stroke_width=2,
        stroke_color=color,
        fill_color=color,
        fill_opacity=0.22,
    )
    lbl = Text(label, font_size=19, color=color).next_to(bar, LEFT, buff=0.22)
    return VGroup(bar, lbl)


def tick_track(labels, cell_w=0.128, cell_h=0.46, gap=0.012, palette=None):
    """A row of small rectangles, one per k-mer position, coloured by origin."""
    palette = palette or COLOR_OF
    cells = VGroup()
    for lab in labels:
        col = palette.get(lab, NOHIT)
        cells.add(
            Rectangle(
                width=cell_w,
                height=cell_h,
                stroke_width=0.6,
                stroke_color=col,
                fill_color=col,
                fill_opacity=0.9 if lab in palette else 0.25,
            )
        )
    cells.arrange(RIGHT, buff=gap)
    return cells


class Section:
    """Mixin giving each numbered section its own banner.

    Deliberately not a Scene subclass, so `manim -a` does not try to render it.
    """

    number = ""
    title = ""
    cmd = None

    def open_section(self):
        num = Text(self.number, font_size=26, color=HILITE, weight=BOLD)
        ttl = Text(self.title, font_size=26, color=FG).next_to(num, RIGHT, buff=0.28)
        row = VGroup(num, ttl)
        if self.cmd:
            code = Text(self.cmd, font_size=18, color=MUTED, font="Menlo")
            code.next_to(ttl, RIGHT, buff=0.5).align_to(ttl, DOWN)
            row.add(code)
        if row.width > SAFE_W:  # long command lines must not run off frame
            row.scale_to_fit_width(SAFE_W)
        row.to_corner(UL, buff=0.45)
        rule = Line(LEFT * 6.9, RIGHT * 6.9, stroke_width=1.2, color=DIM)
        rule.next_to(row, DOWN, buff=0.18).set_x(0)
        header = VGroup(row, rule)
        self.play(FadeIn(header, shift=DOWN * 0.2), run_time=0.6)
        return header


# --- 0 - the problem -------------------------------------------------------


class S0Problem(Scene):
    def construct(self):
        title = Text("KLASSIFY", font_size=66, color=FG, weight=BOLD)
        sub = Text(
            "finding recombination breakpoints with unique k-mers",
            font_size=26,
            color=MUTED,
        ).next_to(title, DOWN, buff=0.35)
        self.play(Write(title, run_time=1.2), FadeIn(sub, shift=UP * 0.2))
        self.wait(0.9)
        self.play(
            title.animate.scale(0.42).to_edge(UP, buff=0.45),
            FadeOut(sub, shift=UP * 0.2),
        )

        # Parental genome: the homologous chromosomes of one polyploid group.
        names = ["SoChr01A", "SoChr01B", "SoChr01C", "SoChr01D"]
        cols = [A_COL, B_COL, C_COL, D_COL]
        bars = VGroup(*[chrom_bar(n, c, width=5.0) for n, c in zip(names, cols)])
        bars.arrange(DOWN, buff=0.40, aligned_edge=LEFT)
        bars.set_x(0.9).set_y(1.15)
        gl = Text("parental genome", font_size=20, color=MUTED)
        gl.next_to(bars, UP, buff=0.3)
        self.play(
            LaggedStart(*[FadeIn(b, shift=RIGHT * 0.2) for b in bars], lag_ratio=0.15),
            run_time=1.1,
        )
        self.play(FadeIn(gl))
        self.wait(0.3)

        # One F1 read that starts on A and ends on B.
        read = VGroup(
            Rectangle(width=2.2, height=0.28, stroke_width=0, fill_color=A_COL, fill_opacity=0.95),
            Rectangle(width=2.2, height=0.28, stroke_width=0, fill_color=B_COL, fill_opacity=0.95),
        ).arrange(RIGHT, buff=0)
        read.set_x(bars[0][0].get_x()).set_y(-1.55)
        rl = Text("one F1 read", font_size=20, color=FG).next_to(read, LEFT, buff=0.32)
        self.play(FadeIn(read, shift=UP * 0.2), FadeIn(rl))
        self.wait(0.3)

        # Each half docks onto the homolog it came from - no crossing arrows.
        seg_a = Rectangle(width=1.6, height=0.28, stroke_width=0,
                          fill_color=A_COL, fill_opacity=0.95)
        seg_a.move_to(bars[0][0].get_center() + LEFT * 1.1)
        seg_b = Rectangle(width=1.6, height=0.28, stroke_width=0,
                          fill_color=B_COL, fill_opacity=0.95)
        seg_b.move_to(bars[1][0].get_center() + RIGHT * 1.1)
        self.play(TransformFromCopy(read[0], seg_a), run_time=1.1)
        self.play(TransformFromCopy(read[1], seg_b), run_time=1.1)
        self.wait(0.3)

        bp = DashedLine(
            read.get_top() + UP * 0.16, read.get_bottom() + DOWN * 0.16,
            stroke_width=2.6, color=HILITE,
        ).set_x(read[0].get_right()[0])
        bp_lbl = Text("breakpoint", font_size=19, color=HILITE).next_to(bp, DOWN, buff=0.12)
        self.play(Create(bp), FadeIn(bp_lbl))

        cap = caption(
            "half the read belongs to one chromosome, half to another - where exactly does it switch?"
        )
        self.play(FadeIn(cap))
        self.wait(2.2)
        self.play(FadeOut(*self.mobjects))


# --- 1 - build -------------------------------------------------------------


class S1Build(Section, Scene):
    number = "1"
    title = "build a table of unique k-mers"
    cmd = "klassify build parents.genome.fa -o kmers.bc"

    def construct(self):
        self.open_section()

        names = ["SoChr01A", "SoChr01B", "SoChr01C"]
        cols = [A_COL, B_COL, C_COL]
        bars = VGroup(*[chrom_bar(n, c, width=5.2) for n, c in zip(names, cols)])
        bars.arrange(DOWN, buff=0.52, aligned_edge=LEFT)
        bars.set_x(0.6).set_y(1.55)
        self.play(LaggedStart(*[FadeIn(b) for b in bars], lag_ratio=0.15), run_time=0.9)

        # Slide a k = 24 window along the first chromosome.
        target = bars[0][0]
        win = Rectangle(width=0.42, height=0.44, stroke_width=2.4,
                        stroke_color=HILITE, fill_opacity=0)
        win.move_to(target.get_left() + RIGHT * 0.26)
        klab = Text("k = 24", font_size=19, color=HILITE).next_to(win, UP, buff=0.14)
        self.play(FadeIn(win), FadeIn(klab))
        self.play(
            win.animate.move_to(target.get_right() + LEFT * 0.26),
            klab.animate.next_to(target.get_right() + LEFT * 0.26, UP, buff=0.14),
            run_time=1.9,
            rate_func=linear,
        )
        self.play(FadeOut(win), FadeOut(klab))

        # Every k-mer lands in one pool; those seen more than once are dropped.
        # (origin chromosome, occurs in more than one chromosome)
        spec = [
            (0, False), (0, False), (0, True), (0, False), (1, True),
            (1, False), (1, False), (2, False), (2, True), (2, False),
            (0, False), (1, False), (2, False), (0, True), (1, False),
            (2, False), (0, False), (1, True), (2, False), (1, False),
        ]
        pool = VGroup(*[
            Square(side_length=0.32, stroke_width=1.4, stroke_color=cols[o],
                   fill_color=cols[o], fill_opacity=0.85)
            for o, _shared in spec
        ])
        pool.arrange_in_grid(rows=2, buff=0.14).set_x(0).set_y(-0.55)
        pool_lbl = Text("all canonical 24-mers", font_size=19, color=MUTED)
        pool_lbl.next_to(pool, UP, buff=0.24)
        self.play(
            FadeIn(pool_lbl),
            LaggedStart(*[FadeIn(s, scale=0.6) for s in pool], lag_ratio=0.05),
            run_time=1.3,
        )
        self.wait(0.3)

        rule = Text(
            "keep a k-mer only if it occurs in exactly one chromosome",
            font_size=22, color=FG,
        ).set_y(-1.75)
        self.play(FadeIn(rule))

        self.play(
            *[
                s.animate.set_fill(NOHIT, opacity=0.35).set_stroke(NOHIT)
                for s, (_o, shared) in zip(pool, spec) if shared
            ],
            run_time=0.8,
        )
        dropped = VGroup(*[s for s, (_o, sh) in zip(pool, spec) if sh])
        kept = VGroup(*[s for s, (_o, sh) in zip(pool, spec) if not sh])
        self.play(LaggedStart(*[Indicate(s, scale_factor=1.25) for s in kept], lag_ratio=0.05),
                  run_time=1.0)
        self.play(FadeOut(dropped, shift=DOWN * 0.6), FadeOut(pool_lbl), run_time=0.7)

        db = RoundedRectangle(corner_radius=0.1, width=3.6, height=0.95,
                              stroke_width=2, stroke_color=HILITE, fill_opacity=0)
        db_lbl = Text("kmers.bc", font_size=22, color=HILITE, font="Menlo")
        db.set_y(-2.35)
        db_lbl.move_to(db)
        self.play(
            kept.animate.arrange(RIGHT, buff=0.08).scale(0.62).set_x(0).set_y(-1.45),
            FadeOut(rule),
        )
        self.play(FadeIn(VGroup(db, db_lbl), shift=UP * 0.1))

        cap = Text("each surviving k-mer is a fingerprint for exactly one chromosome",
                   font_size=22, color=MUTED).set_y(-3.3)
        self.play(FadeIn(cap))
        self.wait(1.8)
        self.play(FadeOut(*self.mobjects))


# --- 2 - classify ----------------------------------------------------------


class S2Classify(Section, Scene):
    number = "2"
    title = "classify reads by k-mer content"
    cmd = "klassify classify kmers.bc f1_reads.fa"

    def construct(self):
        self.open_section()

        legend = VGroup(
            VGroup(Square(0.18, fill_color=A_COL, fill_opacity=0.9, stroke_width=0),
                   Text("unique to SoChr01A", font_size=16, color=MUTED)).arrange(RIGHT, buff=0.14),
            VGroup(Square(0.18, fill_color=B_COL, fill_opacity=0.9, stroke_width=0),
                   Text("unique to SoChr01B", font_size=16, color=MUTED)).arrange(RIGHT, buff=0.14),
            VGroup(Square(0.18, fill_color=NOHIT, fill_opacity=0.3, stroke_width=0),
                   Text("not unique", font_size=16, color=MUTED)).arrange(RIGHT, buff=0.14),
        ).arrange(RIGHT, buff=0.55).set_x(0).set_y(2.35)
        self.play(FadeIn(legend))

        # A normal read: essentially all one chromosome.
        normal = ["A"] * 74
        for i in (7, 30, 55, 66):
            normal[i] = "."
        for i in (21, 48):
            normal[i] = "B"
        # A chimeric read: a run of A followed by a run of B.
        chim = ["A"] * 36 + ["B"] * 38
        for i in (9, 24, 52, 68):
            chim[i] = "."
        chim[14] = "B"
        chim[60] = "A"

        def read_row(labels, name, y):
            track = tick_track(labels, cell_w=0.115, cell_h=0.40, gap=0.01).set_y(y)
            lbl = Text(name, font_size=19, color=FG).next_to(track, LEFT, buff=0.3)
            return VGroup(track, lbl)

        r1 = read_row(normal, "read 1", 1.5)
        r2 = read_row(chim, "read 2", -0.15)
        self.play(FadeIn(r1[1]),
                  LaggedStart(*[FadeIn(c) for c in r1[0]], lag_ratio=0.012), run_time=1.2)
        self.play(FadeIn(r2[1]),
                  LaggedStart(*[FadeIn(c) for c in r2[0]], lag_ratio=0.012), run_time=1.2)

        def tally(labels, anchor):
            a = sum(1 for x in labels if x == "A")
            b = sum(1 for x in labels if x == "B")
            tot = a + b
            txt = Text(f"A {a * 100 // tot}%    B {b * 100 // tot}%",
                       font_size=20, color=FG, font="Menlo")
            txt.next_to(anchor, DOWN, buff=0.22).align_to(anchor, LEFT)
            return txt, b * 100 // tot

        t1, minor1 = tally(normal, r1[0])
        t2, _minor2 = tally(chim, r2[0])
        self.play(FadeIn(t1), FadeIn(t2))
        self.wait(0.4)

        gates = VGroup(
            Text("at least 300 unique k-mers on the read", font_size=20, color=FG),
            Text("major + minor  >=  50% of them", font_size=20, color=FG),
            Text("minor chromosome  >=  10%", font_size=20, color=FG),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT).set_x(0).set_y(-2.35)
        box = SurroundingRectangle(gates, buff=0.28, color=DIM, stroke_width=1.4)
        gate_ttl = Text("chimeric if", font_size=18, color=MUTED)
        gate_ttl.next_to(box, UP, buff=0.12).align_to(box, LEFT)
        self.play(Create(box), FadeIn(gate_ttl),
                  LaggedStart(*[FadeIn(g) for g in gates], lag_ratio=0.15))
        self.wait(0.7)

        # read 1 fails the minor-fraction gate
        cross = Cross(scale_factor=0.26, color=DROP, stroke_width=5)
        cross.next_to(t1, RIGHT, buff=0.75)
        why1 = Text(f"minor is only {minor1}%", font_size=18, color=DROP)
        why1.next_to(cross, RIGHT, buff=0.24)
        self.play(Indicate(gates[2], color=DROP), FadeIn(cross), FadeIn(why1))
        self.play(r1.animate.set_opacity(0.28), t1.animate.set_opacity(0.28))

        # read 2 clears all three
        check = Text("PASS", font_size=20, color=KEEP, weight=BOLD).next_to(t2, RIGHT, buff=0.75)
        tag = Text("SoChr01A@SoChr01B", font_size=20, color=HILITE, font="Menlo")
        tag.next_to(check, RIGHT, buff=0.32)
        self.play(LaggedStart(*[Indicate(g, color=KEEP) for g in gates], lag_ratio=0.2),
                  run_time=1.2)
        self.play(FadeIn(check), FadeIn(tag))
        self.wait(1.8)
        self.play(FadeOut(*self.mobjects))


# --- 3 - parent reads as control -------------------------------------------


class S3Control(Section, Scene):
    number = "3"
    title = "subtract the parent reads"
    cmd = "klassify regions f1_classify.bam parent_classify.bam"

    def construct(self):
        self.open_section()

        # One reference: F1 pileups above it, parent pileups below it.
        ref_w = 11.0
        ref = Rectangle(width=ref_w, height=0.34, stroke_width=1.6, stroke_color=MUTED,
                        fill_color=MUTED, fill_opacity=0.12).set_y(0.2)
        n_bins = 18
        ticks = VGroup(*[
            Line(UP * 0.17, DOWN * 0.17, stroke_width=1, color=DIM)
            .move_to([ref.get_left()[0] + ref_w * i / n_bins, ref.get_y(), 0])
            for i in range(1, n_bins)
        ])
        ref_lbl = Text("parental reference, 10 kb bins", font_size=17, color=MUTED)
        ref_lbl.next_to(ref, DOWN, buff=0.12).align_to(ref, RIGHT)
        self.play(FadeIn(ref), Create(ticks, lag_ratio=0.03), FadeIn(ref_lbl))

        def bin_x(i):
            return ref.get_left()[0] + ref_w * (i + 0.5) / n_bins

        def pileup(x, n, color, y0, up=True):
            d = 1 if up else -1
            return VGroup(*[
                Rectangle(width=0.5, height=0.1, stroke_width=0,
                          fill_color=color, fill_opacity=0.85)
                .move_to([x, y0 + d * j * 0.155, 0])
                for j in range(n)
            ])

        x_real, x_art = bin_x(4), bin_x(12)
        f1_real = pileup(x_real, 8, A_COL, 0.5, up=True)
        f1_art = pileup(x_art, 7, A_COL, 0.5, up=True)
        p_art = pileup(x_art, 6, C_COL, -0.1, up=False)

        f1_lbl = Text("F1 chimeric reads", font_size=18, color=A_COL).set_y(2.15).set_x(-4.6)
        p_lbl = Text("parent chimeric reads (control)", font_size=18, color=C_COL)
        p_lbl.set_y(-1.5).set_x(-4.0)

        self.play(
            LaggedStart(*[FadeIn(s, shift=DOWN * 0.15) for s in f1_real], lag_ratio=0.06),
            LaggedStart(*[FadeIn(s, shift=DOWN * 0.15) for s in f1_art], lag_ratio=0.06),
            FadeIn(f1_lbl),
            run_time=1.3,
        )
        self.wait(0.3)
        self.play(
            LaggedStart(*[FadeIn(s, shift=UP * 0.15) for s in p_art], lag_ratio=0.06),
            FadeIn(p_lbl),
            run_time=1.0,
        )
        self.wait(0.4)

        formula = MathTex(
            r"\text{ratio}=\frac{\text{F1 depth}}{\text{parent}_1+\text{parent}_2+1}",
            font_size=32, color=FG,
        ).set_y(-2.35).set_x(-4.1)
        self.play(FadeIn(formula))

        keep_box = SurroundingRectangle(f1_real, buff=0.12, color=KEEP, stroke_width=2.2)
        keep_txt = VGroup(
            Text("ratio 8", font_size=19, color=KEEP, font="Menlo"),
            Text("real crossover", font_size=18, color=KEEP),
        ).arrange(DOWN, buff=0.1).next_to(keep_box, LEFT, buff=0.3)
        drop_box = SurroundingRectangle(VGroup(f1_art, p_art), buff=0.12,
                                        color=DROP, stroke_width=2.2)
        drop_txt = VGroup(
            Text("ratio 1", font_size=19, color=DROP, font="Menlo"),
            Text("assembly artefact", font_size=18, color=DROP),
        ).arrange(DOWN, buff=0.1).next_to(drop_box, RIGHT, buff=0.3)

        self.play(Create(keep_box), FadeIn(keep_txt))
        self.play(Create(drop_box), FadeIn(drop_txt))
        self.wait(0.9)
        self.play(FadeOut(drop_box), FadeOut(f1_art), FadeOut(p_art), FadeOut(drop_txt),
                  FadeOut(p_lbl), run_time=0.8)

        cap = Text("keep the bins with 5-100x F1 support that the parents do not share",
                   font_size=21, color=FG).set_y(-3.35).set_x(0)
        self.play(FadeIn(cap))
        self.wait(1.8)
        self.play(FadeOut(*self.mobjects))


# --- 4 - breakpoint refinement ---------------------------------------------


class S4Breakpoint(Section, Scene):
    number = "4"
    title = "find the switch point inside the read"
    cmd = "klassify breakpoint kmers.bc regions.fasta"

    def construct(self):
        self.open_section()
        n = len(LABELS)

        track = tick_track(LABELS, cell_w=0.125, cell_h=0.42, gap=0.01).set_y(2.15)
        tlbl = Text("k-mer hits along one chimeric read", font_size=19, color=MUTED)
        tlbl.next_to(track, UP, buff=0.2)
        self.play(FadeIn(tlbl),
                  LaggedStart(*[FadeIn(c) for c in track], lag_ratio=0.008), run_time=1.4)
        self.wait(0.3)

        note = Text("the run of A is not clean, and neither is the run of B",
                    font_size=20, color=FG).next_to(track, DOWN, buff=0.3)
        noise = VGroup(*[
            Circle(radius=0.11, stroke_width=2, color=HILITE).move_to(track[i])
            for i in (11, 25, 44, 61)
        ])
        self.play(FadeIn(note), LaggedStart(*[Create(c) for c in noise], lag_ratio=0.12))
        self.wait(1.1)
        self.play(FadeOut(note), FadeOut(noise))

        question = Text("score every possible cut:   A-hits to the left  +  B-hits to the right",
                        font_size=21, color=FG).next_to(track, DOWN, buff=0.28)
        self.play(FadeIn(question))
        self.wait(0.5)

        axes = Axes(
            x_range=[0, n - 1, 20],
            y_range=[0, max(AB) + 8, 20],
            x_length=10.2,
            y_length=3.0,
            tips=False,
            axis_config={"stroke_width": 1.6, "color": DIM, "include_numbers": False},
        )
        axes.set_x(track.get_x()).set_y(-1.15)
        y_lbl = Text("score", font_size=17, color=MUTED).next_to(axes.y_axis, LEFT, buff=0.12)
        x_lbl = Text("cut position i", font_size=17, color=MUTED).next_to(axes.x_axis, DOWN, buff=0.14)
        self.play(Create(axes), FadeIn(y_lbl), FadeIn(x_lbl))

        curve = VMobject(stroke_width=3.2, stroke_color=HILITE)
        curve.set_points_as_corners([axes.c2p(i, AB[i]) for i in range(n)])

        # Sweep the cut across the read while the score curve draws itself.
        cut = Line(track.get_top() + UP * 0.08, track.get_bottom() + DOWN * 0.08,
                   stroke_width=2.6, color=HILITE)
        cut.move_to([track[0].get_x(), track.get_y(), 0])
        self.play(FadeIn(cut))
        self.play(
            cut.animate.move_to([track[n - 1].get_x(), track.get_y(), 0]),
            Create(curve),
            run_time=3.2,
            rate_func=linear,
        )

        peak = Dot(axes.c2p(IDX, AB[IDX]), radius=0.075, color=KEEP)
        peak_lbl = Text(f"best cut: i = {IDX}    score {AB[IDX]}",
                        font_size=19, color=KEEP, font="Menlo")
        peak_lbl.next_to(peak, UP, buff=0.2)
        self.play(
            cut.animate.move_to([track[IDX].get_x() + 0.07, track.get_y(), 0]),
            FadeIn(peak, scale=0.5),
            FadeOut(question),
            run_time=0.8,
        )
        self.play(FadeIn(peak_lbl))
        self.wait(0.5)

        # The other orientation (B on the left, A on the right) scores far worse.
        alt = VMobject(stroke_width=2.2, stroke_color=MUTED)
        alt.set_points_as_corners([axes.c2p(i, BA[i]) for i in range(n)])
        alt_lbl = Text(f"cutting the other way round only reaches {max(BA)}",
                       font_size=18, color=MUTED)
        alt_lbl.set_x(2.4).set_y(1.05)
        self.play(Create(alt), FadeIn(alt_lbl), run_time=1.2)
        self.wait(0.9)
        self.play(FadeOut(alt), FadeOut(alt_lbl))

        # Support on each side has to clear 30 k-mers.
        left_brace = Brace(VGroup(*track[: IDX + 1]), UP, buff=0.08, color=A_COL)
        right_brace = Brace(VGroup(*track[IDX + 1:]), UP, buff=0.08, color=B_COL)
        left_t = Text(f"{C_LEFT} A-hits", font_size=18, color=A_COL).next_to(left_brace, UP, buff=0.06)
        right_t = Text(f"{C_RIGHT} B-hits", font_size=18, color=B_COL).next_to(right_brace, UP, buff=0.06)
        self.play(FadeOut(tlbl))
        self.play(GrowFromCenter(left_brace), GrowFromCenter(right_brace),
                  FadeIn(left_t), FadeIn(right_t))
        gate = Text("both sides clear the 30 k-mer minimum", font_size=19, color=KEEP).set_y(-3.3)
        self.play(FadeIn(gate))
        self.wait(1.3)

        # Cut the read in two.
        self.play(
            FadeOut(axes), FadeOut(curve), FadeOut(peak), FadeOut(peak_lbl),
            FadeOut(x_lbl), FadeOut(y_lbl), FadeOut(gate),
            FadeOut(left_brace), FadeOut(right_brace), FadeOut(left_t), FadeOut(right_t),
        )
        left_part = VGroup(*track[: IDX + 1])
        right_part = VGroup(*track[IDX + 1:])
        self.play(
            left_part.animate.shift(LEFT * 0.55 + DOWN * 1.9),
            right_part.animate.shift(RIGHT * 0.55 + DOWN * 1.9),
            cut.animate.shift(DOWN * 1.9).set_opacity(0.35),
            run_time=1.0,
        )
        id_left = Text("read|SoChr01A|0-9184", font_size=19, color=A_COL, font="Menlo")
        id_right = Text("read|SoChr01B|9184-18320", font_size=19, color=B_COL, font="Menlo")
        id_left.next_to(left_part, DOWN, buff=0.35)
        id_right.next_to(right_part, DOWN, buff=0.35)
        self.play(FadeIn(id_left), FadeIn(id_right))
        cap = caption("the read is cut midway between the last A k-mer and the first B k-mer")
        self.play(FadeIn(cap))
        self.wait(2.0)
        self.play(FadeOut(*self.mobjects))


# --- 5 - pair the regions --------------------------------------------------


class S5Pair(Section, Scene):
    number = "5"
    title = "map the halves back and pair them up"
    cmd = "klassify cluster-pairs f1_classify.roi.bam"

    def construct(self):
        self.open_section()

        chr_l = chrom_bar("SoChr01B", A_COL, width=9.2, height=0.34)
        chr_r = chrom_bar("SoChr01F", B_COL, width=9.2, height=0.34)
        chrs = VGroup(chr_l, chr_r).arrange(DOWN, buff=2.15, aligned_edge=LEFT)
        chrs.set_x(0.7).set_y(0.9)
        self.play(FadeIn(chr_l), FadeIn(chr_r))

        xa = chr_l[0].get_left()[0] + 3.0
        xb = chr_r[0].get_left()[0] + 5.8

        left_halves, right_halves = VGroup(), VGroup()
        for j in range(6):
            jitter = (j - 2.5) * 0.09
            left_halves.add(
                Rectangle(width=1.15, height=0.11, stroke_width=0,
                          fill_color=A_COL, fill_opacity=0.9)
                .move_to([xa + jitter, chr_l[0].get_y() - 0.42 - j * 0.16, 0])
            )
            right_halves.add(
                Rectangle(width=1.15, height=0.11, stroke_width=0,
                          fill_color=B_COL, fill_opacity=0.9)
                .move_to([xb + jitter, chr_r[0].get_y() + 0.42 + j * 0.16, 0])
            )

        self.play(
            LaggedStart(*[FadeIn(m, shift=UP * 0.2) for m in left_halves], lag_ratio=0.08),
            LaggedStart(*[FadeIn(m, shift=DOWN * 0.2) for m in right_halves], lag_ratio=0.08),
            run_time=1.4,
        )
        note = Text("the two halves of each read land in two different places",
                    font_size=21, color=FG).set_y(-2.5)
        self.play(FadeIn(note))
        self.wait(1.2)

        site_l = SurroundingRectangle(left_halves, buff=0.1, color=HILITE, stroke_width=2)
        site_r = SurroundingRectangle(right_halves, buff=0.1, color=HILITE, stroke_width=2)
        sl = Text("site, 6 reads", font_size=17, color=HILITE).next_to(site_l, LEFT, buff=0.25)
        sr = Text("site, 6 reads", font_size=17, color=HILITE).next_to(site_r, RIGHT, buff=0.25)
        self.play(Create(site_l), Create(site_r), FadeIn(sl), FadeIn(sr))

        link = DashedLine(site_l.get_bottom(), site_r.get_top(), stroke_width=2.4, color=HILITE)
        weight = Text("paired by their shared reads", font_size=18, color=HILITE)
        weight.next_to(link.get_center(), RIGHT, buff=0.4)
        self.play(Create(link), FadeIn(weight))
        self.wait(1.0)

        out = VGroup(
            Text("SoChr01B:71411-81028", font_size=23, color=A_COL, font="Menlo"),
            Text("SoChr01F:81751-88094", font_size=23, color=B_COL, font="Menlo"),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        out_box = SurroundingRectangle(out, buff=0.28, color=KEEP, stroke_width=2)
        out_lbl = Text("f1_classify.roi.paired.regions", font_size=17, color=MUTED, font="Menlo")
        out_lbl.next_to(out_box, UP, buff=0.14)
        out_grp = VGroup(out_box, out, out_lbl).set_x(0).set_y(-2.6)

        self.play(FadeOut(note), FadeOut(sl), FadeOut(sr), FadeOut(weight))
        self.play(FadeIn(out_grp, shift=UP * 0.25))
        self.wait(2.2)
        self.play(FadeOut(*self.mobjects))


# --- 6 - recap -------------------------------------------------------------


class S6Recap(Scene):
    def construct(self):
        title = Text("the whole pipeline", font_size=34, color=FG).to_edge(UP, buff=0.8)
        self.play(FadeIn(title))

        steps = [
            ("1", "build", "unique k-mers, one chromosome each"),
            ("2", "classify", "reads that carry two chromosomes"),
            ("3", "regions", "drop what the parents also show"),
            ("4", "breakpoint", "cut each read at the switch"),
            ("5", "cluster-pairs", "pair the two landing sites"),
        ]
        rows = VGroup()
        for num, cmd, desc in steps:
            n = Text(num, font_size=22, color=HILITE, weight=BOLD)
            c = Text(cmd, font_size=22, color=FG, font="Menlo")
            d = Text(desc, font_size=21, color=MUTED)
            row = VGroup(n, c, d).arrange(RIGHT, buff=0.45)
            rows.add(row)
        rows.arrange(DOWN, buff=0.46, aligned_edge=LEFT).set_x(-0.5).set_y(0.15)

        # keep the command and description columns aligned
        cmd_x = max(r[1].get_left()[0] for r in rows)
        desc_x = max(r[2].get_left()[0] for r in rows)
        for r in rows:
            r[1].shift(RIGHT * (cmd_x - r[1].get_left()[0]))
            r[2].shift(RIGHT * (desc_x - r[2].get_left()[0]))
            r[1].set_y(r[0].get_y())
            r[2].set_y(r[0].get_y())

        self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.25) for r in rows], lag_ratio=0.22),
                  run_time=2.4)
        self.wait(1.0)

        cite = VGroup(
            Text("Zhu et al. (2026)   Nature   doi:10.1038/s41586-026-10863-3",
                 font_size=19, color=MUTED),
            Text("github.com/tanghaibao/klassify", font_size=19, color=MUTED, font="Menlo"),
        ).arrange(DOWN, buff=0.18).to_edge(DOWN, buff=0.75)
        self.play(FadeIn(cite))
        self.wait(2.4)
        self.play(FadeOut(*self.mobjects))
