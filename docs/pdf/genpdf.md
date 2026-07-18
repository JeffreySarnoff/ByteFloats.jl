# genpdf.md — Building a publication-quality PDF from `docs/` (ByteFloats.jl)

A cleaned, reordered, and hardened version of the process actually used to produce
`ByteFloats.jl-0.1.0-documentation.pdf`. Dead ends from the original run have been
removed, fixes that were discovered iteratively are folded into a single up-front
transformation pass, and every stage now has explicit pass/fail criteria.

The pipeline has five stages:

```
0. Preflight        environment + version checks (fail fast)
1. Resolve          Documenter → fully expanded LaTeX source
2. Transform        one idempotent post-processing pass over the .tex
3. Compile          clean latexmk build + log triage
4. Validate         render → inspect → remediate → repeat to fixpoint
5. Package          PDF + editable source + build note
```

---

## Stage 0 — Preflight

Fail here, not twenty minutes in. Every item below caused (or would have caused)
a mid-pipeline failure in the original run.

```bash
# 0.1  Julia version: read the package's compat BEFORE installing anything.
#      (Original run installed Julia 1.11 first; Project.toml requires "1.12".)
grep 'julia *=' Project.toml
#  → install exactly a matching release, e.g.
curl -sL https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.4-linux-x86_64.tar.gz \
  | tar xz -C /opt
export PATH=/opt/julia-1.12.4/bin:$PATH

# 0.2  TeX toolchain — all four are required (minted needs shell-escape + Pygments).
which xelatex latexmk pygmentize || exit 1

# 0.3  Fonts — Documenter's LaTeX style hardcodes DejaVu via fontspec.
fc-list | grep -qi "DejaVu Sans Mono" || exit 1

# 0.4  Renderer/inspector for the validation stage.
which pdftoppm pdfinfo qpdf || exit 1
pip install pdfplumber pypdf pillow --break-system-packages -q
```

Also read `docs/make.jl` now and note three things that shape later stages:
the `pages = [...]` ordering (the PDF must preserve it), `warnonly` /
`checkdocs` settings (they tell you which warnings are expected), and whether
`remotes = nothing` is set (required for non-git or shallow-clone builds under
Documenter 1.x).

---

## Stage 1 — Resolve the documentation with Documenter

Do **not** hand-convert the Markdown. `@autodocs`, `@docs`, `@ref`, and the
navigation structure must be resolved by Documenter itself, with the package
loaded, or the API Reference and cross-references will be wrong. Documenter's
LaTeX writer with `platform = "none"` emits the fully resolved `.tex` without
attempting to compile it — that file is the single source of truth for
everything downstream.

```bash
# 1.1  Instantiate the docs environment (dev the package from the repo root).
julia --project=docs -e '
  import Pkg
  Pkg.develop(Pkg.PackageSpec(path = pwd()))
  Pkg.instantiate()'
```

```julia
# 1.2  docs/make_latex.jl — identical pages/settings to docs/make.jl,
#      only the format differs.
using Documenter, ByteFloats

makedocs(;
    sitename = "ByteFloats.jl",
    modules  = [ByteFloats],
    authors  = "ByteFloats.jl contributors",
    format   = Documenter.LaTeX(platform = "none", version = v"0.1.0"),
    remotes  = nothing,
    build    = "build_latex",
    pages    = [                      # ← copy verbatim from docs/make.jl
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "User Guide" => "user_guide.md",
        "User Examples" => "user_examples.md",
        "Technical Guide" => "technical_guide.md",
        "Technical Examples" => "technical_examples.md",
        "Adding Operations" => "new_operations.md",
        "API Reference" => "reference.md",
    ],
    checkdocs = :none,
    warnonly  = [:cross_references, :missing_docs],
)
```

```bash
julia --project=docs docs/make_latex.jl
# → docs/build_latex/{ByteFloats.jl-0.1.0.tex, documenter.sty, preamble.tex, custom.sty}
```

**Gate 1.** The build log must end at `LaTeXWriter: creating the LaTeX file`
with only the warnings `warnonly` predicted. The `.tex` must contain zero
literal `@docs` / `@autodocs` / `@ref` strings.

### 1.3  Survey the resolved source (drives every later decision)

One analysis script, run once, replaces the original run's scattered probes.
Everything the Transform and Compile stages need to know is measurable here:

```python
import re, collections
tex = open('ByteFloats.jl-0.1.0.tex').read()

blocks = re.findall(r'\\begin\{minted\}.*?\\end\{minted\}', tex, re.S)
print("code blocks:", len(blocks),
      "longest (lines):", max(b.count('\n') for b in blocks),
      "longest line (chars):", max(len(l) for b in blocks for l in b.split('\n')))
print("display math:", tex.count('\\['), "tables:", tex.count('\\begin{tabulary}'),
      "admonitions:", tex.count('admonition'), "figures:", tex.count('includegraphics'))
print("non-ascii:", ''.join(sorted(set(c for c in tex if ord(c) > 127))))

# duplicate labels (identical section titles on different pages collide)
labels = re.findall(r'\\label\{(\d+)\}', tex)
print("duplicate labels:", [l for l,c in collections.Counter(labels).items() if c > 1])
```

For this repository the survey established: 103 code blocks, longest 36 lines
(→ every block fits one page, so blocks may be made *unbreakable*), longest
code line 101 chars (→ widen margins + enable `breaklines`), no display math,
2 tables, 8 admonitions, no figures, a Unicode inventory needing four fallback
glyphs, and 4 duplicate labels. **If the numbers differ on a future run,
re-derive the decisions** — in particular, a code block or docstring taller
than ~40 lines must remain breakable (rules 7/8) instead of being wrapped.

---

## Stage 2 — Transform: one idempotent pass over the `.tex`

The original run discovered these edits across four compile-inspect cycles.
They are all decidable from the Stage 1.3 survey, so apply them together,
**in this order** (later steps match text produced by earlier ones), in a
single script run against a pristine copy of the resolved `.tex`.

**2.1 Flatten `\part` wrappers.** Documenter emits each nav page as a Part
containing exactly one identically titled chapter; keeping them produces eight
near-blank part-divider pages (a mandatory-rule violation: blank pages from
faulty pagination). Delete every `\part{...}` line; chapters carry the
structure. *Report this reorganization in the build note.*

**2.2 Uniquify duplicate labels.** For each label hash appearing twice,
rename the second occurrence (suffix `-tech` here). Verify first that no
`\hyperlinkref{hash}` references it (true for all four in this document);
if references exist, they must be disambiguated by source position instead.

**2.3 Glue colon-introduced code.** Insert `\nopagebreak[4]` between any
paragraph ending in `:` and an immediately following `\begin{minted}` so the
introducing sentence cannot strand at a page bottom (rule 9):

```python
tex = re.sub(r'(:\s*\n\n+)(\\begin\{minted\})', r'\1\\nopagebreak[4]\n\2', tex)
```

**2.4 Group adjacent code blocks.** Documenter renders a `jlcon` example and
its `text` output as *separate* minted environments; a page break between them
strands the output ("`true`") alone at a page top. Detect chains of minted
blocks separated only by whitespace, **check the chain's total line count
against the page budget (~45 code lines)** — all 27 chains here total ≤ 25
lines — and wrap each qualifying chain in one outer
`\begin{minipage}{\linewidth} ... \end{minipage}`. Chains exceeding the budget
are left unwrapped (splitting *between* their blocks is a legal break).

> ⚠ **fancyvrb constraint** (cost one failed compile originally): nothing may
> follow `\end{minted}` on the same line — not even `%`. Place the closing
> `\end{minipage}` on the next line.

**2.5 Wrap API-reference docstrings.** Each docstring is a *definition block*
(rule 7): `\hypertarget{...}{name} -- {Kind.}` followed by an
`adjustwidth` body. Left alone, headers strand at page bottoms and bodies split.
Estimate each docstring's height (code lines counted directly; prose at ~95
chars/line); wrap those ≤ 38 estimated lines — all 72 here — as one unbreakable
minipage with `\nopagebreak` between header and body. For any docstring over
the threshold, glue only the header to the body's first element and let the
body break between its complete elements. The estimate is conservative and
backstopped: Stage 3 treats any `Overfull \vbox` as a hard failure.

**2.6 Table micro-fixes.** From the survey: change the two tables' column spec
`{R R R}` → `{L L L}` (right-aligned prose columns defeat hyphenation and
inflate minimum widths) and add `\allowbreak` after `\_` inside table cells so
long code identifiers can wrap. Sizing lives in `custom.sty` (2.7).

**2.7 Write `custom.sty`** (loaded by the Documenter preamble after
`documenter.sty`; keep all overrides here so the generated files stay pristine):

```latex
\usepackage{etoolbox}

% Geometry: 1.1in side margins so ~100-char code lines fit at \small mono.
\setulmarginsandblock{1.15in}{1.0in}{*}
\setlrmarginsandblock{1.1in}{1.1in}{*}
\setheaderspaces{0.6in}{*}{*}
\checkandfixthelayout

% Mandatory pagination rules.
\raggedbottom                       % short pages beat bad breaks (rule 9)
\interlinepenalty=10000             % rule 1: paragraphs are unbreakable units
\clubpenalty=10000 \widowpenalty=10000
\displaywidowpenalty=10000 \brokenpenalty=10000
% memoir already emits \nobreak after headings; with unbreakable paragraphs
% this makes heading + first element an indivisible unit (rules 3 and 5).

% Code: unbreakable, wrapped safely, visible continuation marker.
\setminted{breaklines=true, breakindent=1.5em,
  breaksymbolleft={\tiny\ensuremath{\hookrightarrow}},
  fontsize=\small, bgcolor=codeblock-background}
\BeforeBeginEnvironment{minted}{\par\medskip\noindent\begin{minipage}{\linewidth}}
\AfterEndEnvironment{minted}{\end{minipage}\par\medskip}
% (Nested inside the Stage-2.4/2.5 outer minipages: legal and unbreakable.)

% Tables: unbreakable, footnotesize, tightened column padding.
\BeforeBeginEnvironment{tabulary}{\par\medskip\noindent
  \begin{minipage}{\linewidth}\footnotesize\setlength{\tabcolsep}{3pt}}
\AfterEndEnvironment{tabulary}{\end{minipage}\par\medskip}

% Running headers/footers; folio-only on chapter openings.
\pagestyle{ruled}
\makeevenfoot{ruled}{}{\thepage}{}  \makeoddfoot{ruled}{}{\thepage}{}
\makeevenhead{ruled}{\small\leftmark}{}{} \makeoddhead{ruled}{}{}{\small\rightmark}
\copypagestyle{plain}{ruled} \makeevenhead{plain}{}{}{} \makeoddhead{plain}{}{}{}
\copypagestyle{chapter}{ruled} \makeevenhead{chapter}{}{}{}
\makeoddhead{chapter}{}{}{} \makeheadrule{chapter}{\textwidth}{0pt}

% Unicode fallbacks — derive the list from the Stage 1.3 inventory by checking
% each character against DejaVu coverage; these four were absent here.
\usepackage{newunicodechar}
\newunicodechar{⟺}{\ensuremath{\Longleftrightarrow}}
\newunicodechar{≲}{\ensuremath{\lesssim}}
\newunicodechar{⟨}{\ensuremath{\langle}}
\newunicodechar{⟩}{\ensuremath{\rangle}}

% Absorb long unbreakable inline-code tokens instead of overfilling lines.
\emergencystretch=3em

% PDF metadata + numbered, open bookmarks.
\hypersetup{pdftitle={ByteFloats.jl 0.1.0 — Documentation},
  pdfauthor={ByteFloats.jl contributors},
  bookmarksnumbered=true, bookmarksopen=true}
```

Two rules are *not* automated and rely on the compile/validate loop:
residual inline-code margin overflows (fixed with targeted `\allowbreak`
insertions at unambiguous points — after `/` between names, after `(`)
and rule 4 (a section opening too low for heading + intro + first
subsection + its first element), remediated with `\needspace` or
`\clearpage` at the flagged spot. Neither occurred in the final document,
but the checks must run every build.

---

## Stage 3 — Compile and triage the log

**Always build from a clean state.** The original run lost a cycle to a stale
`latexmk` database that reported "up-to-date" while the TOC was unpopulated.

```bash
rm -f main.pdf main.aux main.toc main.out main.fdb_latexmk
latexmk -xelatex -shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

Triage, in severity order — do not proceed to Stage 4 with items 1–3 open:

| Log signal | Meaning | Remedy |
|---|---|---|
| any compile error | often the fancyvrb line rule (2.4) | fix the transform |
| `Overfull \vbox` | some unbreakable unit exceeds the page | unwrap that unit; let it break at element boundaries (rules 7/8) |
| `Missing character` | glyph absent from DejaVu | add a `newunicodechar` fallback |
| multiply-defined labels | duplicate slugs missed in 2.2 | uniquify |
| `Overfull \hbox` > ~10 pt | text may enter the margin | `\allowbreak` at an unambiguous point; for table cells, verify against the rendered page — cell-internal overfulls absorbed by intercolumn space (the residual 4.3 pt here) are harmless |
| empty/short `main.toc` | stale aux state | you skipped the clean step |

**Gate 3.** Zero errors, zero `Overfull \vbox`, zero missing characters, no
duplicate labels, and `main.toc` lists every chapter and section.

---

## Stage 4 — Validate: render, inspect, remediate, repeat

Source-level directives are asserted, not trusted (this is where the original
run caught the split code-output pairs and the stranded docstring headers).
Render **every** page and run the full battery; after any remediation,
recompile and rerun **everything** until a pass with zero findings.

```bash
pdftoppm -png -r 200 main.pdf pages/p     # ≥200 dpi, every page
```

Automated checks (pdfplumber for text/font geometry, PIL for pixels):

1. **Split paragraphs** — belt and suspenders: `Overfull \vbox = 0` proves no
   paragraph exceeded a page (they cannot split under `\interlinepenalty`);
   additionally, every page's last content line must end a complete element
   (sentence punctuation, code line, list item) and every page's first line
   must start one (not a lowercase continuation word). Review each flagged
   boundary by eye — code comments legitimately trip the heuristic.
2. **Stranded headings** — the last content line of a page must not be
   heading-styled. Calibrate detection against a known page first: with
   DejaVu `Scale=MatchLowercase`, pdfplumber reports body ≈ 7.9 pt, sections
   ≈ 9.4 pt bold, chapters ≈ 16 pt (naïve "large font" thresholds miss every
   section). Check **both** patterns: bold sectioning heads *and* the API
   reference's `Name -- {Kind.}` docstring headers.
3. **Split code / tables / admonitions** — no code-background band may touch
   the bottom content edge of page *n* and the top of page *n+1*. When the
   scan fires, confirm against the source whether it is one split block
   (violation) or two adjacent blocks (legal break — but if they are an
   input/output pair, glue them; see 2.4).
4. **Rule 4** — for each section whose first subsection lands on a later page,
   require that the section heading sat high on its page (long intro content);
   flag sections opening in the bottom third.
5. **Blank pages** — no nearly-blank page after the front matter.
6. **Margins** — rightmost ink column per page must respect the 1.1 in margin.
7. **Function** — bookmarks (pypdf outline) mirror the chapter/section tree
   with correct pages; every internal `/GoTo` link target exists among the
   named destinations (68/68 here); TOC page numbers match heading pages;
   `qpdf --check` clean; text extraction succeeds on every page (selectable,
   no U+FFFD or private-use glyphs).

**Gate 4.** All seven checks clean on the same build.

Final result for this document: 53 pages, two remediation cycles
(code-pair gluing; docstring wrapping — both now folded into Stage 2),
then a fully clean pass.

---

## Stage 5 — Package

Deliver three things:

- **`ByteFloats.jl-0.1.0-documentation.pdf`** — the validated PDF.
- **Editable source bundle** — `main.tex` (post-transform), `custom.sty`,
  `documenter.sty`, `preamble.tex`, `make_latex.jl`, and a README giving the
  compile command and summarizing the transforms, so the document can be
  edited at either level: regenerate from `docs/src` via Documenter, or edit
  the resolved LaTeX directly.
- **Build note** — toolchain and versions; reorganizations (part-flattening,
  label uniquification); special handling (docstring/code-pair wrapping,
  table resizing, Unicode fallbacks, dropped web-navigation chrome); and
  anything that could not be represented statically (none here — the source
  docs contain no interactive content).

---

## Appendix — what was cut from the original run, and why

- **Julia 1.11 install** → replaced by the Stage 0 compat check.
- **Four compile-fix-compile cycles** (duplicate labels, table alignment and
  sizing, code-pair gluing, docstring wrapping) → all decidable from the
  Stage 1.3 survey; folded into the single Stage 2 pass. The compile loop now
  exists only for what genuinely needs typeset output: residual hbox overflow
  points and rule-4 placement.
- **Stale-`latexmk` TOC confusion** → replaced by the unconditional clean
  build plus the `main.toc` gate.
- **Incremental `sed` edits to `custom.sty`** (one silently failed to match)
  → the style file is written complete, once; every transform script is
  idempotent and asserts its match counts (`assert n == expected`), so a
  no-op edit fails loudly instead of passing silently.
- **Ad-hoc single-page inspections** → replaced by the fixed Stage 4 battery
  over all pages, so a fix in one place cannot silently regress another.