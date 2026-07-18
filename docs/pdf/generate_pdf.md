# generate_pdf.md — Publication-quality PDF from any Julia package's `docs/`

A package-agnostic pipeline for turning a Documenter.jl documentation tree into
one coherent, strictly paginated PDF: title page, table of contents, numbered
sections, bookmarks, working cross-references, syntax-highlighted Julia code,
and no split paragraphs, stranded headings, or broken elements.

The pipeline is measurement-driven: nothing is hardcoded to a particular
package. A survey of the resolved LaTeX source decides which transforms apply;
every transform is conditional on what the survey found; and a render-and-
inspect loop asserts the result page by page.

```
0. Preflight        environment + package compatibility checks (fail fast)
1. Resolve          Documenter → fully expanded LaTeX source
2. Survey           measure the source; fill in the decision table
3. Transform        one idempotent, conditional pass over the .tex
4. Compile          clean latexmk build + log triage
5. Validate         render every page → inspect → remediate → fixpoint
6. Package          PDF + editable source + build note
```

Throughout, "the pagination rules" means: paragraphs never split across pages;
page breaks fall only between complete elements; headings and captions stay
with their content; sections don't open without room for their first
subsection's first element; no widows, orphans, stranded headings, or blank
pages from faulty pagination; short pages are always preferred to bad breaks;
code, tables, figures, equations, admonitions, and definition blocks split
only when they cannot fit on one page — and then only deliberately (readable
size reduction, landscape, or labeled logical parts with repeated headers).

---

## Stage 0 — Preflight

Check everything the pipeline will need before spending time on any of it.

```bash
# 0.1  Julia version: the package's compat decides which Julia to install.
grep -A20 '^\[compat\]' Project.toml | grep '^julia'
# Install a matching release from julialang-s3 and put it on PATH.

# 0.2  Docs environment: read docs/Project.toml for Documenter's major version
#      and for plugins that change the build (DocumenterCitations, Literate,
#      DemoCards, DocumenterMermaid, ...). Plugins must be replicated in the
#      LaTeX build or their content will be missing or broken.
cat docs/Project.toml

# 0.3  TeX toolchain: XeLaTeX or LuaLaTeX (fontspec is mandatory for
#      Documenter's LaTeX style), latexmk, and Pygments for minted.
which xelatex latexmk pygmentize

# 0.4  Fonts: documenter.sty hardcodes DejaVu Sans / DejaVu Sans Mono.
fc-list | grep -qi "DejaVu Sans Mono"

# 0.5  Validation tooling.
which pdftoppm pdfinfo qpdf
pip install pdfplumber pypdf pillow fonttools -q
```

Then read `docs/make.jl` and record, without changing anything yet:

- the `pages = [...]` tree — the PDF must preserve this ordering and nesting;
- `modules`, `sitename`, `authors` — reused verbatim;
- `checkdocs` / `warnonly` — the set of warnings that are *expected*;
- whether the build runs doctests or `@example` blocks — these execute the
  package, so the package must load and any binary/system dependencies must
  be present (test with `julia --project=docs -e 'using <Package>'` after
  Stage 1.1 and resolve failures before proceeding);
- any HTML-only machinery (`assets`, `analytics`, custom writers, `size_threshold`
  tuning) — irrelevant to LaTeX, drop it from the LaTeX build script;
- `remotes` / repo detection — for local clones without a full git remote,
  the LaTeX build needs `remotes = nothing` (Documenter ≥ 1) or it errors.

---

## Stage 1 — Resolve the documentation with Documenter

Never hand-convert the Markdown. `@docs`, `@autodocs`, `@ref`, `@example`,
`@index`, and the navigation tree must be resolved by Documenter with the
package loaded, or generated API material and cross-references will be wrong
or missing. The LaTeX writer with `platform = "none"` emits the resolved
`.tex` without trying to compile it; that file is the source of truth for
every later stage.

```bash
# 1.1  Instantiate the docs environment against the local package.
julia --project=docs -e '
  import Pkg
  Pkg.develop(Pkg.PackageSpec(path = pwd()))
  Pkg.instantiate()'
```

```julia
# 1.2  docs/make_latex.jl — a copy of docs/make.jl with only these changes:
#      format → Documenter.LaTeX(platform = "none", version = v"<pkg version>"),
#      build  → "build_latex",
#      remotes = nothing        (if building from a local/shallow clone),
#      drop deploydocs and HTML-only options; KEEP pages, modules, plugins,
#      checkdocs, and warnonly exactly as the package defines them.
```

```bash
julia --project=docs docs/make_latex.jl
# → docs/build_latex/{<SiteName>.tex, documenter.sty, preamble.tex, custom.sty}
```

**Gate 1.** The log ends at `LaTeXWriter: creating the LaTeX file` with only
the warnings that `warnonly` predicted; doctests (if enabled) passed; and the
emitted `.tex` contains zero literal `@docs`/`@autodocs`/`@ref`/`@example`
strings. If Documenter errors on a specific page construct the LaTeX writer
cannot handle (rare: raw HTML blocks, some plugin output), fix it at the
source level (`@raw html` → static alternative) rather than in the `.tex`,
and record it for the build note.

---

## Stage 2 — Survey the resolved source

One script, run once. Its output fills the decision table that drives every
transform in Stage 3. Measure at least:

```python
import re, collections
tex = open('SITE.tex').read()

# Code blocks: count, per-block line counts, longest raw line.
blocks = re.findall(r'\\begin\{minted\}.*?\\end\{minted\}', tex, re.S)
lens   = sorted((b.count('\n') for b in blocks), reverse=True)

# Structure: parts/chapters, tables, figures, math, admonitions, verbatim.
counts = {k: tex.count(k) for k in (
    '\\part{', '\\chapter{', '\\begin{tabulary}', '\\begin{longtable}',
    '\\includegraphics', '\\[', '\\begin{equation', '\\begin{align',
    'admonition', '\\begin{verbatim}')}

# Labels: duplicates collide (identical section titles on different pages).
labels = re.findall(r'\\label\{([^}]+)\}', tex)
dups   = [l for l,c in collections.Counter(labels).items() if c > 1]

# Docstring entries (Documenter's API-reference shape).
docstrings = re.findall(
    r'\\hypertarget\{[^}]*\}\{[^\n]*?\}  -- \{[A-Za-z]+\.\}', tex)

# Unicode inventory, to be checked against actual font coverage (Stage 3.7).
nonascii = sorted(set(c for c in tex if ord(c) > 127))
```

Derive the **page budget** from the geometry you will set in `custom.sty`
(Letter/A4, ~1.1 in side margins, `\small` mono ⇒ roughly 45 code lines or
50 body lines per page; recompute if you change geometry or fonts). Then fill
the decision table:

| Measurement | Decision |
|---|---|
| every `\part` contains exactly one chapter with the same title | flatten parts to chapters (kills near-blank divider pages); otherwise **keep parts** — the package's grouping is real structure |
| duplicate labels exist | uniquify; if any duplicate is the target of `\hyperlinkref`, disambiguate references by source position too |
| longest code block ≤ page budget | all code blocks become unbreakable |
| some code blocks > page budget | those blocks stay breakable (splitting an oversized element is permitted); optionally split them at the *source* level into logically complete parts, each labeled |
| adjacent code blocks separated only by whitespace ("input/output pairs") | glue chains whose combined length fits the budget; leave longer chains breakable between blocks |
| tables present | wrap unbreakable if the whole table fits a page at a readable size (≥ footnotesize); a table taller than a page gets the deliberate treatment: landscape page, or logical parts with the header row repeated and continuations labeled — never scale below comfortable readability |
| figures present | keep each `\includegraphics` and its caption in one unbreakable unit; oversized figures scale to `\linewidth`/page height, not split |
| display math present | keep standard display-math penalties tight (`\predisplaypenalty`, `\postdisplaypenalty` high) so equations stay with their introducing/explaining text |
| docstrings present | wrap each header+body as one unbreakable unit if its estimated height fits the budget; otherwise glue only the header to the body's first element |
| non-ASCII characters present | fallback-map exactly those absent from the document fonts (Stage 3.7) |
| longest raw code line at chosen geometry | if wider than the text block: widen margins first, then rely on `breaklines` with a visible continuation marker; landscape only for pathological cases |

---

## Stage 3 — Transform: one idempotent, conditional pass

Apply the decision table in a single script against a pristine copy of the
resolved `.tex`. Order matters (later steps match text produced by earlier
ones): **3.1 flatten parts → 3.2 uniquify labels → 3.3 glue colon-introduced
code → 3.4 group code chains → 3.5 wrap docstrings → 3.6 table/figure fixes →
3.7 write `custom.sty`.** Every regex substitution asserts its match count
against the survey (`assert n == expected`) so a silently failing edit stops
the build instead of shipping.

**3.1 Flatten `\part` wrappers** — only under the survey condition above.
Report it as a reorganization in the build note.

**3.2 Uniquify duplicate labels** — rename second and later occurrences with
a deterministic suffix; patch references only if the survey found any.

**3.3 Glue colon-introduced code** — a paragraph ending `:` followed by a
code block must not separate from it:

```python
tex = re.sub(r'(:\s*\n\n+)(\\begin\{minted\})', r'\1\\nopagebreak[4]\n\2', tex)
```

**3.4 Group adjacent code blocks** — for each whitespace-separated chain
whose total lines fit the page budget, wrap the chain in one outer
`\begin{minipage}{\linewidth} ... \end{minipage}`.

> ⚠ **fancyvrb constraint:** nothing may follow `\end{minted}` on its line —
> not even `%`. Put the closing `\end{minipage}` on the next line.

**3.5 Wrap docstrings** — estimate height (code lines counted directly; prose
at ≈ chars-per-line for your geometry); wrap fitting docstrings whole, with
`\nopagebreak` between the `\hypertarget` header and the body; for oversized
ones glue only the header. The estimate is deliberately conservative and
backstopped by Stage 4's hard `Overfull \vbox = 0` gate.

**3.6 Table and figure fixes** — per the decision table: column-spec changes
(right-aligned `R` prose columns defeat hyphenation — prefer `L`),
`\allowbreak` after `\_` and other unambiguous points inside wide cells,
figure+caption keep-together wrappers, and the multi-page-table treatment
where the survey demanded it.

**3.7 Write `custom.sty`.** All overrides live here so the Documenter-
generated files stay pristine. The pagination core is package-independent:

```latex
\usepackage{etoolbox}

% Geometry — pick to fit the measured longest code line, then recompute
% the page budget used by Stages 2/3.
\setulmarginsandblock{1.15in}{1.0in}{*}
\setlrmarginsandblock{1.1in}{1.1in}{*}
\setheaderspaces{0.6in}{*}{*}
\checkandfixthelayout

% Mandatory pagination rules.
\raggedbottom                       % short pages beat bad breaks
\interlinepenalty=10000             % paragraphs are unbreakable units
\clubpenalty=10000 \widowpenalty=10000
\displaywidowpenalty=10000 \brokenpenalty=10000
\predisplaypenalty=10000            % equations stay with their intro text
% memoir emits \nobreak after headings; with unbreakable paragraphs this
% makes heading + first element an indivisible unit automatically.

% Code blocks: unbreakable, safe wrapping, visible continuation marker.
\setminted{breaklines=true, breakindent=1.5em,
  breaksymbolleft={\tiny\ensuremath{\hookrightarrow}},
  fontsize=\small, bgcolor=codeblock-background}
\BeforeBeginEnvironment{minted}{\par\medskip\noindent\begin{minipage}{\linewidth}}
\AfterEndEnvironment{minted}{\end{minipage}\par\medskip}
% If the survey found code blocks taller than a page, do NOT install this
% hook globally; instead wrap only the fitting blocks in Stage 3.4-style
% explicit minipages and let the oversized ones break.

% Tables (adjust size per survey; never below comfortable readability).
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

% Unicode fallbacks — GENERATED, not hardcoded: check every survey character
% against the document fonts and map only the misses to math-mode or
% alternate-font equivalents.
\usepackage{newunicodechar}
%% \input{unicode-fallbacks.tex}   (written by the coverage script below)

\emergencystretch=3em               % absorb long inline-code tokens

\hypersetup{pdftitle={<SiteName> <version> — Documentation},
  pdfauthor={<authors>}, bookmarksnumbered=true, bookmarksopen=true}
```

Generate the fallback list programmatically (this is what makes the step
portable across packages with different symbol inventories):

```python
from fontTools.ttLib import TTFont
cmaps = [set(TTFont(p).getBestCmap()) for p in (SANS_PATH, MONO_PATH)]
missing = [c for c in nonascii if not any(ord(c) in cm for cm in cmaps)]
# map each miss to \ensuremath{...} or a \newfontfamily fallback; if no
# reasonable mapping exists, that's a finding for the build note.
```

Two rules stay in the compile/validate loop because they need typeset output:
residual inline-code margin overflows (targeted `\allowbreak` at unambiguous
points) and rule-4 section placement (`\needspace` or `\clearpage` at the
flagged spot).

---

## Stage 4 — Compile and triage the log

Always build from a clean state — a stale latexmk database can report
"up-to-date" over an unpopulated TOC.

```bash
rm -f main.pdf main.aux main.toc main.out main.fdb_latexmk
latexmk -xelatex -shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

| Log signal | Meaning | Remedy |
|---|---|---|
| compile error | often the fancyvrb line rule (3.4) or a plugin's LaTeX | fix the transform / handle the construct at source level |
| `Overfull \vbox` | an unbreakable unit exceeds the page | **hard failure**: unwrap that unit and let it break between its complete elements, or apply the oversized-element treatment |
| `Missing character` | glyph absent despite Stage 3.7 | extend the fallback map |
| multiply-defined labels | duplicates missed in 3.2 | uniquify |
| `Overfull \hbox` beyond ~10 pt | text may enter the margin | `\allowbreak` at an unambiguous point; verify table-cell overfulls against the rendered page — ones absorbed by intercolumn space are harmless |
| empty/short `.toc` | stale aux state | you skipped the clean step |

**Gate 4.** Zero errors, zero `Overfull \vbox`, zero missing characters, no
duplicate labels, and the `.toc` lists every chapter and section.

---

## Stage 5 — Validate: render, inspect, remediate, fixpoint

Source-level directives are asserted, never trusted. Render **every** page at
≥ 200 dpi and run the full battery; after any remediation, recompile and rerun
**everything** until a build passes with zero findings.

```bash
pdftoppm -png -r 200 main.pdf pages/p
```

1. **Split paragraphs** — `Overfull \vbox = 0` proves none (paragraphs cannot
   split under `\interlinepenalty`); additionally each page's last content
   line must end a complete element and each page's first line must start one
   (not a lowercase continuation word). Eyeball flagged boundaries — code
   comments legitimately trip the heuristic.
2. **Stranded headings** — no page's last content line may be heading-styled.
   **Calibrate detection first** on a page with known headings: with fontspec
   `Scale=MatchLowercase`, extracted sizes differ from nominal (e.g. DejaVu
   body ≈ 7.9 pt, sections ≈ 9.4 pt bold), so naïve "large font" thresholds
   miss everything. Check both sectioning heads *and* docstring headers
   (`Name -- {Kind.}`).
3. **Split blocks** — no code/admonition background band touching the bottom
   content edge of page *n* and the top of *n+1*; on a hit, check the source:
   one split block is a violation; two adjacent blocks is a legal break —
   unless they form an input/output pair, which goes back to 3.4. Same scan
   for tables and figures (caption on one page, body on the next = violation).
4. **Rule 4** — a section whose first subsection lands on a later page must
   have opened high on its own page; flag sections opening in the bottom third.
5. **Blank pages** — none after the front matter.
6. **Margins** — rightmost/leftmost ink respects the margins on every page;
   for landscape pages, check rotated.
7. **Function** — bookmarks mirror the full pages tree with correct targets;
   every internal `/GoTo` link resolves to an existing named destination; TOC
   page numbers match heading pages; `qpdf --check` is clean; text extracts
   on every page with no U+FFFD or private-use glyphs; figures render
   (non-blank image regions where `\includegraphics` was emitted).

**Gate 5.** All checks clean on the same build.

---

## Stage 6 — Package

- **The PDF**, named `<SiteName>-<version>-documentation.pdf`.
- **Editable source bundle**: the transformed `main.tex`, `custom.sty` (and
  generated `unicode-fallbacks.tex`), the untouched `documenter.sty` /
  `preamble.tex`, `make_latex.jl`, and a README with the compile command —
  so the document can be edited at either level: regenerate from `docs/src`
  via Documenter, or edit the resolved LaTeX directly.
- **Build note**: toolchain and versions; every reorganization relative to
  `docs/make.jl` (part-flattening, label uniquification, any page reordering);
  special handling (wrapping decisions, table/figure treatment, Unicode
  fallbacks, dropped web-only chrome); content that could not be represented
  statically (interactive widgets, animations, `@raw html`) and what replaced
  it; and any expected warnings inherited from the package's `warnonly`.

---

## Appendix A — package-specific knobs to revisit on each new package

- **Julia version and binary deps** — from `[compat]` and whatever `using
  <Package>` needs (BinaryBuilder artifacts usually just work; system libs
  may not).
- **Documenter plugins** — DocumenterCitations needs its `.bib` and style in
  the LaTeX build; Literate/DemoCards must run *before* `makedocs`; Mermaid
  and other HTML-only diagram plugins need static image replacement.
- **Geometry and page budget** — recompute whenever paper size, margins, or
  fonts change; all fit-vs-wrap decisions depend on it.
- **The unbreakable-paragraph rule** — safe for typical documentation prose;
  a package whose docs contain page-length paragraphs will hit `Overfull
  \vbox`, which is the signal to exempt that paragraph (locally reset
  `\interlinepenalty`) rather than weaken the global rule.
- **Global minted hook** — only when the survey shows every code block fits a
  page; otherwise wrap selectively.
- **`@example`/doctest output** — executes on the build machine; outputs that
  embed timings, paths, or RNG draws will differ between runs. Acceptable for
  a PDF snapshot; note it if reproducibility matters.

## Appendix B — failure modes already encountered once (don't rediscover them)

- Installing Julia before reading `[compat]`.
- Trusting a warm latexmk state (empty TOC that "compiled fine").
- `%` or anything else on the `\end{minted}` line (fancyvrb error).
- Text edits via unasserted regex/sed — a pattern that no longer matches
  fails silently; always assert match counts against the survey.
- "Large font" heading detection without calibrating against fontspec's
  scaled sizes.
- Treating every code-background band at a page boundary as a violation —
  adjacent independent blocks may legally break there; check the source.
- Table-cell `Overfull \hbox` warnings that never reach the page margin —
  verify against the render before "fixing" readability away.