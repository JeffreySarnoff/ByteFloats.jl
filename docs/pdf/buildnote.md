# Build note — ByteFloats.jl-0.1.0.pdf (regenerated 2026-07-24)

Produced by the pipeline in `generate_pdf.md` from the current `docs/src`
(11 pages: Home … Internal Reference, including the new Julia Compatibility
page and the External/Internal reference split). 65 pages, Letter.

## Toolchain

Fedora TeX Live 2025 (`texlive-scheme-medium` + latexmk/xetex/minted), XeLaTeX,
Pygments, DejaVu Sans / DejaVu Sans Mono, Documenter 1.x LaTeX writer
(`platform = "none"` via `docs/make_latex.jl`), validation with
pdftoppm/pdfplumber/pypdf/qpdf at 200 dpi.

## Reorganizations relative to docs/make.jl

- The 11 one-chapter `\part` wrappers were flattened to chapters (each part
  contained exactly one chapter of the same name; kills near-blank dividers).
- 8 duplicate section labels uniquified (identical heading titles on different
  pages); none were link targets, so no references needed patching.

## Special handling

- All 161 code blocks fit the 45-line page budget → globally unbreakable;
  27 whitespace-adjacent block chains (input/output pairs) glued whole;
  73 colon-introduced blocks glued to their introducing paragraph.
- All 116 docstrings wrapped as single unbreakable units, separated by
  breakable glue (without it, consecutive docstring boxes formed unbreakable
  stacks — measured as a 4016 pt Overfull \vbox before the fix).
- 11 tabulary tables: prose columns switched R → L (right alignment defeats
  hyphenation); each table kept whole at footnotesize.
- Long mono tokens: `\allowbreak` inserted after 284 escaped underscores and
  34 slash-joined `\texttt` pairs; one paragraph with two 23-char camelCase
  tokens absorbs its residue via local `\emergencystretch`.
- Unicode: the 62 non-ASCII characters in the source are all covered by
  DejaVu Sans (prose) and DejaVu Sans Mono (code) — no fallback mappings.

## Verification

Gate 4: zero errors, zero `Overfull \vbox`, zero missing characters, no
duplicate labels, full TOC; worst residual `Overfull \hbox` 5.98 pt (absorbed
by intercolumn/margin space, verified on the rendered pages). Gate 5 at
200 dpi over all 65 pages: no blank pages, no ink outside margins, no
code/table background band crossing a page boundary, no heading- or
docstring-header-styled last lines, bookmarks mirror the full page tree,
`qpdf --check` clean, text extracts on every page with no U+FFFD.

## Expected warnings

`checkdocs = :none`, `warnonly = [:cross_references, :missing_docs]` inherited
from `docs/make.jl`; the resolution pass emitted no warnings.

## Reproducing

```
julia --project=docs docs/make_latex.jl     # resolve → docs/build_latex/
# then apply the transform + custom.sty and compile with
# latexmk -xelatex -shell-escape (see generate_pdf.md for the full pipeline)
```

The transformed `main.tex`, `custom.sty`, `transform.py`, and `validate.py`
for this build are archived in the session's build bundle
(`ByteFloats-pdf-source-bundle.tar.gz`); regenerate from `docs/src` for any
content change rather than editing the LaTeX.
