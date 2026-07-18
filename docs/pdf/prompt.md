Prepare a publication-quality PDF from the documentation in this `../docs` directory of the Julia repository.

Use the repository documentation as the authoritative source. 
Inspect the complete documentation structure, including:

* `docs/make.jl`;
* `docs/Project.toml`;
* the Documenter.jl navigation structure;
* all referenced Markdown files;
* relevant source docstrings included through `@docs`, `@autodocs`, or similar directives;
* local figures, tables, mathematical expressions, examples, and code listings.

Build or otherwise resolve the Documenter.jl documentation before preparing the PDF so that generated API material and cross-references are included correctly.

## Required document organization

Produce one coherent PDF with:

* a title page;
* a table of contents;
* numbered sections and subsections;
* PDF bookmarks;
* working internal cross-references;
* syntax-highlighted Julia code;
* clearly formatted notes, warnings, examples, equations, tables, and figure captions;
* page numbers and consistent running headers or footers.

Preserve the logical ordering established by `docs/make.jl` unless a small reorganization is clearly necessary to make the PDF readable. Report any such reorganization.

## Mandatory pagination rules

Intelligently paginate the document. Short or partially filled pages are acceptable and preferable to awkward or misleading page breaks.

Apply these rules strictly:

1. **Never split a paragraph across pages.**
   When an entire paragraph does not fit in the remaining space, move the complete paragraph to the next page.

2. **A subsection may span multiple pages, but page breaks within it may occur only between complete structural elements**, such as:

   * between paragraphs;
   * before or after a code listing;
   * before or after a table together with its caption;
   * before or after a figure together with its caption;
   * before an equation together with its introductory text;
   * after an equation together with its immediately associated explanation;
   * between complete list items, when the list remains understandable.

3. **Never separate a heading from the content it introduces.**

   * Keep each subsection heading with at least its first complete paragraph or other first substantive element.
   * Keep each code-example heading with the associated code listing.
   * Keep each table or figure caption with its table or figure.

4. **Do not begin a section unless the same page has enough remaining space for:**

   * the section heading;
   * any introductory text that precedes the first subsection;
   * the first subsection heading; and
   * at least the first complete paragraph or substantive element of that subsection.

   Otherwise, start the section on the next page.

5. **Do not leave an isolated heading at the bottom of a page.**

6. **Avoid widows and orphans.**
   Because paragraphs must remain intact, resolve these primarily by moving whole paragraphs or other complete elements to the next page.

7. **Do not split code blocks, equations, tables, figures, admonitions, signatures, or definition blocks unless an element is too large to fit on one page.**

8. When a single table or code listing is longer than one page:

   * first attempt a reasonable landscape page, smaller but readable type, or a deliberate division into logically complete parts;
   * repeat table headers on continuation pages;
   * label continuations clearly;
   * never reduce the text below a comfortably readable size merely to avoid a page break.

9. Prefer a short page over:

   * splitting a paragraph;
   * separating a heading from its content;
   * separating introductory text from the equation, table, figure, or code that it introduces;
   * starting a section without enough room for its first subsection.

## Typography and readability

Use a restrained technical-document style appropriate for Julia package documentation.

* Use a comfortably readable body font size and line spacing.
* Use monospaced type that clearly distinguishes `0/O`, `1/l/I`, punctuation, and Julia operators.
* Prevent code lines from running into the margins.
* Wrap code only when doing so remains unambiguous; otherwise use a wider or landscape layout.
* Preserve Unicode mathematical and Julia symbols correctly.
* Keep heading levels visually distinct but not oversized.
* Use adequate whitespace around headings, equations, tables, code listings, and admonitions.
* Ensure links remain recognizable in both color and grayscale printing.
* Do not use decorative layouts that reduce technical readability.

## Content handling

* Do not silently omit any substantive documentation.
* Resolve Documenter.jl directives rather than printing unresolved directives in the PDF.
* Preserve the meaning and formatting of mathematical expressions.
* Include relevant API documentation generated from source docstrings.
* Remove web-navigation elements that have no useful printed equivalent.
* Convert interactive-only content into an appropriate static representation where possible.
* Identify any content that cannot be represented faithfully in a static PDF.
* Correct obvious layout defects, but do not rewrite the technical content unless necessary to repair a clear formatting or rendering problem.

## Validation

Do not assume that source-level pagination directives worked correctly.

After generating the PDF:

1. Render every PDF page to an image at no less than 200 dpi.
2. Inspect every rendered page for:

   * split paragraphs;
   * stranded headings;
   * sections beginning without their first subsection;
   * broken code indentation;
   * clipped or overflowing text;
   * unreadably small tables or code;
   * separated captions;
   * malformed equations;
   * missing glyphs;
   * blank or nearly blank pages caused by faulty pagination;
   * incorrect cross-references or table-of-contents entries.
3. Revise the layout and regenerate the PDF wherever a mandatory pagination rule is violated.
4. Repeat rendering and inspection until the final PDF satisfies all requirements.
5. Verify that the PDF opens correctly and that its table of contents, bookmarks, page numbers, links, and selectable text work.

Provide:

* the final PDF for download;
* the editable intermediate source used to generate it;
* a brief build note identifying the toolchain and any documentation elements that required special handling.
