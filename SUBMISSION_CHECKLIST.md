# SIADS submission checklist

## 1. Archive the code (do this first, the paper must cite the DOI)

1. Create a public GitHub repository, e.g. `comparator-cusp`, and push the
   contents of this package. Do **not** include the manuscript source or PDF:
   SIAM holds copyright on the accepted version, and the repository is code and
   data only.
2. In `CITATION.cff` and `.zenodo.json`, replace `XXXX-XXXX-XXXX-XXXX` by your
   ORCID and `USER` by your GitHub user name.
3. On zenodo.org, log in with GitHub, enable the repository under
   Settings -> GitHub, then create a GitHub release tagged `v1.0.0`. Zenodo
   mints a version DOI and a concept DOI automatically.
4. Cite the **concept** DOI in the paper, not the version DOI: it always
   resolves to the latest release, which matters if a referee asks for changes.

## 2. Replace the review-stage note

In the manuscript, the section *Code and data availability* begins with
`\emph{Review-stage note.}`. Replace that first sentence by the repository URL
and the concept DOI, and delete the sentence that announces the replacement.
Keep the rest of the paragraph, which describes what is verified and at what
precision.

## 3. Switch the class to submission mode

The preamble is

    \documentclass[review,onefignum,onetabnum]{siamonline250211}

`review` gives double spacing and line numbers, which is what SIADS wants for a
new submission. Keep it. Check that the compiled PDF has line numbers and that
`\headers` produces a short running title that fits.

## 4. Files to upload

- `comparator_cusp_siads.tex`
- `comparator_cusp_siads.bib`
- `figures/singular_set.pdf`, `fold_vs_normal_form.pdf`,
  `normal_form_error.pdf`, `experimental_slice.pdf`
- the compiled PDF
- the SIAM class and style files if the system asks for them

Confirm that the four figure PDFs are the current ones: `experimental_slice.pdf`
must have the shaded bistable wedge with the label inside it, and
`fold_vs_normal_form.pdf` must have no theorem number in its legend.

## 5. Metadata for the submission form

- **Title.** Global fold and cusp geometry of an asymmetric two-gene
  cross-repressive module
- **MSC codes.** 92C42, 37G10, 34C23, 58K05, 37N25
- **Keywords.** genetic toggle switch, mutual repression, fold, cusp, Whitney
  singularities, imperfect pitchfork, center manifold, topological degree,
  competitive systems
- **Section.** Applications of dynamical systems to the life sciences, or
  Bifurcation theory; either is defensible, the first is closer to the framing.

## 6. Suggested referees

Give names working on gene-circuit bifurcation analysis or on singularities of
planar maps, avoiding anyone with a recent joint publication. Natural fits from
the reference list are the authors of the two papers the work is measured
against, provided you are comfortable with that.

## 7. Related-work declaration

SIAM asks whether related manuscripts are under consideration. Two items to
declare if they are live at submission time:

- the corrigendum to Gamermann et al. (2012) intended for *Journal of
  Biological Systems*, which corrects statements no longer discussed here;
- the companion paper on stochastic and information-theoretic aspects of the
  same circuit, announced at the end of the Introduction.

Neither overlaps with the mathematical content of this manuscript, but
declaring them costs nothing and pre-empts a question.

## 8. Cover letter

One page. State the result in two sentences, then the delta with respect to
Richard et al. (2023) and Clement et al. (2026): the factorized singular set for
arbitrary leaks, the complete Whitney classification with a unique cusp, and
the closed-form expansion of the fold about it. Mention that the code is
archived and that every number in the paper is machine-verified.
