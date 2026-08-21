# Editorial notes — Bibliography and research integrity

## Main finding

The source ends with `\nocite{*}` and loads both `refs.bib` and `companionII.bib`. Consequently the rendered paper lists 69 references although only 40 keys are explicitly cited in the Paper-I text. Many uncited entries concern Kerr, Thakurta, McVittie or generic mathematical tools and appear to be carry-over from the former combined manuscript or Paper II.

Remove `\nocite{*}` and include only references actually cited and relevant to Paper I. An editor may read the present list as padding or as evidence that the split is incomplete.

## Missing references

Add the direct relativistic-brachistochrone literature listed in `05_literature_codex.md`, especially the 1997, 1998 and 2002 Giannoni/Piccione papers. These are more important for the novelty argument than several uncited Kerr or cosmological-black-hole references currently present.

## Metadata and consistency

- Verify every title, author list, year, volume, pages/article number and DOI against the version of record.
- Prefer persistent DOI/arXiv links to ordinary web URLs.
- Give an exact version/commit/archive for SageMath/abelfunctions and the code release.
- Ensure the Zenodo DOI resolves to the exact snapshot underlying the PDF, not merely a moving project record.
- The bibliography currently reaches reference 69 but occupies only the upper part of the last page; this is another visible sign that it was generated wholesale rather than curated.

## AI disclosure

The acknowledgement identifies the models and uses, states critical review and author responsibility, and is broadly aligned with IOP's current generative-AI disclosure policy. Keep the disclosure. It can be shortened slightly, but do not remove model/version and purpose. Because IOP prohibits relying on generative AI to generate unverified reference lists, retain evidence that all bibliography metadata was independently checked.

## Companion-paper citation

Keep the Paper-II citation only if it is publicly accessible as a preprint at submission. If it is not yet available, describe it as a companion manuscript in text and follow the journal's guidance for unpublished material rather than treating it as an ordinary retrievable reference.
