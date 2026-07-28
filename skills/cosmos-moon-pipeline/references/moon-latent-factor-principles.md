# MOON Latent-Factor Principles

This note is separate from [the data-to-PKN mapping guide](moon-data-pkn-mapping-principles.md)
on purpose. Latent-factor analysis helps decide which signatures are worth
interpreting, but it is not the core rule for mapping data onto the PKN. Its
maintained source is
[`cosmosR/agent-docs/moon-latent-factor-principles.md`](https://github.com/saezlab/cosmosR/blob/master/agent-docs/moon-latent-factor-principles.md);
prefer that source when its current version differs.

Use this file when a workflow starts from MOFA or another factor-analysis model.
It is based on section 2.3 of the COSMOS+ preprint.

## Factor Selection

- Treat the latent factor as the object to mechanistically explain, not
  individual samples or conditions.
- Check whether the factor model is stable before interpreting a factor. In the
  NCI60 example, MOFA was run while varying the maximum number of factors, and
  interpretation focused on a stable solution.
- Prioritize factors that explain variance in the omic views that the network
  analysis will connect.
- Use metadata enrichment on factor coordinates to orient interpretation. The
  section 2.3 example modeled MOFA Z coordinates as a function of binarized
  metadata classes such as tissue of origin.
- Interpret feature-weight signs relative to factor direction. If a factor is
  negatively associated with a metadata class, positive feature weights do not
  automatically mean higher abundance in that class.

## Prior-Knowledge Readiness

- Before building a MOON network, check whether the factor weights are
  consistent with regulatory prior knowledge.
- For RNA factors, score TF-target footprints with `decoupleR` ULM or a similar
  method.
- For ligand/receptor interpretation, score LR sets from the relevant RNA and/or
  protein features.
- In the NCI60 section 2.3 example, factor 4 was selected partly because it
  captured variance across RNA, proteomic, and metabolomic views and was
  especially consistent with TF-target and LR prior knowledge.

## Interpretation Discipline

- A factor is a statistical axis of variation; the biological label is inferred,
  not built in.
- A high TF, LR, or MOON score means the factor aligns with a prior-knowledge
  pattern. It does not prove that the corresponding mechanism is active without
  further evidence.
- Once the factor-specific signatures are selected, use
  [the data-to-PKN mapping guide](moon-data-pkn-mapping-principles.md) to
  decide how each data type should enter MOON.
