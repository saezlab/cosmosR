# Decompression, Extraction, Translation, and Score Explanation

Use one identifier domain at a time. Most post-processing errors arise from mixing a compressed score table with an uncompressed PKN or a translated display table with internal-ID inputs.

## Decompress scores before using an original PKN

`decompress_moon_result()` does not replace virtual node IDs. It appends `source_original`. Use the final compressed PKN actually scored, then build a fresh score table for original IDs:

```r
mapped_scores <- cosmosR::decompress_moon_result(
  moon_res,
  meta_network_compressed_list = compressed,
  meta_network = pkn_final
)

scores_original <- mapped_scores[, c("source_original", "score", "level")]
names(scores_original)[1] <- "source"
stopifnot(!anyDuplicated(scores_original$source))
```

Keep virtual IDs in a separate provenance column if needed. Do not pass the unmodified `mapped_scores` to an original-ID reducer: its `source` column still contains virtual node names.

Use an uncompressed PKN that is cleaned, coverage-pruned, and coherently filtered in the same original-ID domain as `scores_original`. Do not reduce against the package-wide raw PKN. If compression made it impossible to project coherence-filtered edge identities safely, retain an uncompressed scoring run instead.

When the RNA coherence loop was used and protected IDs make the original mapping valid, apply the same final score-based filter to the uncompressed PKN before extraction:

```r
pkn_original_final <- cosmosR::filter_incohrent_TF_target(
  scores_original, tf_reg_net, pkn_uncompressed, rna_input
)
```

Record this projected filter explicitly. It keeps a raw PKN from reintroducing TF-target edges that the scored workflow rejected; it does not make the RNA-only helper protein-aware.

`decompress_solution_network()` is a different legacy helper for CARNIVAL-style `Activity`/`NodeType` tables. Do not use it on raw MOON output.

## Extract a solution network

Use the single-threshold extractor for its explicit structural pruning:

```r
solution <- cosmosR::reduce_solution_network(
  decoupleRnival_res = scores_original,
  meta_network = pkn_original_final,
  cutoff = cutoff,
  upstream_input = upstream_input,
  RNA_input = rna_input,
  n_steps = n_steps
)
```

It retains nodes with `abs(score) > cutoff`, keeps sign-consistent edges, limits forward distance from upstream inputs, removes non-input pure roots, and removes non-level-0 pure leaves. Its output schemas are:

- `SIF`: `source`, `target`, `interaction`, `consistency`;
- `ATT`: `nodes`, `score`, `level`, `RNA_input`.

Audit the result independently. The single-threshold reducer can retain a sign-consistent cycle disconnected from level-0 observations. Require every retained component to have at least one valid directed route from an intended upstream seed to a level-0 node.

Use `reduce_solution_network_double_thresh()` only after independently auditing connectivity. It keeps edges incident to a primary-threshold node with both endpoints above the secondary threshold. The current code performs its seed-to-level-0 restriction only after it has removed an incoherent edge, so fully coherent disconnected components can survive. Require `primary_thresh >= secondary_thresh` and run a separate path audit.

The double-threshold `ATT` schema is different: `source`, `score`, `level`, `type`, `RNA_input`. Normalize it explicitly before combining it with a single-threshold output.

## Explain one MOON score

Use an identity-matched scoring PKN and score table:

```r
explanation <- cosmosR::get_moon_scoring_network(
  upstream_node = node_id,
  meta_network = pkn_final,
  moon_scores = moon_res,
  keep_upstream_node_peers = FALSE
)
```

This derives the number of downstream layers from the node’s score level, retains controllable and observable neighbors, and optionally omits same-level peers. It is appropriate for explaining which scored downstream nodes were eligible to contribute to one result. It is not a reduced global solution network, a pathway analysis, or a causal validation.

## Translate display copies last

Keep internal `source`/`target` IDs in all calculation and audit tables. Once the network is final:

```r
display_res <- cosmosR::translate_res(
  SIF = solution$SIF,
  ATT = solution$ATT,
  HMDB_mapper_vec = HMDB_mapper_vec
)
```

`translate_res()` edits the first two SIF columns and first ATT column positionally. It changes metabolite and enzyme labels and is not reversible; unmapped IDs can lose prefixes or collide. `translate_column_HMDB()` is useful for a standalone display column, but also requires an explicit named mapping vector.

Write both forms:

1. raw SIF/ATT with stable internal identifiers and the mapping/provenance table;
2. translated SIF/ATT for human-facing figures and Cytoscape-style export.

Never join translated labels back to the PKN without a retained one-to-one mapping check.
