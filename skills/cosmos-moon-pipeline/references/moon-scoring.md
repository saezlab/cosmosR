# MOON Scoring and Consistency Filtering

Use MOON to score upstream nodes from downstream activities over the PKN. It is a layer-wise signed regulon analysis, not a path optimizer and not experimental proof of causality.

## Respect the current `moon()` contract

```r
moon_res <- cosmosR::moon(
  upstream_input = upstream_input,
  downstream_input = downstream_input,
  meta_network = pkn_scoring,
  n_layers = n_layers,
  statistic = "ulm"
)
```

Require `source`, `target`, and one sign column (`interaction` or `sign`). The current code removes outgoing edges whose source is a downstream input, scores successive upstream layers with decoupleR, and returns `source`, `score`, and `level`.

Interpret inputs precisely:

- `downstream_input` determines the first scoring layer and nonzero entries are appended as level `0` at the default cutoff;
- `upstream_input` only filters *overlapping scored nodes* whose inferred sign conflicts with the supplied sign. It does not force a node into the result or constrain propagation as a seed;
- `n_layers` is an upper limit. Sparse layers can stop early, so inspect `table(moon_res$level)`, maximum realized level, and overlap with expected upstream nodes;
- `ulm` produces a score from decoupleR’s univariate linear model. Do not label it a p-value. In the current implementation, `n_perm` is relevant to `norm_wmean`, ignored by ULM, and not honored by `wmean`.

Do not use `return_levels` as a behavior switch; the current implementation always returns `level`.

## Apply transcript-level coherence deliberately

The package helper requires:

```r
tf_reg_net  # source, target, mor
rna_input   # named numeric RNA values
```

It removes a TF-target edge when `sign(TF_score * RNA_target * mor) < 0`. Missing RNA, zero values, and `NA` values do not cause removal. The helper does not inspect the PKN edge sign, does not enforce protein support, and does not implement the manuscript’s pre-score consistency pass.

Use a pre-score check when data and the biological question support it. With matched proteomics, do not pretend that the current helper enforces protein-gated RNA evidence: explicitly pre-prune unsupported TF-target edges or implement a separately documented rule. Keep raw RNA evidence separate from TF activity scores.

Run a bounded post-score loop and retain the PKN used for every score:

```r
pkn_final <- pkn_scoring
converged <- FALSE

for (i in seq_len(max_iterations)) {
  pkn_used_for_score <- pkn_final
  moon_res <- cosmosR::moon(
    upstream_input = upstream_input,
    downstream_input = downstream_input,
    meta_network = pkn_used_for_score,
    n_layers = n_layers,
    statistic = "ulm"
  )
  pkn_final <- cosmosR::filter_incohrent_TF_target(
    moon_res, tf_reg_net, pkn_used_for_score, rna_input
  )

  if (identical(pkn_final, pkn_used_for_score)) {
    converged <- TRUE
    break
  }
}

if (!converged) {
  # Score the final PKN explicitly and report non-convergence.
  moon_res <- cosmosR::moon(
    upstream_input = upstream_input,
    downstream_input = downstream_input,
    meta_network = pkn_final,
    n_layers = n_layers,
    statistic = "ulm"
  )
}
```

Record the number and identities of removed edges per iteration. A stable edge count is only a convergence signal if the filter is monotonic; compare the actual edge table whenever practical.

## Validate before interpretation

Before extracting a network, verify all of the following:

- `moon_res` has exactly one row per scored node and contains `source`, `score`, and `level`;
- every required level-0 node has been retained according to the intended downstream cutoff;
- the realized depth and upstream input coverage are sufficient for the stated question;
- `pkn_final` is the same PKN whose scores will be explained;
- TF-target filtering was either not applicable, applied and converged, or stopped at a recorded cap;
- compressed virtual IDs have not been mixed with an uncompressed PKN.

Use `get_moon_scoring_network()` later to inspect the local evidence underlying one score. It returns an explanatory topology based on levels; it does not prove a unique causal path.
