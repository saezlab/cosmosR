# PKN Preparation and Safe Compression

Prepare a separate, disposable working PKN. Preserve the raw PKN and its provenance unchanged.

## Normalize the PKN contract

`meta_network_cleanup()` is strict in practice. Give it only numeric signed edges:

```r
pkn <- raw_pkn[, c("source", "interaction", "target")]
pkn$source <- as.character(pkn$source)
pkn$target <- as.character(pkn$target)
pkn$interaction <- as.numeric(pkn$interaction)

stopifnot(
  !anyNA(pkn),
  all(pkn$interaction %in% c(-1, 1))
)

pkn <- cosmosR::meta_network_cleanup(pkn)
pkn <- pkn[order(pkn$source, pkn$target, pkn$interaction), ]
row.names(pkn) <- NULL
```

The cleaner removes self interactions, collapses duplicate source-target pairs by mean interaction, and drops pairs whose resulting interaction is not exactly `-1` or `1`. Do not pass annotations or character metadata: it averages all remaining columns.

Use `interaction` as the working sign column. `moon()` also accepts `sign`, but `keep_controllable_neighbours()`, `keep_observable_neighbours()`, and both solution reducers require `interaction`.

## Validate inputs before pruning

Require each input to be a named numeric vector with unique, non-missing, finite values. Do not silently repair duplicate names by averaging unless the data-generating transformation justifies that aggregation.

For metabolomics, map feature identifiers to the PKN’s HMDB identity before using `prepare_metab_inputs()`. That function only adds `Metab__` and one or more compartment suffixes; it does not map identifiers, validate coverage, or deduplicate names.

```r
metab_input <- cosmosR::prepare_metab_inputs(metab_input, c("c", "m"))
```

Valid compartments are `r`, `c`, `e`, `x`, `m`, `l`, `n`, and `g`. The function duplicates each input score across requested compartments. Record this expansion rather than interpreting it as independent measurements.

If RNA expression is an appropriate context filter, run `cosmosR:::filter_pkn_expressed_genes_fast(names(rna_input), pkn)`. It retains metabolite and orphan-reaction nodes by design. Do not use it to infer biological activity; it only removes unsupported gene-associated mechanisms.

## Prune and stabilize input coverage

Follow the generic COSMOS_basic order, working on a copy:

```r
n_steps <- 6L # A starting search radius, not a universal default.

upstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_input, pkn)
pkn <- cosmosR:::keep_controllable_neighbours(pkn, n_steps, names(upstream_input))

downstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_input, pkn)
pkn <- cosmosR:::keep_observable_neighbours(pkn, n_steps, names(downstream_input))

upstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_input, pkn)
downstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_input, pkn)
```

Repeat membership filtering until the retained input lengths are stable when a pruning step can remove nodes. Log the IDs removed at each step. The forward and reverse ego-neighborhood filters define a feasible region, but do not guarantee that every retained node lies on one directed upstream-to-downstream path. Their radii are conceptually separate from `moon()`’s `n_layers`, even if a tutorial reuses the same value.

Stop and report a coverage problem if either retained input set is biologically inadequate or empty. Do not continue with a numerically valid but scientifically unrelated network.

## Compress only with safeguards

Conceptually, `compress_same_children()` merges non-input parents with exactly the same direct signed children to prevent redundant paths from being repeatedly counted. If TF-target coherence will be used, protect TF-regulon sources and targets as well as measured inputs. Values supplied only for protection are irrelevant; the function uses their names.

```r
protected_parents <- unique(c(
  names(upstream_input), names(downstream_input),
  tf_reg_net$source, tf_reg_net$target
))
protection_vector <- setNames(rep(0, length(protected_parents)), protected_parents)

# FALSE sorts first, so protected parents occur after non-protected parents.
pkn <- pkn[order(
  pkn$source %in% protected_parents,
  pkn$source, pkn$target, pkn$interaction
), ]
row.names(pkn) <- NULL
```

Keep the full compression return value:

```r
compressed <- cosmosR::compress_same_children(
  pkn,
  sig_input = protection_vector,
  metab_input = downstream_input
)
pkn_scoring <- cosmosR::meta_network_cleanup(compressed$compressed_network)
```

The current implementation serializes child edges in their existing row order and can rename a protected seed if it appears before an otherwise identical non-seed parent. The deterministic order above makes groups that include a protected parent remain uncompressed, while allowing safe groups without protected parents to compress.

After compression, require all protected IDs to remain unchanged and require none to occur in `names(compressed$duplicated_signatures)`. If either check fails, do not use the compressed PKN. Use the stable uncompressed PKN instead and record that the safety check prevented compression.

Keep these objects together:

- `pkn_uncompressed`: cleaned, expression-filtered, pruned PKN in original IDs;
- `pkn_scoring`: the actual compressed PKN, or `pkn_uncompressed` when compression is skipped;
- `compressed`: compression mapping, only when compression passed checks;
- the exact input vectors retained after the final membership filter.

Never discard `pkn_uncompressed`: it is needed for original-ID output and audit.
