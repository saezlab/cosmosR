# Manuscript-Style Pathway Control Analysis

The COSMOS+ manuscript describes pathway control as Algorithm 7: for each high-scoring controller, test whether pathways are over-represented among its downstream reachable gene nodes. This is not an exported `cosmosR` MOON function and is not part of the generic COSMOS_basic pipeline.

## Define the analysis honestly

Use a final, original-ID, score-consistent solution network. For each controller with a predeclared score threshold:

1. collect nodes reachable downstream within a predeclared number of directed steps;
2. remove the controller itself, non-gene nodes, metabolites, virtual nodes, and IDs absent from the pathway annotation;
3. intersect the result with an explicit eligible-gene background from the matched filtered PKN;
4. run over-representation analysis for a declared pathway resource;
5. record effect size, raw p-value, adjusted p-value, controller, pathway, target count, overlap, universe, and step limit;
6. reshape results to a controller-by-pathway matrix only after retaining the long-form table.

The manuscript used `keep_controllable_neighbours()` with a two-step limit and `piano::runGSAhyper()` with NABA and KEGG collections. The current package does not export this complete workflow, and `piano` is not a package dependency. Select and document an external ORA implementation rather than claiming cosmosR performed the test.

## Build a defensible universe

Do not use all human genes or all pathway-resource genes by default. Define the universe as the gene nodes eligible to appear downstream in the same cleaned, expression-filtered, pruned PKN and that are representable in the pathway resource. Keep the identifier conversion used to construct it.

Do not use `extract_nodes_for_ORA()` directly: it expects a classic CARNIVAL SIF/ATT schema, not current MOON outputs.

## Interpret the result conservatively

A significant controller-pathway association means the pathway is enriched among a graph-defined downstream neighborhood under the chosen PKN and background. It does not establish that the controller regulates the pathway in the sample.

Control false discoveries across the full controller-by-pathway family, not separately within each controller unless that is the predeclared scientific question. Report pathway size, overlap, step limit, and network version alongside adjusted results. Sensitivity-check meaningful findings across score and neighborhood thresholds.
