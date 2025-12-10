# Toy Input Network

This signaling network is the reduced COSMOS network solution obtained
in the cosmos test on 786-O NCI60 data. Here, this network solution is
reused as an exemplary input prior knowledge network (PKN).

## Usage

``` r
data(toy_network)
```

## Format

An object of class “`data.frame`” with 19 rows (interactions) and three
variables:

- `source`:

  Source node, either metabolite or protein

- `interaction`:

  Type of interaction, 1 = Activation, -1 = Inhibition

- `target`:

  Target node, either metabolite or protein

## Source

The network data are available on github:
<https://github.com/saezlab/COSMOS_MSB/tree/main/results/COSMOS_result/COSMOS_res_session.RData>.
The toy_network is the combined network of the COSMOS network solutions
CARNIVAL_Result2 and CARNIVAL_Result_rerun subsequently reduced to 19
exemplary nodes.

## References

Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) *Molecular
Systems Biology*. **17**, e9730.

## Examples

``` r
data(toy_network)
```
