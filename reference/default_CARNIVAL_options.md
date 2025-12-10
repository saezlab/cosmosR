# Setting Default CARNIVAL Options

Returns the default CARNIVAL options as a list. You can modify the
elements of the list and then use it as an argument in
[`run_COSMOS_metabolism_to_signaling`](run_COSMOS_metabolism_to_signaling.md)
or
[`run_COSMOS_signaling_to_metabolism`](run_COSMOS_signaling_to_metabolism.md).
If you choose CPLEX or CBC, you must modify then the solverPath field
and point to the CPLEX/CBC executable (See Details).

## Usage

``` r
default_CARNIVAL_options(solver = NULL)
```

## Arguments

- solver:

  one of \`cplex\` (recommended, but require 3rd party tool), \`cbc\`
  (also require 3rd party tool) or \`lpSolve\` (only for small networks)

## Value

returns a list with all possible options implemented in CARNIVAL. see
the documentation on
[`runCARNIVAL`](https://rdrr.io/pkg/CARNIVAL/man/runCARNIVAL.html).

## Details

COSMOS is dependent on CARNIVAL for exhibiting the signalling pathway
optimisation. CARNIVAL requires the interactive version of IBM Cplex,
Gurobi or CBC-COIN solver as the network optimiser. The IBM ILOG Cplex
is freely available through Academic Initiative
[here](https://www.ibm.com/products/ilog-cplex-optimization-studio).
Gurobi license is also free for academics, request a license following
instructions
[here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).
The [CBC](https://projects.coin-or.org/Cbc) solver is open source and
freely available for any user, but has a significantly lower performance
than CPLEX or Gurobi. Obtain CBC executable directly usable for cosmos
[here](https://ampl.com/products/solvers/open-source/#cbc).
Alternatively for small networks, users can rely on the freely available
[lpSolve
R-package](https://cran.r-project.org/web/packages/lpSolve/index.html),
which is automatically installed with the package.

## Examples

``` r
# load and change default options: 
my_options = default_CARNIVAL_options(solver = "cplex")
 
my_options$solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
my_options$threads = 2
my_options$timelimit = 3600*15
```
