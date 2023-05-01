# CS 520 Spring 23 - Project Choice #2

Authors: Yunhan Huang and Jiasong Dong

Primal-Dual Interior-Point Method Solver for LP Problems

## Run test(s)

After cloning the github repository, run:

`julia run-test.jl`

to run the tests required by professor on course website ("lp_afiro","lp_brandy","lp_fit1d","lp_adlittle","lp_agg","lp_ganges","lp_stocfor1", "lp_25fv47", "lpi_chemcom")

Or, to run solver on one problem:

`julia run-test.jl [problem name]`

where [problem name] is the problem name on [LPNetLib](https://netlib.org/lp/data/index.html), no "LPNetLib/" prefix needed.

*We also provide the jupyter notebook here, if that is prefered*

## Comparing Results

The test compares result with the result of GLPK optimizer, which sometimes is inaccurate or incorrect since there may be more parameters that need to be configurate. The actual result of these LP problems can be found at [Problem Summary](https://netlib.org/lp/data/readme) provided by NetLib.ORG (we also use the summary here to compare with in our report).
