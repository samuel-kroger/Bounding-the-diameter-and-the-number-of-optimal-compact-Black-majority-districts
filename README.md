# Bounding the number and the diameter of compact Black-majority districts
This GitHub repository contains all the code, data, and results used in the paper "[Bounding the number and the diameter of compact Black-majority districts](https://optimization-online.org/2024/08/bounding-the-number-and-the-diameter-of-optimal-compact-black-majority-districts/)," authored by Samuel Kroger, Hamidreza Validi, Illya V. Hicks, and Tyler Perini. 

This paper proposes mixed integer programming (MIP) formulations and techniques to explore the maximum number of Black-majority congressional districts for multiple states of the United States. Furthermore, we generalize the diameter-based compactness criterion of Garfinkel and Nemhauser (*Management Science*, 1970) and provide a framework for optimizers to capture compactness in constraints rather than the objective function. To alleviate the solving process, we propose fixing procedures and symmetry-breaking constraints. Our proposed MIP formulations provide (i) an upper bound on the number of Black-majority districts and (ii) lower bounds for the diameter of districts. We finally run a state-of-the-art districting package, GerryChain, to provide feasible bounds for the maximum number of Black-majority districts. This provides the best-known existing optimality gap that can be closed by both districters and optimizers in the future. 

The following figures show two county-level maps of MS that are obtained by our proposed MIP formulation (left) and 10,000 iterations of Gerry Chain's short bursts module (right).

A MIP-generated map of MS/county with one Black-majority and diameter of 6             |  A GerryChain-generated map of MS/county with one Black-majority district and diameter of 7
:-------------------------:|:-------------------------:
![](readme_images/MS_map_MIP_opt.png?raw=true "A MIP-generated map of MS/county with one Black-majority and diameter of 6")   |  ![](readme_images/short_bursts_opt.png?raw=true "A GerryChain-generated map of MS/county with one Black-majority district and diameter of 7")

# Requirements

To run the code, you will need installations of [Gurobi](https://www.gurobi.com/) and [GerryChain](https://gerrychain.readthedocs.io/en/latest/).

The input data is provided by [Daryl DeFord](https://www.math.wsu.edu/faculty/ddeford/).

# How to run the code yourself?

All of the code works off running comp_experiment.py and can be executed with the following line:
```
C:\src\comp_expermient.py
```
At the bottom of comp_experiment.py you will find all the run options. Uncomment a line in comp_experiment.py to generate the corresponding table.
