# Finance Algorithm Functions

This is a quick overview of the functions used for carrying out CLA on the standard sum-to-one constrained portfolio problem that commonly occurs in finance.

## Utility Functions

These functions aren't actually part of the algorithm themselves, but they carry out calculations or operations the results of which the algorithm uses.

1. `max_bounded`: maximum value in $\boldsymbol{x}$ satisfying $x_i < x_{\text{bound}}$
2. `inverse_grow`: adjust the inverse if gaining a row and column
3. `inverse_shrink`: adjust the inverse if removing row and column $i$

## KKT Check Functions

These two functions carry out a KKT check on a solution to a problem with the corresponding input variables. They make use of a vector $\boldsymbol{d}$ calculated by the ALGORITHM FUNCTIONS that contain the value of the bound a given weight is moving towards (itself calculated using the derivative).

1. `kkt_full`: perform a full KKT check on a constrained CLA problem
2. `kkt_partitioned`: equivalent KKT checks for constrained CLA

## Algorithm Functions

These four functions make up the key part of the algorithm; finding the starting solution, handling the case when an asset moves to its bound, handling the case when an asset become free, and then calculating the turning points themselves through CLA.

1. `starting_solution`: return starting solution for CLA
2. `becomes_free`: handle the CLA case where an asset becomes free
3. `move_to_bound`: handle CLA case where an asset moves to its bound
4. `calculate_turningpoints`: return portfolios for CLA turning points
