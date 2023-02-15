# Genetics Algorithm Functions

This is a quick overview of the functions used for carrying out CLA on the standard two-sum-to-half constrained portfolio problem which most commonly occurs in genetics.

## Utility Functions

These functions aren't actually part of the algorithm themselves, but they carry out calculations or operations the results of which the algorithm uses.

1. `max_conditional`: maximum value in $x$ satisfying $x_i < x_{\text{bound}}$ and $i\in X$.
2. `subindex`: given $X$, $Y$ which index $V$, returns $Z = \lbrace j : Y_j \in X \rbrace$.

This code also depends on the following functions from the implementation in the finance code base:

1. `max_bounded`: maximum value in $x$ satisfying $x_i < x_{\text{bound}}$
2. `inverse_grow`: adjust the inverse if gaining a row and column
3. `inverse_shrink`: adjust the inverse if removing row and column $i$

## KKT Check Functions

These two functions carry out a KKT check on a solution to a problem with the corresponding input variables. They make use of a vector d calculated by the ALGORITHM FUNCTIONS which contains the value of the bound a given weight is moving towards (itself calculated using the derivative).

1. `kkt_full_gen`: perform a full KKT check on a genetics CLA problem
2. `kkt_partitioned_gen`: equivalent KKT checks for genetics CLA

## Algorithm Functions

These four functions make up the key part of the algorithm; finding the starting solution, handling the case when an asset moves to its bound, handling the case when an asset become free, and then calculating the turning points themselves through CLA.

1. `starting_solution_gen`: return starting solution for CLA genetics problems
2. `multiplier_update`: updating of lagrangian multipliers ($\lambda$, $\gamma$, and $\delta$)
3. `move_to_bound_gen`: handle the genetics CLA case where an asset moves to its bound
4. `becomes_free_gen`: handle the genetics CLA case where an asset becomes free
5. `calculate_turningpoints_gen`: return portfolio weights at genetics CLA turning points
