# FINANCE ALGORITHM FUNCTIONS

## UTILITY FUNCTIONS

These functions aren't actually part of the algorithm themselves, but they carry out calculations or operations the results of which the algorithm uses.

1. max_conditional
2. subindex

This code also depends on the following functions from the standard implementation:

- MAX_BOUNDED maximum value in x satisfying x(i) < x_bound
- INVERSE_SHRINK adjust the inverse if gaining a row and column
- INVERSE_SHRINK adjust the inverse if removing row and column i

## KKT CHECK FUNCTIONS

These two functions carry out a KKT check on a solution to a problem with the corresponding input variables. They make use of a vector d calculated by the ALGORITHM FUNCTIONS which contains the value of the bound a given weight is moving towards (itself calculated using the derivative).

1. kkt_full_gen
2. kkt_partitioned_gen

## ALGORITHM FUNCTIONS

These four functions make up the key part of the algorithm; finding the starting solution, handling the case when an asset moves to its bound, handling the case when an asset become free, and then calculating the turning points themselves through CLA.

1. starting_solution_gen
2. multiplier_update
3. move_to_bound_gen
4. becomes_free_gen
5. calculate_turningpoints_gen
