# FINANCE ALGORITHM FUNCTIONS

## UTILITY FUNCTIONS

These functions aren't actually part of the algorithm themselves, but they carry out calculations or operations the results of which the algorithm uses.

1. max_bounded
2. inverse_grow
3. inverse_shrink

## KKT CHECK FUNCTIONS

These two functions carry out a KKT check on a solution to a problem with the corresponding input variables. They make use of a vector d calculated by the ALGORITHM FUNCTIONS which contains the value of the bound a given weight is moving towards (itself calculated using the derivative).

1. kkt_full
2. kkt_partitioned

## ALGORITHM FUNCTIONS

These four functions make up the key part of the algorithm; finding the starting solution, handling the case when an asset moves to its bound, handling the case when an asset become free, and then calculating the turning points themselves through CLA.

1. starting_solution
2. becomes_free
3. move_to_bound
4. calculate_turningpoints
