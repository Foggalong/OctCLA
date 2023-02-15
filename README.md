# OctCLA

Implementations of Markowitz' Critical Line Algorithm in [GNU Octave](https://www.gnu.org/software/octave/) with an example. Code is MATLAB compatible and released under the [MIT License](LICENSE).

## Finance

The finance codebase explores the standard portfolio optimization problem

$$\min w^T\mu,\ \max w^T\Sigma w\ \text{ subject to }\ w\geq0,\ w^Te = 1,\ l\leq w\leq u,$$

where $w$ is the vector of weights, $\Sigma$ is the covariance matrix, $\mu$ is the vector of expected returns, $e$ is the vector of ones, $l$ is a lower bound on weights, and $u$ is an upper bound on weights.

This problem is studied most frequently in the context of financial investment where each $w_i$ is the amount invested in asset $i$. The $w\geq0$ corresponds to a block on short selling in that context.

An overview of the functions involved can be found in [`finance/`](finance/).

## Genetics

The genetics codebase explores the related portfolio optimization problem

$$\min w^T\mu,\ \max w^T\Sigma w\ \text{ subject to }\ w\geq0,\ w^T_{\mathcal{S}}e_{\mathcal{S}}^{\phantom{T}} = \frac{1}{2}, w^T_{\mathcal{D}}e_{\mathcal{D}}^{\phantom{T}} = \frac{1}{2},\ l\leq w\leq u,$$

where $\mathcal{S}$ and $\mathcal{D}$ are two index sets which partition the complete set of asset indices.

This problem arises in genetics in the context of managing breeding programs. The two sum-to-half constraints ensure that half our genetic product is coming from male candidates (sires) and half from female candidates (dams).

While some of the required functions are reused from [`finance/`](finance/), many have differences (but are similar) and others are completely unrelated. An overview of all the functions involved can be found in [`genetics/`](genetics/).

## References

The following works proved useful when creating this implementation.

1. H. M. Markowitz and G. P. Todd, "Portfolio Choice: Program Description," in _Mean-Variance Analysis in Portfolio Choice and Capital Markets_, pp. 302–336, John Wiley & Sons, Feb. 2000.
2. A. Niedermayer and D. Niedermayer, "Applying Markowitz’s Critical Line Algorithm," in _Handbook of Portfolio Construction_, pp. 383–400, Springer, Jan. 2010.
3. D. Bailey and M. López de Prado, "An Open-Source Implementation of the
Critical-Line Algorithm for Portfolio Optimization," _Algorithms_, vol. 6, pp. 169–
196, Mar. 2013.
