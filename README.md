# OctCLA

Two implementations of algorithms that explore the Pareto frontier of given bi-objective optimization problems, both building on Markowitz' Critical Line Algorithm. These are written in [GNU Octave](https://www.gnu.org/software/octave/) and are MATLAB compatible. Code is released under the [MIT License](LICENSE) and example usage is included.

## Finance

The "finance" implementation explores the standard portfolio optimization problem

$$
    \min_{\boldsymbol{w}} \left(\boldsymbol{w}^T\boldsymbol{\mu}\right),\ \max_{\boldsymbol{w}} \left(\boldsymbol{w}^T\Sigma\boldsymbol{w}\right)\ \text{ subject to }\ \boldsymbol{w}\geq0,\ \boldsymbol{w}^T\boldsymbol{e} = 1,\ \boldsymbol{l}\leq \boldsymbol{w}\leq \boldsymbol{u},
$$

where $\boldsymbol{w}$ is the vector of weights, $\Sigma$ is the covariance matrix, $\boldsymbol{\mu}$ is the vector of expected returns, $\boldsymbol{e}$ is the vector of ones, $\boldsymbol{l}$ is a lower bound on weights, and $\boldsymbol{u}$ is an upper bound on weights.

This problem is studied in the context of financial investment where each $w_i$ is the amount invested in asset $i$. The $\boldsymbol{w}\geq0$ corresponds to a block on short selling in that context.

An overview of the functions involved can be found in [`finance/`](finance/README.md).

## Genetics

The "genetics" implementation explores a similar optimal contribution selection (OCS) problem

$$
    \min_{\boldsymbol{w}} \left(\boldsymbol{w}^T\boldsymbol{\mu}\right),\ \max_{\boldsymbol{w}} \left(\boldsymbol{w}^T\Sigma\boldsymbol{w}\right)\ \text{ subject to }\ \boldsymbol{w}\geq0,\ \boldsymbol{w}^T_{\mathcal{S}}\boldsymbol{e}_{\mathcal{S}}^{} = \frac{1}{2}, \boldsymbol{w}^T_{\mathcal{D}}\boldsymbol{e}_{\mathcal{D}}^{} = \frac{1}{2},\ \boldsymbol{l}\leq \boldsymbol{w}\leq \boldsymbol{u},
$$

where $\mathcal{S}$ and $\mathcal{D}$ are two index sets which partition the complete set of asset indices.

OCS problems arise in the context of breeding programme management. The two sum-to-half constraints ensure that half our genetic product is coming from male candidates (sires) and half from female candidates (dams).

While some of the required functions are reused from [`finance/`](finance/), many have differences (but are similar) and others are completely unrelated. An overview of all the functions involved can be found in [`genetics/`](genetics/README.md).

## References

This code was written as an integral part of Josh Fogg's thesis, _New Frontiers for CLA_, submitted to the University of Edinburgh. Each implementation is woven into a didactic for deriving and implementing a CLA-like algorithm for the corresponding problem.

The following works were instrumental when creating this implementation.

1. Harry M. Markowitz and G. Peter Todd, "Portfolio Choice: Program Description," in _Mean-Variance Analysis in Portfolio Choice and Capital Markets_, pp. 302-336, John Wiley & Sons, Feb. 2000.
2. Andras Niedermayer and Daniel Niedermayer, "Applying Markowitz’s Critical Line Algorithm," in _Handbook of Portfolio Construction_, pp. 383-400, Springer, Jan. 2010.
3. David Bailey and Marcos López de Prado, "An Open-Source Implementation of the Critical-Line Algorithm for Portfolio Optimization," _Algorithms_, vol. 6, pp. 169-196, Mar. 2013.
