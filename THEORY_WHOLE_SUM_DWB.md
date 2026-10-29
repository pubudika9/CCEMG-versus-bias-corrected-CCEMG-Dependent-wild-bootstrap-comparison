# Whole-sum multiplier bootstrap with $\ell_N=1$

## 1. Whole-sum representation

After analytical bias correction, the leading expansion can be written as

$$
\sqrt{NT}(\widehat\beta_{ABC}-\beta_0)
=\frac1{\sqrt N}\sum_{i=1}^N Z_{i,NT}+o_p(1),
$$

where

$$
Z_{i,NT}
=\sqrt{\frac TN}v_i
+\frac1{\sqrt T}\sigma_{v,i}^{-1}V_i'\varepsilon_i.
$$

Thus the coefficient-heterogeneity and observation-error components are
combined in one unit-level summand. Multipliers must be applied to the
cross-sectional index $i$.

## 2. Feasible corrected individual slopes

The implementation uses the subtraction sign in Equation (44):

$$
\widehat\beta_{ABC}
=\widehat\beta_{MG}
-\frac{\widehat b_0+\widehat b_1}{\sqrt{NT}}.
$$

Write the estimated aggregate corrections as averages of unit terms:

$$
\widehat b_0=\sqrt{\frac TN}\frac1N\sum_i\widehat b_{0i},
\qquad
\widehat b_1=\sqrt{\frac TN}\frac1N\sum_i\widehat b_{1i}.
$$

Define

$$
\widehat\beta_{ABC,i}
=\widehat\beta_i
-\frac{\widehat b_{0i}+\widehat b_{1i}}N.
$$

Then

$$
\widehat\beta_{ABC}
=\frac1N\sum_i\widehat\beta_{ABC,i}.
$$

The feasible centered unit contribution is

$$
\widehat Z_{i,NT}
=\sqrt T
\left(\widehat\beta_{ABC,i}-\widehat\beta_{ABC}\right).
$$

For vector slopes, all quantities are replaced by their vector counterparts.

## 3. Finite-sample correction

The code applies

$$
\widehat Z_{i,NT}^{c}
=\sqrt{\frac N{N-1}}\widehat Z_{i,NT}.
$$

This makes the $\ell_N=1$ variance coincide exactly with the usual
finite-sample-corrected mean-group covariance formula applied to the corrected
individual slopes.

## 4. Independent multiplier bootstrap

Fix $\ell_N=1$ and generate

$$
\xi_i^*\stackrel{iid}{\sim}N(0,1)
$$

conditionally on the simulated panel. The bootstrap statistic is

$$
R_{NT}^*
=\frac1{\sqrt N}\sum_{i=1}^N
\widehat Z_{i,NT}^{c}\xi_i^*.
$$

Its exact conditional covariance is

$$
\widehat\Sigma_{DWB}
=\frac1N\sum_{i=1}^N
\widehat Z_{i,NT}^{c}
\widehat Z_{i,NT}^{c\prime}.
$$

For scalar slopes,

$$
\widehat\Sigma_{DWB}
=\frac T{N-1}\sum_{i=1}^N
(\widehat\beta_{ABC,i}-\widehat\beta_{ABC})^2.
$$

Consequently,

$$
\widehat{\operatorname{Var}}(\widehat\beta_{ABC})
=\frac{\widehat\Sigma_{DWB}}{NT}
=\frac1{N(N-1)}
\sum_{i=1}^N
(\widehat\beta_{ABC,i}-\widehat\beta_{ABC})^2.
$$

The reported ABC standard error is

$$
\widehat{se}(\widehat\beta_{ABC})
=\sqrt{\frac{\widehat\Sigma_{DWB}}{NT}}.
$$

The sample variance of the simulated multiplier statistics is retained only as
a numerical diagnostic; it is not used in the reported standard error.

## 5. Asymmetric basic interval

Let $q_{\alpha/2}^*$ and $q_{1-\alpha/2}^*$ be the two empirical tail
quantiles of the bootstrap statistics. The interval is

$$
\left[
\widehat\beta_{ABC}
-\frac{q_{1-\alpha/2}^*}{\sqrt{NT}},
\quad
\widehat\beta_{ABC}
-\frac{q_{\alpha/2}^*}{\sqrt{NT}}
\right].
$$

Thus the upper-tail bootstrap quantile determines the lower endpoint, and the
lower-tail quantile determines the upper endpoint. The implementation evaluates
the two tails separately.

## 6. Standard CCEMG comparison

The standard mean-group covariance is

$$
\widehat{\operatorname{Var}}(\widehat\beta_{MG})
=\frac1{N(N-1)}
\sum_{i=1}^N
(\widehat\beta_i-\widehat\beta_{MG})^2.
$$

Its normal confidence interval is

$$
\left[
\widehat\beta_{MG}-z_{1-\alpha/2}\widehat{se}_{MG},
\quad
\widehat\beta_{MG}+z_{1-\alpha/2}\widehat{se}_{MG}
\right].
$$

## 7. Validity conditions

The $\ell_N=1$ implementation is appropriate under the working paper's
cross-sectional independence conditions after common-factor removal. A
representative sufficient set of requirements is:

1. the displayed unit-level linear representation holds;
2. the unit influences have uniformly bounded $2+\delta$ moments;
3. the leading unit influences are independent across $i$;
4. a Lindeberg condition holds;
5. $N^{-1}\sum_i E[Z_{i,NT}Z_{i,NT}']\to\Sigma_1$;
6. replacing the infeasible influences by the centered corrected individual
   slopes is valid in the diagonal quadratic form.

Under these conditions,

$$
\widehat\Sigma_{DWB}\xrightarrow{p}\Sigma_1
$$

and the conditional multiplier distribution consistently estimates the
distribution of
$\sqrt{NT}(\widehat\beta_{ABC}-\beta_0)$.

If cross-sectional dependence remains after CCE removal, independent
multipliers are not sufficient. A spatial, network or otherwise substantively
ordered dependence-respecting bootstrap would then be required.

## 8. Heterogeneity exponent

For

$$
\beta_i=\beta_0+N^{-\alpha}\eta_i,
$$

the unit influence contains

$$
\sqrt T N^{-\alpha}\eta_i.
$$

When $T/N\to\rho\in(0,\infty)$, the localized $\sqrt{NT}$ limit is
nondegenerate at $\alpha=1/2$. The $\alpha=0$ experiments remain useful
robustness exercises, but they are not direct simulations of the fixed-variance
localized limit.
