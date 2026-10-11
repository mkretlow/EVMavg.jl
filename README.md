# EVMavg.jl

[![Build Status](https://github.com/mkretlow/EVMavg.jl/workflows/CI/badge.svg)](https://github.com/mkretlow/EVMavg.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia package for robust averaging of numerical data using the **Expected Value Method (EVM)** - an alternative to traditional weighted and unweighted averaging methods.

## Why EVM?

Traditional averaging methods (unweighted mean, inverse-variance weighted average) can be significantly affected by outliers and may underestimate uncertainties when data are inconsistent. The Expected Value Method addresses these limitations by:

- **Automatically handling outliers**: Discrepant data points receive reduced weight without manual intervention
- **Providing realistic uncertainties**: Both internal (precision-based) and external (scatter-based) estimates
- **Supporting asymmetric errors**: Native handling of data with different upper (+) and lower (-) uncertainties
- **Assessing data quality**: Built-in goodness-of-fit measures

## When to Use EVM

**EVM is ideal for:**

- Any scenario with potential outliers or inconsistent data
- Data with asymmetric uncertainties (from likelihood fits, Bayesian analysis, bootstrapping)
- Situations where you need conservative, realistic uncertainty estimates

**Stick with traditional methods when:**

- All data points are highly consistent (no outliers)
- You need simple, fast computation for very large datasets
- Only symmetric uncertainties are present and data quality is high

## Features

- **Three averaging methods**: Traditional weighted/unweighted mean, symmetric EVM, and asymmetric EVM
- **Realistic Uncertainty Estimates**: Provides both internal (precision-based) and external (scatter-based) uncertainties
- **Asymmetric Uncertainty Support**: Handles data with different upper and lower error bars
- **Robust to Outliers**: Outlying data points receive reduced weight automatically
- **Goodness of Fit**: Includes confidence level calculation to assess data consistency

## Installation

```julia
using Pkg
Pkg.add("EVMavg")
```

Or in the Julia REPL package mode (press `]`):

```julia
pkg> add EVMavg
```

## Quick Start

```julia
using EVMavg

# Example:  
data = [5.7, 6.7, 4.3, 6.0, 5.3, 6.16]
sigmas = [0.3, 0.6, 2.0, 1.3, 0.8, 0.34]

# Traditional weighted average
avg_result = avg(data, sigmas)
println("Weighted average: $(avg_result.weighted) ± $(avg_result.std_error)")

# EVM average (more robust)
evm_result = evm(data, sigmas)
σ = max(evm_result.internal_error, evm_result.external_error)
println("EVM average: $(evm_result.mean) ± $σ")
println("Confidence: $(round(evm_result.confidence_level * 100, digits=1))%")

# Weighted average: 5.94 ± 0.48
# EVM average: 5.88 ± 0.52  (internal 0.37, external 0.52)
# Confidence: 97.7%
```

## Methods

### Weighted Average (Traditional Method)

For comparison and cases where traditional methods are appropriate:

```julia
result = avg(data, sigmas)
```

**Returns:** `AvgResult` containing:

- `unweighted`: Simple arithmetic mean
- `weighted`: Inverse-variance weighted average
- `std_error`: Standard error of the weighted average

### Expected Value Method (EVM)

The EVM constructs a mean probability density from all measurements and computes a weighted average that balances measurement precision with data consistency.

#### Symmetric Uncertainties

```julia
result = evm(data, errors; compute_confidence=true)
```

**Arguments:**

- `data`: Vector of measured values
- `errors`: Vector of uncertainties (standard deviations)
- `compute_confidence`: Whether to calculate goodness of fit (default: `true`)

**Returns:** `EVMResult` containing:

- `mean`: EVM average value
- `internal_error`: Uncertainty from measurement precision
- `external_error`: Uncertainty from measurement spread
- `weights`: Normalized weight for each measurement
- `confidence_level`: Goodness of fit measure (0-1)

**Which uncertainty to report?**

Report the larger of `internal_error` and `external_error`. Following Birch & Singh (2014),
the quoted uncertainty should be the more conservative of the two estimates:

- **Internal uncertainty**: Reflects measurement precision (from the input errors)
- **External uncertainty**: Reflects measurement scatter (consistency of the data)

#### Asymmetric Uncertainties

For measurements with different upper and lower error bars (e.g., from likelihood fits):

```julia
result = evm_asymmetric(data, lower_errors, upper_errors; compute_confidence=true)
```

**Arguments:**

- `data`: Vector of measured values
- `lower_errors`: Vector of lower (negative direction) uncertainties
- `upper_errors`: Vector of upper (positive direction) uncertainties
- `compute_confidence`: Whether to calculate goodness of fit (default: `true`)

**Returns:** `EVMAsymmetricResult` containing:

- `mean`: EVM average value
- `lower_error`: Final lower uncertainty to report
- `upper_error`: Final upper uncertainty to report
- `internal_lower`: Internal lower uncertainty (diagnostic)
- `internal_upper`: Internal upper uncertainty (diagnostic)
- `external_error`: External uncertainty, symmetric (diagnostic)
- `weights`: Normalized weight for each measurement
- `confidence_level`: Goodness of fit measure (0-1)

**Which uncertainty to report?**

Always report `upper_error` and `lower_error` - these are the final uncertainties, chosen
conservatively from the two sources:

- **Internal uncertainties**: Reflect measurement precision (from input errors), separately for each direction
- **External uncertainty**: Reflects measurement scatter (consistency of data), symmetric
- **Final uncertainties**: Following Birch & Singh (2014), the estimate whose associated
  (asymmetric) normal distribution has the larger variance. The variance of the internal split
  normal, (1 − 2/π)(σ_int+ − σ_int-)² + σ_int+ σ_int-, is compared with σ_ext². Either the
  internal pair (+σ_int+ −σ_int-) or the symmetric external uncertainty (±σ_ext) is reported;
  the two are never mixed.

The other fields are provided for diagnostic purposes to understand uncertainty sources.

**Example:**

```julia
# Measurements with asymmetric errors
data = [10.2, 10.5, 10.1]
lower_errors = [0.3, 0.4, 0.2]  # -σ
upper_errors = [0.5, 0.6, 0.4]  # +σ

result = evm_asymmetric(data, lower_errors, upper_errors)
println("Result: $(result.mean) +$(result.upper_error) -$(result.lower_error)")

# Result: 10.25635530596962 +0.2885762483615959 -0.1752029094931697
```

## Advantages Over Traditional Methods

### 1. More Realistic Uncertainties

When data disagree, the traditional weighted average is pulled toward the discrepant point, and its
internal uncertainty (1/√Σwᵢ) ignores the scatter entirely:

```julia
# Discrepant data
data = [10.2, 10.5, 14.3, 10.1, 10.4]  # One clear outlier
errors = [0.2, 0.2, 0.3, 0.2, 0.2]

evm_result = evm(data, errors)
avg_result = avg(data, errors)

# Unweighted average: 11.1
# Weighted average  : 10.7  (internal uncertainty 1/√Σwᵢ = 0.09)
# EVM average       : 10.5 ± 1.0 (external uncertainty, reflects the spread)
# W.avg, w/o outlier: 10.3
```

The EVM external uncertainty reflects the actual scatter in the data, providing a more realistic estimate even when data are discrepant.

### 2. Robust to Outliers

Outlying data points automatically receive reduced weight without manual intervention:

```julia
println("EVM weights: ", round.(evm_result.weights, digits=3))
# EVM weights: [0.256, 0.213, 0.061, 0.213, 0.256]
# The outlier at 14.3 receives a much lower weight than the other data points
```

### 3. Interpreting Confidence Levels

The confidence level α tests whether the data are consistent with the EVM model, i.e. with
the mean probability density M(x). It is not related to the confidence interval of the
mean, which is always one standard deviation (68%).

- **Low values** (a few percent or less): the data are unlikely to be realizations of M(x).
  This may indicate outliers, systematic errors, underestimated uncertainties, or
  measurements that do not sample the same underlying value.
- **Moderate to high values**: no evidence against the model.
- **Very high values** (close to 100%): be cautious. With few data points the χ²
  approximation behind α is poor and the test may simply not be able to draw a conclusion;
  very high values can also hint at overestimated uncertainties.

The test only compares the number of points below and above the mean with the expected
fractions, so it is a fairly weak test, especially for small data sets.

## Method Details

### Mathematical Background

For each measurement (μᵢ ± σᵢ), the EVM assumes a Gaussian probability density:

```
p(x; μᵢ, σᵢ) = (1/√(2πσᵢ²)) exp(-(x-μᵢ)²/(2σᵢ²))
```

The mean probability density is:

```
M(x) = (1/n) Σᵢ p(x; μᵢ, σᵢ)
```

The EVM average is computed as:

```
x_EVM = Σᵢ wᵢ μᵢ,  where wᵢ = M(μᵢ) / Σⱼ M(μⱼ)
```

For asymmetric uncertainties (μᵢ +aᵢ −bᵢ), p is replaced by a split Gaussian (Eq. 1 of the paper):

```
A(x; μ, a, b) = √(2/(π(a+b)²)) exp(-(x-μ)²/(2b²)),  x ≤ μ
A(x; μ, a, b) = √(2/(π(a+b)²)) exp(-(x-μ)²/(2a²)),  x > μ
```

### Uncertainties

**Internal uncertainty** (from measurement precision):
```
σ_int = √(Σᵢ wᵢ² σᵢ²)
```

For asymmetric uncertainties: σ_int+ = √(Σᵢ wᵢ² aᵢ²) and σ_int- = √(Σᵢ wᵢ² bᵢ²).

**External uncertainty** (from measurement spread):
```
σ_ext = √(Σᵢ wᵢ (μᵢ - x_EVM)²)
```

### Confidence Level

A modified χ² test (Eqs. 6 and 7 of the paper) compares the observed numbers of
measurements below (N_l) and above (N_h) the mean with the expected ones:

```
Q = (N_l - n p_l)² / (n p_l) + (N_h - n p_h)² / (n p_h)
p_l = ∫_{-∞}^{x_EVM} M(x) dx,   p_h = 1 - p_l
α = 1 - erf(√(Q/2))
```

p_l is computed exactly from the cumulative distribution functions of the individual
(split) Gaussians. For large n, Q approaches a χ² distribution with one degree of freedom,
so α is only approximate for small data sets.

### Validation

The test suite reproduces the EVM results of Table I in Birch & Singh (2014) for the
first 2⁺ state in ⁴⁸Ti: τ = 5.83(86) ps for the lifetimes and B(E2)↑ = 0.075(11) e²b².

## Best Practices

- **Report the larger uncertainty**: Quote the larger of the internal and external uncertainty
- **Inspect weights**: Large weight disparities may indicate outliers
- **Use confidence level**: Low values suggest investigating data quality
- **Remove systematic outliers**: Suspect measurements should be excluded before averaging


## References

- Birch, M., Singh, B., 2014. Method of Best Representation for Averages in Data Evaluation. *Nuclear Data Sheets*, **120**, 106-108. DOI: [10.1016/j.nds.2014.07.019](https://doi.org/10.1016/j.nds.2014.07.019)


## Contributing

Contributions are welcome! Please feel free to submit a pull request. For major changes, please open an issue first to discuss what you would like to change.

## License

This package is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@article{birch2014method,
  title={Method of Best Representation for Averages in Data Evaluation},
  author={Birch, M. and Singh, B.},
  journal={Nuclear Data Sheets},
  volume={120},
  pages={106--108},
  year={2014},
  doi={10.1016/j.nds.2014.07.019}
}
```