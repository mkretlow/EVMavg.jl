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
println("EVM average: $(evm_result.mean) ± $(evm_result.external_error)")
println("Confidence: $(round(evm_result.confidence_level * 100, digits=1))%")

# EVM automatically down-weights outliers and provides realistic uncertainties
# Weighted average: 5.93 ± 0.48
# EVM average: 5.89 ± 0.87  (more realistic given the scatter)
# Confidence: 65.4%
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

Always report `upper_error` and `lower_error` - these are the final uncertainties that 
conservatively combine both sources:

- **Internal uncertainties**: Reflect measurement precision (from input errors)
- **External uncertainty**: Reflects measurement scatter (consistency of data)
- **Final uncertainties**: Take max(internal, external) variance for each direction

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

**Which uncertainty to report?**

- Use **external error** when measurements show scatter (typical case)
- Use **internal error** only when measurements are highly consistent
- Generally, report the larger of the two uncertainties

### Weighted Average

For comparison, the package also provides traditional weighted averaging:

```julia
result = avg(data, sigmas)
```

**Returns:** `AvgResult` containing:

- `unweighted`: Simple arithmetic mean
- `weighted`: Inverse-variance weighted average
- `std_error`: Standard error of the weighted average

## Advantages Over Traditional Methods

### 1. More Realistic Uncertainties

Traditional weighted averaging can significantly underestimate uncertainties when data disagree:

```julia
# Discrepant data
data = [10.2, 10.5, 14.3, 10.1, 10.4]  # One clear outlier
errors = [0.2, 0.2, 0.3, 0.2, 0.2]

evm_result = evm(data, errors)
avg_result = avg(data, errors)

# Unweighted average: 11.1
# Weighted average  : 10.7 ± 0.8 (unrealistic - ignores scatter)
# EVM average       : 10.5 ± 1.0 (realistic - reflects actual spread)
# W.avg, w/o outlier: 10.3
```

The EVM external uncertainty reflects the actual scatter in the data, providing a more realistic estimate even when data are discrepant.

### 2. Robust to Outliers

Outlying data points automatically receive reduced weight without manual intervention:

```julia
println("EVM weights: ", evm_result.weights)
# Outlier at 14.3 will have much lower weight than other data points
```

### 3. Interpreting Confidence Levels
The confidence level is a goodness-of-fit measure, but it's specifically testing whether the data are consistent with the EVM model (the mean probability density M(x)). Interpretation:  

**High confidence** (e.g., 0.8-0.95 or 80-95%):

- The data points are distributed around the EVM mean in a way that's consistent with the expected probability distribution
- The data "look like" they could have come from the model
- Good data quality and consistency

**Medium confidence** (e.g., 0.5-0.8 or 50-80%):

- Acceptable fit, some scatter but not alarming
- Data are reasonably consistent

**Low confidence** (e.g., <0.5 or <50%):

- Data distribution is unusual given the model
- Could indicate: outliers, systematic errors, underestimated uncertainties, or that measurements aren't actually sampling the same underlying value

**Very high confidence** (e.g., >0.95 or >95%):

- Could be suspicious! May indicate overestimated uncertainties or too few data points for the test to be reliable


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

### Uncertainties

**Internal uncertainty** (from measurement precision):
```
σ_int = √(Σᵢ wᵢ² σᵢ²)
```

**External uncertainty** (from measurement spread):
```
σ_ext = √(Σᵢ wᵢ (μᵢ - x_EVM)²)
```

### Confidence Level

Uses a modified χ² test comparing observed vs. expected numbers of measurements above/below the mean. Values near 0 or 1 may indicate issues; typical good fits have 0.5-0.95 confidence.

The Statistical Test:  

- Count how many points fall above vs. below the EVM mean
- Compare to the expected ratio (based on the mean probability density)
- Compute a χ² statistic Q
- Convert to confidence level using 1 - erf(√(Q/2))

## Best Practices

- **Report external error**: Unless data is extremely consistent, external error is more realistic
- **Inspect weights**: Large weight disparities may indicate outliers
- **Use confidence level**: Values < 0.5 suggest investigating data quality
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