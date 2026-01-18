"""
    EVMavg

A Julia package for robust averaging of numerical data using the 
Expected Value Method (EVM), an alternative to traditional weighted 
and unweighted averaging.

Unlike conventional averaging methods, EVM provides:
- **Robust outlier handling**: Outlying data points are automatically down-weighted
- **Realistic uncertainties**: Both internal (precision) and external (scatter) estimates
- **Asymmetric error support**: Native handling of data with different upper/lower uncertainties
- **Data quality assessment**: Includes goodness-of-fit measures

# When to Use EVM vs. Traditional Methods

**Use EVM when:**
- Data show significant scatter or contain potential outliers
- You need realistic uncertainties that reflect actual data spread
- Working with asymmetric uncertainties (e.g., from likelihood fits, Bayesian inference, bootstrapping)

**Use traditional weighted averaging when:**
- All data points are highly consistent
- Outliers have been carefully removed beforehand
- You need simple, fast computation for large datasets
- Symmetric uncertainties only

# References
- Birch, M., Singh, B., 2014. Method of Best Representation for Averages in Data Evaluation.
  Nuclear Data Sheets, 120, 106-108. [DOI: 10.1016/j.nds.2014.07.019](https://doi.org/10.1016/j.nds.2014.07.019)

"""
module EVMavg

export avg, evm, evm_asymmetric, AvgResult, EVMResult, EVMAsymmetricResult

using Statistics
using SpecialFunctions: erf

"""
    AvgResult

Structure containing results from the `avg` function.

# Fields
- `unweighted::Float64`: Unweighted average
- `weighted::Float64`: Weighted average (inverse variance weighting)
- `std_error::Float64`: Standard error of the weighted average
"""
struct AvgResult
    unweighted::Float64
    weighted::Float64
    std_error::Float64
end

"""
    EVMResult

Structure containing results from the `evm` function.

# Fields
- `mean::Float64`: EVM average value
- `internal_error::Float64`: Internal uncertainty (from data precision)
- `external_error::Float64`: External uncertainty (from data scatter)
- `weights::Vector{Float64}`: Normalized weight for each data point
- `confidence_level::Float64`: Approximate confidence level (goodness of fit)
"""
struct EVMResult
    mean::Float64
    internal_error::Float64
    external_error::Float64
    weights::Vector{Float64}
    confidence_level::Float64
end

"""
    EVMAsymmetricResult

Structure containing results from the `evm_asymmetric` function.

# Fields
- `mean::Float64`: EVM average value
- `lower_error::Float64`: Lower (negative direction) uncertainty
- `upper_error::Float64`: Upper (positive direction) uncertainty
- `internal_lower::Float64`: Internal lower uncertainty
- `internal_upper::Float64`: Internal upper uncertainty
- `external_error::Float64`: External uncertainty (symmetric)
- `weights::Vector{Float64}`: Normalized weights for each measurement
- `confidence_level::Float64`: Approximate confidence level (goodness of fit)
"""
struct EVMAsymmetricResult
    mean::Float64
    lower_error::Float64
    upper_error::Float64
    internal_lower::Float64
    internal_upper::Float64
    external_error::Float64
    weights::Vector{Float64}
    confidence_level::Float64
end

function Base.show(io::IO, result::EVMResult)
    println(io, "EVMResult:")
    println(io, "  Mean: $(result.mean)")
    println(io, "  Internal error: $(result.internal_error)")
    println(io, "  External error: $(result.external_error)")
    println(io, "  Confidence level: $(round(result.confidence_level * 100, digits=2))%")
end

function Base.show(io::IO, result::AvgResult)
    println(io, "AvgResult:")
    println(io, "  Unweighted average: $(result.unweighted)")
    println(io, "  Weighted average: $(result.weighted)")
    println(io, "  Standard error: $(result.std_error)")
end

function Base.show(io::IO, result::EVMAsymmetricResult)
    println(io, "EVMAsymmetricResult:")
    println(io, "  Mean: $(result.mean) +$(result.upper_error) -$(result.lower_error)")
    println(io, "  Internal errors: +$(result.internal_upper) -$(result.internal_lower)")
    println(io, "  External error: ±$(result.external_error)")
    println(io, "  Confidence level: $(round(result.confidence_level * 100, digits=2))%")
end

"""
    avg(data::AbstractVector, sigmas::AbstractVector) -> AvgResult

Compute unweighted and weighted averages with standard error.

The weighted average uses inverse variance weighting: w_i = 1/σ_i².

# Arguments
- `data::AbstractVector`: Data values
- `sigmas::AbstractVector`: Uncertainties (standard deviations) for each data point

# Returns
- `AvgResult`: Structure containing unweighted average, weighted average, and standard error

# Throws
- `ArgumentError`: If arrays have different lengths or contain invalid values

# Examples
```julia
data = [10.2, 10.5, 10.1, 10.4]
sigmas = [0.3, 0.2, 0.4, 0.25]
result = avg(data, sigmas)
println(result.weighted)  # Weighted average
println(result.std_error)  # Standard error
```
"""
function avg(data::AbstractVector{<:Real}, sigmas::AbstractVector{<:Real})
    n = length(data)
    
    # Input validation
    if length(sigmas) != n
        throw(ArgumentError("data and sigmas must have the same length"))
    end
    if n == 0
        throw(ArgumentError("data cannot be empty"))
    end
    if any(sigmas .<= 0)
        throw(ArgumentError("all sigmas must be positive"))
    end
    if any(isnan.(data)) || any(isnan.(sigmas))
        throw(ArgumentError("data and sigmas cannot contain NaN values"))
    end
    
    # Handle single measurement
    if n == 1
        return AvgResult(data[1], data[1], sigmas[1])
    end
    
    # Calculate weights
    w = 1 ./ sigmas.^2
    sw = sum(w)
    
    # Unweighted average
    unweighted_avg = mean(data)
    
    # Weighted average
    weighted_avg = sum(data .* w) / sw
    
    # Standard error of weighted average
    variance = sum((data .- unweighted_avg).^2) / (n - 1)
    effective_n = sw^2 / sum(w.^2)
    std_err = sqrt(variance / effective_n)
    
    return AvgResult(unweighted_avg, weighted_avg, std_err)
end

"""
    evm(data::AbstractVector, errors::AbstractVector; compute_confidence=true) -> EVMResult

Compute the Expected Value Method (EVM) average for data with uncertainties.

The EVM constructs a mean probability density from the data and computes
a weighted average that is robust to outliers and provides realistic uncertainties.

# Arguments
- `data::AbstractVector`: Data values
- `errors::AbstractVector`: Uncertainties for each data point
- `compute_confidence::Bool=true`: Whether to compute the confidence level (goodness of fit)

# Returns
- `EVMResult`: Structure containing EVM average, uncertainties, weights, and confidence level

# Details
The method assumes each data point follows a Gaussian distribution with the given
value and uncertainty. The final result is a weighted combination where weights
are proportional to the expected frequency of each data point given the full dataset.

Internal uncertainty reflects the precision of individual data points, while external 
uncertainty reflects the spread of the data. The larger of the two should typically be reported.

# Throws
- `ArgumentError`: If arrays have different lengths or contain invalid values

# Examples
```julia
data = [5.7, 6.7, 4.3, 4.2, 8.29, 6.0, 5.3, 6.16]
sigmas = [0.3, 0.6, 2.0, 1.5, 1.24, 1.3, 0.8, 0.34]
result = evm(data, sigmas)

println("EVM average: \$(result.mean) ± \$(result.external_error)")
println("Confidence: \$(round(result.confidence_level * 100, digits=1))%")
```

# References
Birch, M., Singh, B., 2014, Nuclear Data Sheets, 120, 106-108.
"""
function evm(data::AbstractVector{<:Real}, 
             errors::AbstractVector{<:Real}; 
             compute_confidence::Bool=true)
    
    n = length(data)
    
    # Input validation
    if length(errors) != n
        throw(ArgumentError("data and errors must have the same length"))
    end
    if n == 0
        throw(ArgumentError("data cannot be empty"))
    end
    if any(errors .<= 0)
        throw(ArgumentError("all errors must be positive"))
    end
    if any(isnan.(data)) || any(isnan.(errors))
        throw(ArgumentError("data and errors cannot contain NaN values"))
    end
    
    # Handle single measurement
    if n == 1
        return EVMResult(data[1], errors[1], errors[1], [1.0], 1.0)
    end
    
    # Compute mean probability density at each measured point
    M = zeros(n)
    for i in 1:n
        for j in 1:n
            M[i] += gaussian_pdf(data[i], data[j], errors[j])
        end
        M[i] /= n
    end
    
    # Compute normalized weights
    sum_M = sum(M)
    weights = M ./ sum_M
    
    # EVM average
    evm_mean = sum(weights .* data)
    
    # Internal uncertainty (from measurement precision)
    internal_error = sqrt(sum(weights.^2 .* errors.^2))
    
    # External uncertainty (from measurement spread)
    external_error = sqrt(sum(weights .* (data .- evm_mean).^2))
    
    # Confidence level (goodness of fit)
    confidence = compute_confidence ? calculate_confidence(data, evm_mean, M, n) : NaN
    
    return EVMResult(evm_mean, internal_error, external_error, weights, confidence)
end

"""
    evm_asymmetric(data::AbstractVector, lower_errors::AbstractVector, upper_errors::AbstractVector; 
                   compute_confidence=true) -> EVMAsymmetricResult

Compute the Expected Value Method (EVM) average for data with asymmetric uncertainties.

This function implements the asymmetric uncertainty treatment described in Equation 1 of 
Birch & Singh (2014), using a piecewise Gaussian probability density function that has
different widths above and below the mean.

# Arguments
- `data::AbstractVector`: Data values
- `lower_errors::AbstractVector`: Lower (negative direction) uncertainties
- `upper_errors::AbstractVector`: Upper (positive direction) uncertainties
- `compute_confidence::Bool=true`: Whether to compute the confidence level (goodness of fit)

# Returns
- `EVMAsymmetricResult`: Structure containing EVM average, asymmetric uncertainties, weights, and confidence level

# Details
For each data point μᵢ with asymmetric uncertainty +aᵢ/-bᵢ, the probability density is:

    A(x; μ, a, b) = √(2/π(a+b)²) × exp(-(x-μ)²/(2σ²))

where σ = b for x ≤ μ and σ = a for x > μ.

The internal uncertainties are calculated separately for upper and lower bounds:
- σ_int+ = √(Σᵢ wᵢ² aᵢ²)
- σ_int- = √(Σᵢ wᵢ² bᵢ²)

The external uncertainty is symmetric and reflects the overall spread:
- σ_ext = √(Σᵢ wᵢ(μᵢ - x_EVM)²)

**Final reported uncertainties** (`upper_error` and `lower_error`):
For each direction, the final uncertainty is the square root of the maximum variance:
- upper_error = √(max(σ_int+², σ_ext²))
- lower_error = √(max(σ_int-², σ_ext²))

This conservative approach ensures the uncertainty reflects both data precision
and data scatter. The internal uncertainties tell you about precision limitations,
while external uncertainty tells you about consistency. The final errors incorporate both.

# Examples
```julia
# Data with asymmetric errors (e.g., from likelihood fits)
data = [10.2, 10.5, 10.1]
lower_errors = [0.3, 0.4, 0.2]  # -σ
upper_errors = [0.5, 0.6, 0.4]  # +σ

result = evm_asymmetric(data, lower_errors, upper_errors)
println("Result: \$(result.mean) +\$(result.upper_error) -\$(result.lower_error)")
```

# References
Birch, M., Singh, B., 2014, Nuclear Data Sheets, 120, 106-108, Equation 1.
"""
function evm_asymmetric(data::AbstractVector{<:Real}, 
                       lower_errors::AbstractVector{<:Real},
                       upper_errors::AbstractVector{<:Real}; 
                       compute_confidence::Bool=true)
    
    n = length(data)
    
    # Input validation
    if length(lower_errors) != n || length(upper_errors) != n
        throw(ArgumentError("data, lower_errors, and upper_errors must have the same length"))
    end
    if n == 0
        throw(ArgumentError("data cannot be empty"))
    end
    if any(lower_errors .<= 0) || any(upper_errors .<= 0)
        throw(ArgumentError("all errors must be positive"))
    end
    if any(isnan.(data)) || any(isnan.(lower_errors)) || any(isnan.(upper_errors))
        throw(ArgumentError("data and errors cannot contain NaN values"))
    end
    
    # Handle single measurement
    if n == 1
        return EVMAsymmetricResult(data[1], lower_errors[1], upper_errors[1], 
                                  lower_errors[1], upper_errors[1], 0.0, [1.0], 1.0)
    end
    
    # Compute mean probability density at each measured point
    M = zeros(n)
    for i in 1:n
        for j in 1:n
            M[i] += asymmetric_gaussian_pdf(data[i], data[j], 
                                           lower_errors[j], upper_errors[j])
        end
        M[i] /= n
    end
    
    # Compute normalized weights
    sum_M = sum(M)
    weights = M ./ sum_M
    
    # EVM average
    evm_mean = sum(weights .* data)
    
    # Internal uncertainties (from measurement precision)
    # As defined in Equation 4 of Birch & Singh (2014)
    internal_upper = sqrt(sum(weights.^2 .* upper_errors.^2))
    internal_lower = sqrt(sum(weights.^2 .* lower_errors.^2))
    
    # External uncertainty (from measurement spread)
    # Symmetric, as defined in Equation 5
    external_error = sqrt(sum(weights .* (data .- evm_mean).^2))
    
    # Determine final uncertainties (larger variance wins)
    # Convert to variances, compare, then back to standard deviations
    upper_var_int = internal_upper^2
    lower_var_int = internal_lower^2
    ext_var = external_error^2
    
    final_upper = sqrt(max(upper_var_int, ext_var))
    final_lower = sqrt(max(lower_var_int, ext_var))
    
    # Confidence level (goodness of fit)
    confidence = compute_confidence ? calculate_confidence(data, evm_mean, M, n) : NaN
    
    return EVMAsymmetricResult(evm_mean, final_lower, final_upper,
                              internal_lower, internal_upper, external_error,
                              weights, confidence)
end

"""
    gaussian_pdf(x, μ, σ)

Compute the Gaussian probability density function at x.

# Arguments
- `x`: Point at which to evaluate the PDF
- `μ`: Mean of the distribution
- `σ`: Standard deviation

# Returns
- Probability density at x
"""
function gaussian_pdf(x::Real, μ::Real, σ::Real)
    return 1.0 / (sqrt(2π) * σ) * exp(-(x - μ)^2 / (2 * σ^2))
end

"""
    asymmetric_gaussian_pdf(x, μ, lower_error, upper_error)

Compute the asymmetric Gaussian probability density function at x.

Implements Equation 1 from Birch & Singh (2014):
- For x ≤ μ: uses lower_error (b) as the width
- For x > μ: uses upper_error (a) as the width
- Normalization factor ensures the PDF integrates to 1

# Arguments
- `x`: Point at which to evaluate the PDF
- `μ`: Mean of the distribution
- `lower_error`: Width below the mean (b in the paper)
- `upper_error`: Width above the mean (a in the paper)

# Returns
- Probability density at x
"""
function asymmetric_gaussian_pdf(x::Real, μ::Real, 
                                lower_error::Real, upper_error::Real)
    # Normalization factor from Equation 1
    norm = sqrt(2 / (π * (upper_error + lower_error)^2))
    
    if x <= μ
        # Use lower error (b) for x below mean
        return norm * exp(-(x - μ)^2 / (2 * lower_error^2))
    else
        # Use upper error (a) for x above mean
        return norm * exp(-(x - μ)^2 / (2 * upper_error^2))
    end
end

"""
    calculate_confidence(data, mean, M, n)

Calculate the approximate confidence level for the EVM result.

Uses a modified Chi-Square test based on the number of data points
above and below the mean compared to expected probabilities.

# Arguments
- `data`: Data values
- `mean`: EVM average
- `M`: Mean probability density at each point
- `n`: Number of data points

# Returns
- Confidence level (0 to 1)
"""
function calculate_confidence(data::AbstractVector, mean::Real, 
                             M::AbstractVector, n::Int)
    # Count measurements below and above mean
    n_low = sum(data .< mean)
    n_high = n - n_low
    
    # Expected probabilities (approximate integration)
    sum_M = sum(M)
    p_low = sum(M[data .< mean]) / sum_M
    p_high = 1.0 - p_low
    
    # Avoid division by zero
    if n * p_low < 1e-10 || n * p_high < 1e-10
        return NaN
    end
    
    # Modified Chi-Square statistic
    Q = (n_low - n * p_low)^2 / (n * p_low) + 
        (n_high - n * p_high)^2 / (n * p_high)
    
    # Convert to confidence level using error function
    confidence = 1.0 - erf(sqrt(Q / 2))
    
    return clamp(confidence, 0.0, 1.0)
end

end # module
