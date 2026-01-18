"""
Basic usage examples for EVMavg.jl

This script demonstrates the main features of the package.
"""

using EVMavg, Statistics

# Example 1: Simple EVM averaging
println("=" ^ 60)
println("Example 1: Basic EVM Averaging")
println("=" ^ 60)

# Data and sigmas
data = [5.7, 6.7, 4.3, 6.0, 5.3, 6.16, 5.48]
sigmas = [0.3, 0.6, 2.0, 1.3, 0.8, 0.34, 0.545]

result = evm(data, sigmas)

println("\nMeasurements:")
for (i, (t, e)) in enumerate(zip(data, sigmas))
    println("  $i: $t ± $e")
end

println("\nEVM Results:")
println(result)

println("\nWeights for each measurement:")
for (i, w) in enumerate(result.weights)
    println("  Measurement $i: $(round(w, digits=4))")
end

println("\nRecommendation:")
if result.external_error > result.internal_error
    println("  Report external error: $(result.mean) ± $(round(result.external_error, digits=2))")
else
    println("  Report internal error: $(result.mean) ± $(round(result.internal_error, digits=2))")
end

# Example 2: Comparison with weighted average
println("\n" * "=" ^ 60)
println("Example 2: EVM vs. Weighted Average")
println("=" ^ 60)

# Data with an outlier
data = [10.1, 10.2, 10.0, 10.3, 10.1, 14.5]  # Last value is outlier
sigmas = [0.2, 0.2, 0.2, 0.2, 0.2, 0.3]

evm_result = evm(data, sigmas)
avg_result = avg(data, sigmas)

println("\nData with outlier:")
for (i, (d, s)) in enumerate(zip(data, sigmas))
    println("  $i: $d ± $s")
end

println("\nWeighted Average:")
println("  Result: $(round(avg_result.weighted, digits=2)) ± $(round(avg_result.std_error, digits=2))")

println("\nEVM Average:")
println("  Result: $(round(evm_result.mean, digits=2)) ± $(round(evm_result.external_error, digits=2))")
println("  Outlier weight: $(round(evm_result.weights[end], digits=4))")
println("  Typical weight: $(round(mean(evm_result.weights[1:end-1]), digits=4))")


println("\n" * "=" ^ 60)
println("Example 3: Asymmetric Uncertainties")
println("=" ^ 60)

# Measurements with asymmetric errors (common in likelihood fits)
# For example, particle lifetime measurements where systematic
# uncertainties are different in each direction
data_asym = [12.5, 13.2, 12.8, 13.0]
lower_errors = [0.8, 1.0, 0.6, 0.7]  # -σ (smaller)
upper_errors = [1.5, 2.0, 1.2, 1.4]  # +σ (larger)

println("\nMeasurements with asymmetric errors:")
for (i, (d, le, ue)) in enumerate(zip(data_asym, lower_errors, upper_errors))
    println("  $i: $d +$ue -$le")
end

result_asym = evm_asymmetric(data_asym, lower_errors, upper_errors)

println("\nAsymmetric EVM Results:")
println(result_asym)

println("\nUnderstanding the uncertainties:")
println("  Internal upper: $(round(result_asym.internal_upper, digits=3)) (from measurement precision)")
println("  Internal lower: $(round(result_asym.internal_lower, digits=3)) (from measurement precision)")
println("  External:       $(round(result_asym.external_error, digits=3)) (from data scatter, symmetric)")
println()
println("  Final upper:    $(round(result_asym.upper_error, digits=3)) = max(internal_upper, external)")
println("  Final lower:    $(round(result_asym.lower_error, digits=3)) = max(internal_lower, external)")
println()
println("  ✓ REPORT: $(round(result_asym.mean, digits=2)) +$(round(result_asym.upper_error, digits=2)) -$(round(result_asym.lower_error, digits=2))")

println("\nComparison with symmetric EVM:")
# Using average of lower and upper as symmetric error
avg_errors = (lower_errors .+ upper_errors) ./ 2
result_sym = evm(data_asym, avg_errors)

println("  Symmetric:  $(round(result_sym.mean, digits=2)) ± $(round(result_sym.external_error, digits=2))")
println("  Asymmetric: $(round(result_asym.mean, digits=2)) +$(round(result_asym.upper_error, digits=2)) -$(round(result_asym.lower_error, digits=2))")
println("\nNote: Asymmetric method preserves the different error magnitudes in each direction")
