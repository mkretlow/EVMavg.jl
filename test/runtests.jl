using EVMavg
using Test

@testset "EVMavg.jl" begin
   
    @testset "avg function - basic functionality" begin
        # Simple test case
        data = [10.0, 10.2, 10.1, 10.3]
        sigmas = [0.2, 0.2, 0.2, 0.2]
        
        result = avg(data, sigmas)
        
        @test result isa AvgResult
        @test result.unweighted ≈ 10.15
        @test result.weighted ≈ 10.15
        @test result.std_error > 0
        
        # Single measurement
        result_single = avg([5.0], [0.5])
        @test result_single.unweighted == 5.0
        @test result_single.weighted == 5.0
        @test result_single.std_error == 0.5

        # Another (legacy) test
        data = [1.0, 2.0, 3.0, 4.0]
        sigmas = [0.1, 0.2, 0.3, 0.4]

        result = avg(data, sigmas)

        @test result.unweighted ≈ 2.5 atol=1E-14
        @test result.weighted ≈ 1.4634146341463417 atol=1E-14
        @test result.std_error ≈ 0.941876146945452 atol=1E-14

    end
    
    @testset "avg function - weighted averaging" begin
        # Different uncertainties should give different weights
        data = [10.0, 12.0]
        sigmas = [0.1, 1.0]  # First measurement much more precise
        
        result = avg(data, sigmas)
        
        # Weighted average should be closer to first value
        @test result.weighted < result.unweighted
        @test result.weighted > 10.0
        @test result.weighted < 10.5
    end
    
    @testset "avg function - input validation" begin
        # Mismatched lengths
        @test_throws ArgumentError avg([1.0, 2.0], [0.1])
        
        # Empty arrays
        @test_throws ArgumentError avg(Float64[], Float64[])
        
        # Non-positive sigmas
        @test_throws ArgumentError avg([1.0, 2.0], [0.1, 0.0])
        @test_throws ArgumentError avg([1.0, 2.0], [0.1, -0.1])
        
        # NaN values
        @test_throws ArgumentError avg([1.0, NaN], [0.1, 0.1])
        @test_throws ArgumentError avg([1.0, 2.0], [0.1, NaN])
    end
    
    @testset "evm function - basic functionality" begin
        # Simple test case with similar values
        data = [10.0, 10.2, 10.1, 10.3, 9.9]
        errors = [0.2, 0.2, 0.2, 0.2, 0.2]
        
        result = evm(data, errors)
        
        @test result isa EVMResult
        @test result.mean ≈ 10.1
        @test result.internal_error > 0
        @test result.external_error > 0
        @test length(result.weights) == length(data)
        @test sum(result.weights) ≈ 1.0
        @test all(result.weights .>= 0)
        @test 0 <= result.confidence_level <= 1
        
        # Single measurement
        result_single = evm([5.0], [0.5])
        @test result_single.mean == 5.0
        @test result_single.internal_error == 0.5
        @test result_single.external_error == 0.5
        @test result_single.weights == [1.0]

        # Another (legacy) test
        data = [1.0, 2.0, 3.0, 4.0]
        sigmas = [0.1, 0.2, 0.3, 0.4]

        result = evm(data, sigmas)

        @test result.mean ≈ 1.9269849828944758 atol=1E-14
        @test result.internal_error ≈ 0.096357574983269 atol=1E-14
        @test result.external_error ≈ 1.0559929440171045 atol=1E-14

    end
    
    @testset "evm function - outlier robustness" begin
        # Data with one clear outlier
        data = [10.0, 10.2, 10.1, 10.3, 15.0]  # Last value is outlier
        errors = [0.2, 0.2, 0.2, 0.2, 0.3]
        
        result = evm(data, errors)
        
        # Outlier should have lowest weight
        @test argmin(result.weights) == 5
        @test result.weights[5] < minimum(result.weights[1:4])
        
        # Mean should be pulled less toward outlier than simple average
        simple_avg = sum(data) / length(data)
        @test abs(result.mean - 10.15) < abs(simple_avg - 10.15)
    end
    
    @testset "evm function - uncertainty behavior" begin
        # Consistent data should have small external error
        data_consistent = [10.0, 10.1, 10.0, 10.1]
        errors_consistent = [0.5, 0.5, 0.5, 0.5]
        
        result_consistent = evm(data_consistent, errors_consistent)
        
        # Internal error should be related to measurement errors
        @test result_consistent.internal_error < maximum(errors_consistent)
        
        # Scattered data should have larger external error
        data_scattered = [8.0, 10.0, 12.0, 14.0]
        errors_scattered = [0.5, 0.5, 0.5, 0.5]
        
        result_scattered = evm(data_scattered, errors_scattered)
        
        # External error should reflect the scatter
        @test result_scattered.external_error > result_consistent.external_error
    end
    
    @testset "evm function - confidence level" begin
        # Very consistent data should have high confidence
        data = [10.0, 10.1, 10.0, 10.1, 10.0, 10.1]
        errors = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        
        result = evm(data, errors)
        @test result.confidence_level > 0.95  # Should be reasonably high
        
        # Can disable confidence calculation
        result_no_conf = evm(data, errors; compute_confidence=false)
        @test isnan(result_no_conf.confidence_level)
    end
    
    @testset "evm function - input validation" begin
        # Mismatched lengths
        @test_throws ArgumentError evm([1.0, 2.0], [0.1])
        
        # Empty arrays
        @test_throws ArgumentError evm(Float64[], Float64[])
        
        # Non-positive errors
        @test_throws ArgumentError evm([1.0, 2.0], [0.1, 0.0])
        @test_throws ArgumentError evm([1.0, 2.0], [0.1, -0.1])
        
        # NaN values
        @test_throws ArgumentError evm([1.0, NaN], [0.1, 0.1])
        @test_throws ArgumentError evm([1.0, 2.0], [0.1, NaN])
    end
    
    @testset "evm function - literature example" begin
        # Simplified version of 48Ti example from Birch & Singh (2014)
        # Using subset of data from Table I
        lifetimes = [5.7, 6.7, 4.3, 6.0, 5.3, 6.16]
        errors = [0.3, 0.6, 2.0, 1.3, 0.8, 0.34]
        
        result = evm(lifetimes, errors)
        
        # Result should be in reasonable range
        @test result.mean ≈ 5.881634209952061 atol=1E-14
        @test result.internal_error > 0
        @test result.external_error > 0
        @test sum(result.weights) ≈ 1.0
        
        # High-error measurement should have lower weight
        max_error_idx = argmax(errors)
        @test result.weights[max_error_idx] < maximum(result.weights)
    end
    
    @testset "gaussian_pdf function" begin
        # Test Gaussian PDF calculation
        # PDF at mean should be maximum
        μ, σ = 0.0, 1.0
        @test EVMavg.gaussian_pdf(μ, μ, σ) ≈ 1/(sqrt(2π)*σ)
        
        # PDF should be symmetric around mean
        @test EVMavg.gaussian_pdf(μ+1, μ, σ) ≈ EVMavg.gaussian_pdf(μ-1, μ, σ)
        
        # PDF should decrease away from mean
        @test EVMavg.gaussian_pdf(μ, μ, σ) > EVMavg.gaussian_pdf(μ+σ, μ, σ)
        @test EVMavg.gaussian_pdf(μ+σ, μ, σ) > EVMavg.gaussian_pdf(μ+2σ, μ, σ)
    end
    
    @testset "Result types - show methods" begin
        # Test that show methods work without errors
        data = [10.0, 10.2, 10.1]
        errors = [0.2, 0.2, 0.2]
        
        evm_result = evm(data, errors)
        avg_result = avg(data, errors)
        
        # Should not throw errors
        io = IOBuffer()
        show(io, evm_result)
        @test !isempty(String(take!(io)))
        
        show(io, avg_result)
        @test !isempty(String(take!(io)))
    end
    
    @testset "Comparison: EVM vs weighted average" begin
        # EVM should handle outliers better
        data = [10.0, 10.1, 10.2, 10.0, 10.1, 15.0]  # Clear outlier
        errors = [0.2, 0.2, 0.2, 0.2, 0.2, 0.3]
        
        evm_result = evm(data, errors)
        avg_result = avg(data, errors)
        
        # EVM should be more robust (closer to true value ~10.1)
        @test abs(evm_result.mean - 10.1) < abs(avg_result.weighted - 10.1)
        
        # Both methods should reflect the scatter, but EVM is designed to be robust
        # External error should be positive and reflect data spread
        @test evm_result.external_error > 0.5
        @test avg_result.std_error > 0
    end
    
    @testset "Edge cases" begin
        # Two measurements
        result = evm([10.0, 12.0], [0.5, 0.5])
        @test 10.0 < result.mean < 12.0
        @test sum(result.weights) ≈ 1.0
        
        # Very different uncertainties with different values
        # Precise measurement should get higher weight
        result = evm([10.0, 11.0], [0.01, 10.0])
        @test result.mean ≈ 10.0 atol=0.5  # Should be close to precise measurement
        @test result.weights[1] > result.weights[2]  # Precise measurement dominates
        
        # All identical values
        result = evm([10.0, 10.0, 10.0], [0.5, 0.5, 0.5])
        @test result.mean ≈ 10.0
        @test result.external_error ≈ 0.0 atol=1e-10
    end
    
    @testset "evm_asymmetric function - basic functionality" begin
        # Simple test case with asymmetric errors
        data = [10.0, 10.2, 10.1, 9.9]
        lower_errors = [0.2, 0.3, 0.2, 0.2]
        upper_errors = [0.3, 0.4, 0.3, 0.3]
        
        result = evm_asymmetric(data, lower_errors, upper_errors)
        
        @test result isa EVMAsymmetricResult
        @test result.mean ≈ 10.051786210771155 atol=1E-14
        @test result.lower_error > 0
        @test result.upper_error > 0
        @test result.internal_lower > 0
        @test result.internal_upper > 0
        @test result.external_error >= 0
        @test length(result.weights) == length(data)
        @test sum(result.weights) ≈ 1.0
        @test all(result.weights .>= 0)
        @test 0 <= result.confidence_level <= 1
        
        # Single measurement
        result_single = evm_asymmetric([5.0], [0.3], [0.5])
        @test result_single.mean == 5.0
        @test result_single.lower_error == 0.3
        @test result_single.upper_error == 0.5
        @test result_single.weights == [1.0]
    end
    
    @testset "evm_asymmetric function - asymmetry preservation" begin
        # Data with strongly asymmetric errors
        data = [10.0, 10.1, 10.2]
        lower_errors = [0.1, 0.1, 0.1]  # Small lower errors
        upper_errors = [1.0, 1.0, 1.0]  # Large upper errors
        
        result = evm_asymmetric(data, lower_errors, upper_errors)
        
        # Upper error should be larger than lower error
        @test result.upper_error > result.lower_error
        @test result.internal_upper > result.internal_lower
    end
    
    @testset "evm_asymmetric function - comparison with symmetric" begin
        # When errors are symmetric, should match regular evm
        data = [10.0, 10.2, 10.1, 10.3]
        errors = [0.2, 0.2, 0.2, 0.2]
        
        result_sym = evm(data, errors)
        result_asym = evm_asymmetric(data, errors, errors)
        
        # Results should be similar when errors are symmetric
        @test result_sym.mean ≈ result_asym.mean atol=1e-10
        @test result_sym.weights ≈ result_asym.weights atol=1e-10
        @test result_asym.lower_error ≈ result_asym.upper_error atol=1e-6
    end
    
    @testset "evm_asymmetric function - input validation" begin
        # Mismatched lengths
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [0.1], [0.2, 0.2])
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [0.1, 0.1], [0.2])
        
        # Empty arrays
        @test_throws ArgumentError evm_asymmetric(Float64[], Float64[], Float64[])
        
        # Non-positive errors
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [0.0, 0.1], [0.2, 0.2])
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [0.1, 0.1], [-0.1, 0.2])
        
        # NaN values
        @test_throws ArgumentError evm_asymmetric([1.0, NaN], [0.1, 0.1], [0.2, 0.2])
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [NaN, 0.1], [0.2, 0.2])
        @test_throws ArgumentError evm_asymmetric([1.0, 2.0], [0.1, 0.1], [0.2, NaN])
    end
    
    @testset "asymmetric_gaussian_pdf function" begin
        # Test asymmetric PDF calculation
        μ = 0.0
        lower_err = 0.5  # b in the paper (narrow side)
        upper_err = 1.0  # a in the paper (wide side)
        
        # PDF should be continuous at μ
        pdf_at_mean = EVMavg.asymmetric_gaussian_pdf(μ, μ, lower_err, upper_err)
        pdf_just_below = EVMavg.asymmetric_gaussian_pdf(μ - 1e-10, μ, lower_err, upper_err)
        pdf_just_above = EVMavg.asymmetric_gaussian_pdf(μ + 1e-10, μ, lower_err, upper_err)
        @test pdf_at_mean ≈ pdf_just_below rtol=1e-6
        @test pdf_at_mean ≈ pdf_just_above rtol=1e-6
        
        # PDF should be asymmetric (different values at ±1 from mean)
        pdf_minus_one = EVMavg.asymmetric_gaussian_pdf(μ - 1.0, μ, lower_err, upper_err)
        pdf_plus_one = EVMavg.asymmetric_gaussian_pdf(μ + 1.0, μ, lower_err, upper_err)
        @test pdf_minus_one != pdf_plus_one
        
        # At same distance from mean, wider side has higher PDF 
        # (narrower side decays faster)
        @test pdf_plus_one > pdf_minus_one
        
        # PDF should decrease away from mean on both sides
        @test EVMavg.asymmetric_gaussian_pdf(μ, μ, lower_err, upper_err) > 
              EVMavg.asymmetric_gaussian_pdf(μ - 0.5, μ, lower_err, upper_err)
        @test EVMavg.asymmetric_gaussian_pdf(μ, μ, lower_err, upper_err) > 
              EVMavg.asymmetric_gaussian_pdf(μ + 1.0, μ, lower_err, upper_err)
    end
end
