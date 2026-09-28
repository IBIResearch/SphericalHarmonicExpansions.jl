@testset "type-generic construction and normalization" begin
    @polyvar tx ty tz

    legacy = ylm(3, 2, tx, ty, tz)
    typed = ylm_typed(Float64, 3, 2, tx, ty, tz; scaling=:distributed)
    @test typed == legacy

    setprecision(BigFloat, 128) do
        ybig = ylm_typed(BigFloat, 3, 2, tx, ty, tz; scaling=:distributed)
        rbig = rlylm_typed(BigFloat, 3, 2, tx, ty, tz; scaling=:distributed)
        zbig = zlm_typed(BigFloat, 3, 2, tx, ty, tz; scaling=:distributed)
        @test all(coefficient -> coefficient isa BigFloat,
                  SphericalHarmonicExpansions.coefficients(ybig))
        @test all(coefficient -> coefficient isa BigFloat,
                  SphericalHarmonicExpansions.coefficients(rbig))
        @test all(coefficient -> coefficient isa BigFloat,
                  SphericalHarmonicExpansions.coefficients(zbig))

        y00 = ylm_typed(BigFloat, 0, 0, tx, ty, tz; scaling=:distributed)
        expected = inv(sqrt(BigFloat(4) * BigFloat(pi)))
        @test isapprox(first(SphericalHarmonicExpansions.coefficients(y00)), expected;
                       rtol=16eps(BigFloat))
    end

    exact_binomial = binomial(big(100), big(50))
    @test binomial_as(BigFloat, 100, 50) == BigFloat(exact_binomial)
    @test_nowarn SphericalHarmonicExpansions.ylmCosSinPolynomial(67, tx, ty)
end

@testset "measured precision policy" begin
    float_recommendation = recommend_precision(l=24, tolerance=1e-8)
    @test float_recommendation.status == :supported
    @test float_recommendation.numeric_type === Float64
    @test float_recommendation.precision_bits == 53
    @test float_recommendation.measured_bound <= BigFloat("1e-8")

    big_recommendation = recommend_precision(l=24, tolerance=1e-20)
    @test big_recommendation.status == :supported
    @test big_recommendation.numeric_type === BigFloat
    @test big_recommendation.precision_bits == 128
    @test big_recommendation.measured_bound <= BigFloat("1e-20")

    @test recommend_precision(l=65, tolerance=1e-8).status == :order_not_measured
    @test recommend_precision(l=24, tolerance=5e-9).status == :tolerance_not_measured
    @test recommend_precision(l=24, tolerance=1e-31).status == :not_supported_by_measured_grid
    @test recommend_precision(l=24, tolerance=1e-8, basis=:rlylm).status ==
          :configuration_not_measured
    @test recommend_precision(l=24, tolerance=1e-8, variant=:adaptive).status ==
          :configuration_not_measured
    @test recommend_precision(l=24, tolerance=1e-8, mode=:nearest).status ==
          :unsupported_mode

    float_evaluation = evaluate_with_precision(
        :ylm, 2, 0, 0.3, 0.4, 0.5; tolerance=1e-8,
    )
    @test float_evaluation.status == :supported
    @test float_evaluation.value isa Float64

    big_evaluation = evaluate_with_precision(
        :ylm, 2, 0, 3//10, 4//10, 5//10; tolerance=1e-20,
    )
    @test big_evaluation.status == :supported
    @test big_evaluation.value isa BigFloat
    @test precision(big_evaluation.value) == big_evaluation.precision_bits

    unsupported_evaluation = evaluate_with_precision(
        :ylm, 65, 0, 0.0, 0.0, 1.0; tolerance=1e-8,
    )
    @test unsupported_evaluation.status == :order_not_measured
    @test unsupported_evaluation.value === nothing
end
