##############################################################
## Test error propagation for the different error functions ##
##############################################################
@testset "error propagation" begin

    ε = 100*eps(Float64)
    @polyvar x y z

    ###############
    # Translation #
    ###############
    # Coefficients
    Cspher = SphericalHarmonicCoefficients(zeros(25)); 
    Cspher[3,-2] = 1/sqrt(15)*sqrt((4*pi)/(2*3+1));
    Csolid = SphericalHarmonicCoefficients(zeros(25),1.0,true); 
    Csolid[3,-2] = 1/sqrt(15);

    # Translation
    v = [1,-1,2]

    coeffsTransSpher = SphericalHarmonicExpansions.errorTranslation(Cspher,v)
    coeffsTransSolid = SphericalHarmonicExpansions.errorTranslation(Csolid,v)

    # Correct solution
    CtransSpher = zeros(25)
    CtransSpher[1:11] = sqrt(4pi/3).*[2*sqrt(3), 2, 1, 2, 2/sqrt(5), 1/sqrt(5), 0, 1/sqrt(5), 0, 0, 1/sqrt(35)]
    CtransSpher = SphericalHarmonicCoefficients(CtransSpher)
    CtransSolid = zeros(25)
    CtransSolid[1:11] = [2, 2, 1, 2, 2/sqrt(3), 1/sqrt(3), 0, 1/sqrt(3), 0, 0, 1/sqrt(15)]
    CtransSolid = SphericalHarmonicCoefficients(CtransSolid,1.0,true)

    # Test
    @test isapprox(CtransSpher,coeffsTransSpher)
    @test isapprox(CtransSolid,coeffsTransSolid)

    ##############
    # Quadrature #
    ##############
    # 3-design:
    coordinates = [1. 0. 0.;
                   -1. 0. 0.;
                   0. 1. 0.;
                   0. -1. 0.;
                   0. 0. 1.;
                   0. 0. -1.]

    # Calculating the error of the coefficients up to l=1
    δ = ones(6)
    εquadr = SphericalHarmonicExpansions.errorSphericalQuadrature(δ,coordinates,1);

    # Correct solution
    εcorr = [sqrt(2*pi/3); sqrt(2*pi/3); sqrt(2*pi/3); sqrt(2*pi/3)]

    # Test
    @test isapprox(εquadr.c,εcorr,atol=ε)
end
