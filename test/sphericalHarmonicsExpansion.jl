@testset "spherical harmonics expansion" begin
    ɛ = eps(Float64)
    @polyvar x y z
    @test_throws DomainError SphericalHarmonicCoefficients(-2)
    @test_throws DomainError SphericalHarmonicCoefficients([1,2,3,4,5])
    @test_throws DomainError SphericalHarmonicCoefficients(-2,1.0,true)
    @test_throws DomainError SphericalHarmonicCoefficients([1,2,3,4,5],1.0,true)
    @test_throws DomainError SphericalHarmonicCoefficients([[1,2,3,4],[1,2,3,4]],[1.0],BitArray(undef,(1, 0)))
    @test_throws DomainError SphericalHarmonicCoefficients([[1,2,3,4],[1,2,3,4]],[1.0, 1.0],BitArray(undef,(1,)))
    @test_throws DomainError SphericalHarmonicCoefficients([[1,2,3,4],[1,2,3]],[1.0, 1.0],BitArray(undef,(1, 0)))
    @test_throws DomainError SphericalHarmonicCoefficients([[1,2,3,4],[1,2]],1.0,true)
    c1 = SphericalHarmonicCoefficients(2)
    c2 = SphericalHarmonicCoefficients([1,2,3,4])
    c3 = SphericalHarmonicCoefficients([1+ɛ,2,3,4])
    c4 = SphericalHarmonicCoefficients([3,2,3,4])
    c5 = SphericalHarmonicCoefficients(2,1.0,true)
    c6 = SphericalHarmonicCoefficients([1,2,3,4],2.0,false)
    c7 = SphericalHarmonicCoefficients([[2],[1,2,3,4]],1.0,false)
    @test c2==c2
    @test c2!=c3
    @test isapprox(c2,c3)
    @test !isapprox(c2,c4)
    @test !isapprox(c1,c5)
    @test !isapprox(c2,c6)
    @test isapprox(c2,c7[2])
    @test length(c1.c) == 9
    @test c2.L == 1
    @test c7[1].L != c7[2].L
    @test c2[1] == 1
    @test c2[2] == 2
    @test c2[3] == 3
    @test c2[4] == 4
    @test c2[0,0] == 1
    @test c2[1,-1] == 2
    @test c2[1,0] == 3
    @test c2[1,1] == 4
    c2[1] = 10
    c2[2] = 20
    c2[3] = 30
    c2[4] = 40
    @test c2[1] == 10
    @test c2[2] == 20
    @test c2[3] == 30
    @test c2[4] == 40
    c2[0,0] = 100
    c2[1,-1] = 200
    c2[1,0] = 300
    c2[1,1] = 400
    @test c2[1] == 100
    @test c2[2] == 200
    @test c2[3] == 300
    @test c2[4] == 400
    #Test by using unit vectors
    #l = 10 -> size(C) = l²+2l+1 = 121
    for l in 0:2
        for m in -l:l
            C = SphericalHarmonicCoefficients(2)
            C[l,m] = 1
            @test isapprox(sphericalHarmonicsExpansion(C,x,y,z),sqrt((2*l+1)/(4*pi))*SphericalHarmonicExpansions.zlm(l,m,x,y,z),atol=ɛ)
        end
    end
    #Testing the addition
    C = SphericalHarmonicCoefficients([1;0.5;10.78;0.87456])
    #@test isapprox(sphericalHarmonicsExpansion(C,x,y,z),1*ylm(0,0,x,y,z)+0.5*ylm(1,-1,x,y,z)+10.78*ylm(1,0,x,y,z)+0.87456*ylm(1,1,x,y,z),atol=ɛ)
    # Test: solid & spherical harmonic expansion
    c = ones(4)
    Csolid = SphericalHarmonicCoefficients(c,1.0,true)
    Cspher = SphericalHarmonicCoefficients(c,1.0,false)
    @test isapprox(sphericalHarmonicsExpansion(Csolid,x,y,z),1+x+y+z,atol=ɛ)
    @test isapprox(sphericalHarmonicsExpansion(Cspher,x,y,z),sqrt(1/(4pi))+sqrt(3/(4pi))*x+sqrt(3/(4pi))*y+sqrt(3/(4pi))*z,atol=ɛ)

    # Type-generic complete expansion: low-order Float64 result must agree with
    # the legacy Float64 baseline. BigFloat construction must also be available.
    Ctyped = SphericalHarmonicCoefficients([1.0, 0.5, -0.25, 0.75])
    @test isapprox(
        sphericalHarmonicsExpansion_typed(Float64, Ctyped, x, y, z; scaling=:distributed),
        sphericalHarmonicsExpansion(Ctyped, x, y, z),
        atol=16eps(Float64),
    )
    @test_nowarn sphericalHarmonicsExpansion_bigfloat(Ctyped, x, y, z; precision=128, scaling=:adaptive)
end

@testset "coefficients" begin
    ɛ = eps(Float64)
    # Test: solid ↔ spherical
    csolid = zeros(25)
    csolid[1:4] = [sqrt(1/(4pi)),sqrt(3/(4pi)),sqrt(3/(4pi)),sqrt(3/(4pi))]
    cspher = zeros(25)
    cspher[1:4] = [1,1,1,1]
    Csolid = SphericalHarmonicCoefficients(csolid,0.042,true)
    Cspher = SphericalHarmonicCoefficients(cspher,0.042,false)
    @test isapprox(spherical!(deepcopy(Csolid)),Cspher,atol=ɛ)
    @test isapprox(solid!(deepcopy(Cspher)),Csolid,atol=ɛ)
    @test isapprox(spherical!(deepcopy(Cspher)),Cspher,atol=ɛ)
    @test isapprox(solid!(Csolid),Csolid,atol=ɛ)
    # Test: Normalization
    c1 = ones(9)
    C1 = SphericalHarmonicCoefficients(c1,1.0,true)
    c2 = [1,2,2,2,4,4,4,4,4]
    C2 = SphericalHarmonicCoefficients(c2,2.0,true)
    @test isapprox(normalize!(C1,1/2.0),C2,atol=ɛ)
    normalize!(C1,2.0)
    @test isapprox(normalize(C2,C2.R),C1,atol=ɛ)
    # Test: generate an array of SphericalHarmonicCoefficients
    c3 = reshape([ones(9) for i=1:6],2,3)
    R = ones(Float64,2,3)
    sol = BitArray(R)
    C3 = SphericalHarmonicCoefficients(c3,R,sol)
    for (i,co) in enumerate(C3)
        @test isapprox(C1,co,atol=ɛ)
    end
    # Test: Operations
    c4 = 2*ones(9)
    C4 = SphericalHarmonicCoefficients(c4,1.0,true)
    @test isapprox(C1*2,C4,atol=ɛ)
    @test isapprox(C4/2,C1,atol=ɛ)
    @test isapprox(1.0+C1,C4,atol=ɛ)
    @test isapprox(3.0-C1,C4,atol=ɛ)
    @test isapprox(C4-1.0,C1,atol=ɛ)
    C5 = SphericalHarmonicCoefficients(zeros(9),2.0,true)
    @test isapprox(2*C1,C1+C1,atol=ɛ)
    @test isapprox(C2-C2,C5,atol=ɛ)
    @test_throws DomainError Csolid+Cspher
    @test_throws DomainError Csolid-Cspher
    @test_throws DomainError C4+C5
    @test_throws DomainError C4-C5
    # Test: Save and read coefficients to/from an HDF5 file
    write("Test1.h5",[C2])
    CTest1 = SphericalHarmonicCoefficients("Test1.h5")
    @test isapprox(CTest1[1],C2,atol=ɛ)
    write("Test2.h5",C3)
    CTest2 = SphericalHarmonicCoefficients("Test2.h5")
    for i=1:length(C3)
        @test isapprox(CTest2[i],C3[i],atol=ɛ)
    end
    rm("Test1.h5")
    rm("Test2.h5")
end