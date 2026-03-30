@testset "spherical harmonics" begin
	ɛ = 4*eps(Float64)
    @polyvar x y z
    
    # test OverflowError for large m
    @test_throws OverflowError SphericalHarmonicExpansions.ylmCosSinPolynomial(67, x,y)
    @test_throws OverflowError SphericalHarmonicExpansions.ylmSinSinPolynomial(67, x,y)

	# test expasion of (x²+y²+z²)^n
	for n = 0:10
		@test SphericalHarmonicExpansions.trinomialExpansion(n , x, y, z) == (x^2+y^2+z^2)^n
	end

	#@test_throws(ErrorException(BoundsError), ylm(0,1,x,y,z))
	@test_throws DomainError ylm(0,1,x,y,z)
	# l = 0
	#Y_{0,0} = √{\frac{1}{4π}}
	@test isapprox(ylm(0,0,x,y,z),sqrt(1/(4pi))+0*x,atol=ɛ)
	# l = 1
	#Y_{1,-1} = √{\frac{3}{4π}}\frac{y}{r}
	#Y_{1,0} = √{\frac{3}{4π}}\frac{z}{r}
	#Y_{1,1} = √{\frac{3}{4π}}\frac{x}{r}
	@test isapprox(ylm(1,-1,x,y,z),0+sqrt(3/(4pi))*y,atol=ɛ)
	@test isapprox(ylm(1,0,x,y,z),0+sqrt(3/(4pi))*z,atol=ɛ)
	@test isapprox(ylm(1,1,x,y,z),0+sqrt(3/(4pi))*x,atol=ɛ)
	# l = 2
	#Y_{2,-2} = √{\frac{15}{4π}}\frac{xy}{r²}
	#Y_{2,-1} = √{\frac{15}{4π}}\frac{yz}{r²}
	#Y_{2,0} = √{\frac{5}{16π}}\frac{3z² - 1}{r²}
	#Y_{2,1} = √{\frac{15}{4π}}\frac{zx}{r²}
	#Y_{2,2} = √{\frac{15}{16π}}\frac{x²-y²}{r²}
	@test isapprox(ylm(2,-2,x,y,z),0+sqrt(15/(4pi))*x*y,atol=ɛ)
	@test isapprox(ylm(2,-1,x,y,z),0+sqrt(15/(4pi))*y*z,atol=ɛ)
	@test isapprox(ylm(2,0,x,y,z),sqrt(5/(16pi))*(3*z^2 - 1),atol=ɛ)
	@test isapprox(ylm(2,1,x,y,z),0+sqrt(15/(4pi))*z*x,atol=ɛ)
	@test isapprox(ylm(2,2,x,y,z),sqrt(15/(16pi))*(x^2-y^2),atol=ɛ)
	# l = 3
	#Y_{3,-3} = √{\frac{35}{32π}}\frac{(3x²-y²)y}{r³}
	#Y_{3,-2} = √{\frac{105}{4π}}\frac{xyz}{r³}
	#Y_{3,-1} = √{\frac{21}{32π}}\frac{y(5z²-1)}{r³}
	#Y_{3,0} = √{\frac{7}{16π}}\frac{5z³-3z}{r³}
	#Y_{3,1} = √{\frac{21}{32π}}\frac{x(5z²-1)}{r³}
	#Y_{3,2} = √{\frac{105}{16π}}\frac{(x²-y²)z}{r³}
	#Y_{3,3} = √{\frac{35}{32π}}\frac{(x²-3y²)x}{r³}
	@test isapprox(ylm(3,-3,x,y,z),sqrt(35/(32pi))*(3*x^2-y^2)*y,atol=ɛ)
	@test isapprox(ylm(3,-2,x,y,z),0+sqrt(105/(4pi))*x*y*z,atol=ɛ)
	@test isapprox(ylm(3,-1,x,y,z),sqrt(21/(32pi))*y*(5*z^2-1),atol=ɛ)
	@test isapprox(ylm(3,0,x,y,z),sqrt(7/(16pi))*(5*z^3-3*z),atol=ɛ)
	@test isapprox(ylm(3,1,x,y,z),sqrt(21/(32pi))*x*(5*z^2-1),atol=ɛ)
	@test isapprox(ylm(3,2,x,y,z),sqrt(105/(16pi))*(x^2-y^2)*z,atol=ɛ)
	@test isapprox(ylm(3,3,x,y,z),sqrt(35/(32pi))*(x^2-3*y^2)*x,atol=ɛ)
end

@testset "evaluate r^l*ylm" begin
	ε = 32*eps(Float64)
	@polyvar x y z
	
	@test isapprox(rlylm(0,0,x,y,z),sqrt(1/(4pi))+0*x,atol=ɛ)
	
	@test isapprox(rlylm(1,-1,x,y,z),0+sqrt(3/(4pi))*y,atol=ɛ)
	@test isapprox(rlylm(1,0,x,y,z),0+sqrt(3/(4pi))*z,atol=ɛ)
	@test isapprox(rlylm(1,1,x,y,z),0+sqrt(3/(4pi))*x,atol=ɛ)
	
	@test isapprox(rlylm(2,-2,x,y,z),0+sqrt(15/(4pi))*x*y,atol=ɛ)
	@test isapprox(rlylm(2,-1,x,y,z),0+sqrt(15/(4pi))*y*z,atol=ɛ)
	@test isapprox(rlylm(2,0,x,y,z),sqrt(5/(16pi))*(2*z^2 - x^2 - y^2),atol=ɛ)
	@test isapprox(rlylm(2,1,x,y,z),0+sqrt(15/(4pi))*z*x,atol=ɛ)
	@test isapprox(rlylm(2,2,x,y,z),sqrt(15/(16pi))*(x^2-y^2),atol=ɛ)

	@test isapprox(rlylm(3,-3,x,y,z),sqrt(35/(32pi))*(3*x^2-y^2)*y,atol=ɛ)
	@test isapprox(rlylm(3,-2,x,y,z),0+sqrt(105/(4pi))*x*y*z,atol=ɛ)
	@test isapprox(rlylm(3,-1,x,y,z),sqrt(21/(32pi))*y*(4*z^2-x^2-y^2),atol=ɛ)
	@test isapprox(rlylm(3,0,x,y,z),sqrt(7/(16pi))*(5*z^3-3*z*(x^2+y^2+z^2)),atol=ɛ)
	@test isapprox(rlylm(3,1,x,y,z),sqrt(21/(32pi))*x*(4*z^2-x^2-y^2),atol=ɛ)
	@test isapprox(rlylm(3,2,x,y,z),sqrt(105/(16pi))*(x^2-y^2)*z,atol=ɛ)
	@test isapprox(rlylm(3,3,x,y,z),sqrt(35/(32pi))*(x^2-3*y^2)*x,atol=ɛ)

end

@testset "solid harmonic" begin
	ε = 32*eps(Float64)
	@polyvar x y z

	# test for l = 0,1,2,3

	# l = 0
	@test isapprox(SphericalHarmonicExpansions.zlm(0,0,x,y,z),1+0*x,atol=ε)

	# l = 1
	@test isapprox(SphericalHarmonicExpansions.zlm(1,-1,x,y,z),y+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(1,0,x,y,z),z+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(1,1,x,y,z),x+0,atol=ε)

	# l = 2
	@test isapprox(SphericalHarmonicExpansions.zlm(2,-2,x,y,z),sqrt(3)*x*y+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(2,-1,x,y,z),sqrt(3)*y*z+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(2,0,x,y,z),1/2*(2*z^2-x^2-y^2),atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(2,1,x,y,z),sqrt(3)*x*z+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(2,2,x,y,z),sqrt(3/4)*(x^2-y^2),atol=ε)

	# l = 3
	@test isapprox(SphericalHarmonicExpansions.zlm(3,-3,x,y,z),sqrt(5/8)*(3*x^2*y-y^3),atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,-2,x,y,z),sqrt(15)*x*y*z+0,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,-1,x,y,z),sqrt(3/8)*y*(4*z^2-x^2-y^2),atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,0,x,y,z),1/2*(5*z^3-3*z*(x^2+y^2+z^2)),atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,1,x,y,z),sqrt(3/8)*x*(4*z^2-x^2-y^2),atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,2,x,y,z),sqrt(15/4)*(x^2-y^2)*z,atol=ε)
	@test isapprox(SphericalHarmonicExpansions.zlm(3,3,x,y,z),sqrt(5/8)*(x^3-3*x*y^2),atol=ε)
end
