using Test
using EVMavg

x = [1.0,2.0,3.0,4.0]
s = [0.1,0.2,0.3,0.4]

@testset "Test EVMavg: avg" begin

	uav,wav,mer = avg(x,s)

	@test uav ≈ 2.5 atol=1E-14
    @test wav ≈ 1.4634146341463417 atol=1E-14
    @test mer ≈ 0.941876146945452 atol=1E-14

end

@testset "Test EVMavg: evm" begin

	xb, se_int, se_ext = evm(x,s)

	@test xb ≈ 1.9269849828944758 atol=1E-14
    @test se_int ≈ 0.096357574983269 atol=1E-14
    @test se_ext ≈ 1.0559929440171045 atol=1E-14

end
