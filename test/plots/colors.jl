using Colors
using Juscan.Pl

@testset "expand_palette" begin
  short = expand_palette([colorant"#000000", colorant"#FFFFFF"], 1)
  @test length(short) == 1
  @test short[1] == colorant"#000000"

  expanded = expand_palette([colorant"#000000", colorant"#FFFFFF"], 5)
  @test length(expanded) == 5
  @test eltype(expanded) <: Colorant
end

@testset "get_palette" begin
  pal = get_palette("friendly", 4)
  @test length(pal) == 4
  @test eltype(pal) <: Colorant

  pal2 = get_palette("friendly", 10)  # 插值
  @test length(pal2) == 10

  @test_throws ErrorException get_palette("nonexistent", 5)
end

@testset "get_continuous_colormap" begin
  viridis = get_continuous_colormap("viridis", 265)
  @test length(viridis) == 265
  @test eltype(viridis) <: Colorant

  turbo = get_continuous_colormap("turbo", 50)
  @test length(turbo) == 50
  @test isa(turbo[1], Colorant)
end
