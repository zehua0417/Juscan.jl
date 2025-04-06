#=========================================
filename: colors
author: Lihuax
date: 2025-04-05 12:41
description:
=========================================#

const COLORS_DISCRETE_FRIENDLY = [
  colorant"#0072B2",
  colorant"#56B4E9",
  colorant"#009E73",
  colorant"#F5C710",
  colorant"#E69F00",
  colorant"#D55E00",
]
const COLORS_DISCRETE_SEASIDE =
  [colorant"#8ecae6", colorant"#219ebc", colorant"#023047", colorant"#ffb703", colorant"#fb8500"]
const COLORS_DISCRETE_FRIENDLY_LONG = [
  colorant"#CC79A7",
  colorant"#0072B2",
  colorant"#56B4E9",
  colorant"#009E73",
  colorant"#F5C710",
  colorant"#E69F00",
  colorant"#D55E00",
]
const COLORS_DISCRETE_APPLE = [
  colorant"#ff3b30",
  colorant"#ff9500",
  colorant"#ffcc00",
  colorant"#4cd964",
  colorant"#5ac8fa",
  colorant"#007aff",
  colorant"#5856d6",
]
const COLORS_DISCRETE_IBM =
  [colorant"#5B8DFE", colorant"#725DEE", colorant"#DD227D", colorant"#FE5F00", colorant"#FFB109"]
const COLORS_DISCRETE_RAINBOW = [
  colorant"#FF7777",
  colorant"#FFAB74",
  colorant"#FFE577",
  colorant"#DBF47B",
  colorant"#91E480",
  colorant"#7CC9E5",
  colorant"#7DA8E6",
  colorant"#887DE6",
  colorant"#BC7BE4",
]
const COLORS_DISCRETE_CANDY =
  [colorant"#9b5de5", colorant"#f15bb5", colorant"#fee440", colorant"#00bbf9", colorant"#00f5d4"]

const COLORS_DISCRETE_ALGER =
  [colorant"#000000", colorant"#1A5B5B", colorant"#ACC8BE", colorant"#F4AB5C", colorant"#D1422F"]

const COLOR_SCHEMES = Dict(
  "friendly" => COLORS_DISCRETE_FRIENDLY,
  "seaside" => COLORS_DISCRETE_SEASIDE,
  "friendly_long" => COLORS_DISCRETE_FRIENDLY_LONG,
  "apple" => COLORS_DISCRETE_APPLE,
  "ibm" => COLORS_DISCRETE_IBM,
  "rainbow" => COLORS_DISCRETE_RAINBOW,
  "candy" => COLORS_DISCRETE_CANDY,
  "alger" => COLORS_DISCRETE_ALGER,
)
