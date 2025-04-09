#=========================================
filename: plots
author: Lihuax
date: 2025-04-04 21:16
description:
=========================================#

"""
    expand_palette(base_colors::Vector{<:Colorant}, n::Int)

Return `n` interpolated colors based on the `base_colors` using `cgrad`.
"""
function expand_palette(
  base_colors::Vector{<:Colorant},
  n::Int;
  downsample_strategy::String="evenly",
)
  if n == length(base_colors)
    return base_colors
  elseif n < length(base_colors)
    return downsample_vector(base_colors, n; strategy=downsample_strategy)
  else
    cs = cgrad(base_colors, n, categorical=true)
    return Colorant[c for c in cs.colors]
  end
end

#==============
"""
    downsample_vector(x::Vector, n::Int; strategy::String="evenly")

Return a vector of length `n` downsampled from `x` using the given `strategy`.

Supported strategies:
- "evenly"
- "first"
- "last"
- "middle"
"""
=============#
function downsample_vector(x::Vector, n::Int; strategy::String="evenly")
  len = length(x)
  len <= n && return x

  if strategy == "evenly"
    idx = round.(Int, range(1, len, length=n))
    return x[idx]
  elseif strategy == "first"
    return x[1:n]
  elseif strategy == "last"
    return x[(end - n + 1):end]
  elseif strategy == "middle"
    start_idx = Int(floor((len - n) / 2)) + 1
    return x[start_idx:(start_idx + n - 1)]
  else
    error("Unknown downsample strategy: $strategy")
  end
end

function get_palette(name::String, n::Int; downsample_strategy::String="evenly")
  base = haskey(COLOR_SCHEMES, name) ? COLOR_SCHEMES[name] : error("Palette '$name' not found")
  return expand_palette(base, n; downsample_strategy=downsample_strategy)
end

"""
    get_continuous_colormap(name::String, n::Int=265)

Get `n` evenly spaced colors from a predefined continuous colormap in `ColorSchemes.jl`.
Available names: :viridis, :inferno, :plasma, :magma, :turbo, :cividis, etc.
"""
function get_continuous_colormap(name::String, n::Int=265)
  cs = getfield(ColorSchemes, Symbol(name))
  return get(cs, range(0, 1; length=n))
end

"""
    violin(adata::AnnData, keys; kwargs...)

Draws violin plots for one or more features stored in `adata.obs`.

# Arguments
- `adata::AnnData`: Annotated data object.
- `keys::Union{String, Vector{String}}`: One or more feature names from `adata.obs` to plot.

# Keyword Arguments
- `width::Real=600`: Width of each subplot.
- `height::Real=400`: Height of the plot.
- `jitter::Union{Bool, Real}=0.5`: Jitter width or `true` to apply default jitter.
- `dot_size::Real=2`: Size of jittered dots.
- `title::String="violin plot"`: Title of the full figure.
- `palette_name::String="friendly"`: Color palette name.
- `fill_alpha::Real=1.0`: Transparency of the violin fill.
- `downsample_strategy::String="evenly"`: Strategy for downsampling color palette.
- `savefig::Union{Bool, String}=false`: Whether to save the figure, or path to save.

# Returns
- A `CairoMakie.Figure` or `true` if `savefig=true`.

"""
function violin(
  adata::Muon.AnnData,
  keys::Union{String, Vector{String}};
  width::Real=600,
  height::Real=400,
  jitter::Union{Bool, Real}=0.5,
  dot_size::Real=2,
  title::String="violin plot",
  palette_name::String="friendly",
  fill_alpha::Real=1.0,
  downsample_strategy::String="evenly",
  savefig::Union{Bool, String}=false,
)
  keys = isa(keys, String) ? [keys] : keys
  n = length(keys)

  fig = Figure(size=(width * n, height))
  palette = get_palette(palette_name, n; downsample_strategy=downsample_strategy)

  for (i, key) in enumerate(keys)
    data = adata.obs[:, key]
    ax = Axis(fig[1, i], xticks=([1], [""]), ylabel=key, xtrimspine=true)

    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.rightspinevisible = false
    ax.topspinevisible = false
    CairoMakie.xlims!(ax, 0.5, 1.5)
    CairoMakie.ylims!(
      ax,
      minimum(data) - 0.05 * abs(minimum(data)),
      maximum(data) + 0.05 * abs(maximum(data)),
    )

    CairoMakie.violin!(
      ax,
      fill(1, length(data)),
      data;
      width=0.6,
      side=:both,
      show_median=true,
      color=coloralpha(palette[i], fill_alpha),
    )

    if jitter != false
      jwidth = jitter === true ? 0.3 : jitter
      x_jittered = 1 .+ jwidth .* (rand(length(data)) .- 0.5)
      CairoMakie.scatter!(ax, x_jittered, data; color=RGBAf(0, 0, 0, 0.7), markersize=dot_size)
    end
  end

  Label(fig[0, :], title, fontsize=18, tellwidth=false)

  if savefig != false
    path = savefig === true ? "violin_plot.png" : String(savefig)
    CairoMakie.save(path, fig)
    return true
  end

  return fig
end

"""
    scatter(adata::AnnData, x_key, y_key; kwargs...)

Generates a scatter plot from two features in `adata.obs`, optionally colored by a third.

# Arguments
- `adata::AnnData`: Annotated data object.
- `x_key::String`: Feature name for x-axis.
- `y_key::String`: Feature name for y-axis.

# Keyword Arguments
- `color_key::Union{String, Nothing}=nothing`: Feature used for point colors.
- `title::String="scatter plot"`: Plot title.
- `width::Real=600`: Width of the plot.
- `height::Real=400`: Height of the plot.
- `colormap_name::String="viridis"`: Colormap for `color_key`.
- `ncolors::Int=265`: Number of colors in colormap.
- `downsample_strategy::String="evenly"`: Strategy for colormap downsampling.
- `savefig::Union{Bool, String}=false`: Whether to save the figure.

# Returns
- A `Figure` or `true` if `savefig=true`.

"""
function scatter(
  adata::Muon.AnnData,
  x_key::String,
  y_key::String;
  width::Real=600,
  height::Real=400,
  color_key::Union{String, Nothing}=nothing,
  title::String="scatter plot",
  savefig::Union{Bool, String}=false,
  colormap_name::String="viridis",
  ncolors::Int=265,
  downsample_strategy::String="evenly",
)
  fig = Figure(size=(width, height))
  ax = Axis(fig[1, 1], xlabel=x_key, ylabel=y_key, title=title)
  x_data = adata.obs[:, x_key]
  y_data = adata.obs[:, y_key]

  if color_key !== nothing
    color_data = adata.obs[:, color_key]
    cmap = get_continuous_colormap(colormap_name, ncolors)
    sc = CairoMakie.scatter!(
      ax,
      x_data,
      y_data;
      color=color_data,
      colormap=cmap,
      colorrange=(minimum(color_data), maximum(color_data)),
    )
    Colorbar(fig[1, 2], sc; label=color_key)
  else
    CairoMakie.scatter!(ax, x_data, y_data)
  end

  if savefig != false
    path = savefig === true ? "scatter_plot.png" : String(savefig)
    CairoMakie.save(path, fig)
    return true
  end
  return fig
end

"""
    hvg_scatter(adata::AnnData; savefig=false)

Plots mean expression vs dispersion for genes, highlighting highly variable genes (HVGs).

# Arguments
- `adata::AnnData`: Annotated data object with `highly_variable`, `means`, `variances`, and `variances_norm` in `adata.var`.

# Keyword Arguments
- `savefig::Union{Bool, String}=false`: Whether to save the figure.

# Returns
- A `Figure` or `true` if `savefig=true`.
"""
function hvg_scatter(adata::AnnData; savefig::Union{Bool, String}=false)
  df = adata.var
  means = df.means
  variances = df.variances
  variances_norm = df.variances_norm
  hvg_mask = df.highly_variable

  x = log10.(means .+ 1e-6)
  y_norm = variances_norm
  y_raw = log10.(variances .+ 1e-6)

  fig = Figure(size=(1000, 400))
  ax1 = Axis(fig[1, 1], title="Normalized dispersion")
  ax2 = Axis(fig[1, 2], title="Raw dispersion")

  # Gray“other genes”
  Makie.scatter!(
    ax1,
    x[.!hvg_mask],
    y_norm[.!hvg_mask];
    color=:gray,
    markersize=3,
    transparency=true,
  )
  Makie.scatter!(
    ax2,
    x[.!hvg_mask],
    y_raw[.!hvg_mask];
    color=:gray,
    markersize=3,
    transparency=true,
  )

  # Black HVGs
  Makie.scatter!(ax1, x[hvg_mask], y_norm[hvg_mask]; color=:black, markersize=3)
  Makie.scatter!(ax2, x[hvg_mask], y_raw[hvg_mask]; color=:black, markersize=3)

  ax1.xlabel = "mean expressions of genes"
  ax1.ylabel = "dispersions of genes (normalized)"

  ax2.xlabel = "mean expressions of genes"
  ax2.ylabel = "dispersions of genes (not normalized)"

  if savefig != false
    path = savefig === true ? "hvg_scatter_plot.png" : String(savefig)
    CairoMakie.save(path, fig)
    return true
  end

  fig
end

"""
    plot_variance_ratio(adata::AnnData; key="pca_variance", n=50, savefig=false)

Plots the explained variance ratio for the top principal components (PCA).

# Arguments
- `adata::AnnData`: Annotated data object.
- `key::AbstractString="pca_variance"`: Key in `adata.uns` for PCA singular values.
- `n::Real=50`: Number of PCs to plot.

# Keyword Arguments
- `savefig::Union{Bool, String}=false`: Whether to save the figure.

# Returns
- A `Figure` or `true` if `savefig=true`.

"""
function plot_variance_ratio(
  adata::AnnData;
  key::AbstractString="pca_variance",
  n::Real=50,
  savefig::Union{Bool, String}=false,
)
  s = adata.uns[key]
  var_ratio = (s .^ 2) ./ sum(s .^ 2)
  log_var = log10.(var_ratio[1:n])

  fig = Figure(size=(500, 500))
  ax = Axis(fig[1, 1], xlabel="ranking", ylabel="variance ratio", title="variance ratio")

  xs = 1:n
  Makie.scatter!(ax, xs, log_var; color=:black, markersize=4)
  lines!(ax, xs, log_var; color=:black)

  for i in 1:5
    text!(ax, "PC$i", position=(xs[i], log_var[i]), align=(:left, :center), fontsize=9)
  end

  if savefig != false
    path = savefig === true ? "variance_ratio_plot.png" : String(savefig)
    CairoMakie.save(path, fig)
    return true
  end

  fig
end

"""
    plot_umap(adata; kwargs...)

Plots UMAP embedding colored by a categorical label in `adata.obs`.

# Keyword Arguments
- `color_by::String="clusters_0.6"`: Observation field to color points.
- `key::String="umap"`: Key in `adata.obsm` containing UMAP coordinates.
- `palette_name::String="rainbow"`: Color palette name.
- `width::Real=800`: Width of the plot.
- `height::Real=600`: Height of the plot.
- `downsample_strategy::String="evenly"`: Color downsampling strategy.
- `savefig::Union{Bool, String}=false`: Save to file if `true` or path given.

# Returns
- A `Figure` or `true` if `savefig=true`.

# Throws
- Error if `umap_coords` is not found in `adata.obsm`.

"""
function plot_umap(
  adata;
  color_by::String="clusters_latest",
  key="umap",
  palette_name::String="rainbow",
  width::Real=800,
  height::Real=600,
  downsample_strategy::String="evenly",
  savefig::Union{Bool, String}=false,
)
  umap_coords = adata.obsm[key]
  if umap_coords === nothing
    error("UMAP coordinates not found. Please run umap!(adata) first.")
  end

  labels = adata.obs[!, color_by]
  unique_labels = sort(unique(labels))
  n_labels = length(unique_labels)

  cmap = get_palette(palette_name, n_labels; downsample_strategy=downsample_strategy)

  label_to_idx = Dict(k => i for (i, k) in enumerate(unique_labels))
  color_indices = [label_to_idx[l] for l in labels]

  fig = Figure(size=(width, height))
  ax = Axis(fig[1, 1], title="UMAP colored by $color_by", xlabel="UMAP 1", ylabel="UMAP 2")

  Makie.scatter!(
    ax,
    umap_coords[:, 1],
    umap_coords[:, 2];
    color=cmap[color_indices],
    markersize=6,
    strokewidth=0.2,
  )

  for (idx, label) in enumerate(unique_labels)
    Makie.scatter!(ax, [NaN], [NaN]; color=cmap[idx], label="Cluster $label")
  end
  axislegend(ax; position=:rb)

  if savefig != false
    path = savefig === true ? "umap_plot.png" : String(savefig)
    CairoMakie.save(path, fig)
    return true
  end

  fig
end
