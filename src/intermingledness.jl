# Note: Attractors.jl implementation is better.
# First, it accepts a generic distance function.
# Second, it returns a dictionary instead of a matrix.
# (as it is done for only one distance function, simplifying the API.
# it can be looped at a higher level for more distance functions)

"""
    intermingledness(u0s::StateSpaceSet, labels)

Given a `StateSpaceSet` of "initial conditions", or features, or anything,
and the labels corresponding to what group ID each `u0` belongs in,
return the measure of intermingledness.
It is returned as a matrix, with as columns as the dimension of `u0s` and
as many rows as attractors, since intermingledness is calculated basin.

Intermingledness is normalized by the mean distance in each cluster,
the closer to 1, the more intermingled!
"""
function intermingledness(u0s::StateSpaceSet, labels;
        summarizer = maximum
    )
    ukeys = unique(labels)
    imetric_per_dim = map(1:dimension(u0s)) do dim_idx
        x = u0s[:, dim_idx]
        # separate into the clusters
        xs = [x[findall(isequal(clu_idx), labels)] for clu_idx in ukeys]
        imetric = map(eachindex(xs)) do clu_idx
            x = xs[clu_idx]
            # mean distance of current claster to all clusters (and self)
            ds = mean_distance.(Ref(x), xs)
            # the metric is now the distance of belonging cluster divided
            # by distance to any other cluster
            imetrics = ds[clu_idx] ./ ds
            # Right, but now we still need to return a summarizing number across
            # all clusters. First we drop the same cluster cluster entry (which is 1)
            deleteat!(imetrics, clu_idx)
            # and then summarize
            return summarizer(imetrics)
        end
        return imetric
    end
    return hcat(imetric_per_dim...)
end

function mean_distance(x::AbstractVector, y::AbstractVector)
    c = 0
    d = zero(eltype(x))
    for i in eachindex(x)
        for dim_idx in eachindex(y)
            d += @inbounds abs(x[i] - y[dim_idx])
            c += 1
        end
    end
    return d/c
end


"""
    intermingledness(groups::Dict{Int, <:StateSpaceSet})

Helper method that allows input a dictionary mapping group IDs to initial conditions
or features belonging to that group.
"""
function intermingledness(attractors::Dict{Int, <:StateSpaceSet})
    u0s = eltype(first(values(attractors)))[]
    labels = Int[]
    ukeys = sort!(collect(keys(attractors)))
    for j in ukeys
        A = attractors[j]
        append!(u0s, vec(A))
        append!(labels, fill(j, length(A)))
    end
    intermingledness(StateSpaceSet(u0s), labels)
end

using Statistics, CairoMakie
function plot_intermingledness(im; kw...)
    fig = Figure()
    plot_intermingledness!(fig, im; kw...)
    return fig
end
function plot_intermingledness!(fig, intermingledness;
        cmap = to_color.([:white, :black]), ilims = (0, 1),
        names = string.(axes(intermingledness, 2)),
        unitlen = 40, hidexlabels = false, fixed_size = true,
        add_colorbar = true, xlabel = "diagnostic (state space dimension)",
    )

    N = size(intermingledness, 1) # number of basins
    D = size(intermingledness, 2) # number of diagnostics

    # setup colors
    hmapkw = (colorrange = ilims, colormap = cgrad(cmap, 9, categorical = true), highclip = cmap[end])
    axbgkw = (xgridvisible = false, ygridvisible = false, xminorgridvisible = true,
        yminorgridvisible = true, xminorgridcolor = :grey, yminorgridcolor = :grey,
        xminorticks = IntervalsBetween(2), yminorticks = IntervalsBetween(2),
    )

    # Heatmap of main data
    axhm = Axis(fig[2,1];
        width = fixed_size ? unitlen*D : Auto(),
        height = fixed_size ? unitlen * N : Auto(), axbgkw...,
        ylabel = "group ID", xlabel,
    )
    hm = heatmap!(axhm, intermingledness'; hmapkw...)
    translate!(hm, 0, 0, -100)
    axhm.xticks = (1:length(names), names)
    axhm.xticklabelrotation = π/4
    axhm.yticks = 1:N
    axhm.xgridvisible = false
    # heatmap of means along rows: averages for each basin
    hmeans = mean(intermingledness; dims = 2)
    ax2 = Axis(fig[2,2];
        width = fixed_size ? unitlen : Auto(),
        height = fixed_size ? unitlen * N : Auto(),
    axbgkw...,)
    ax2.xticks = ([1], ["mean"])
    ax2.xticklabelrotation = π/4
    ax2.yticks = 1:N
    ax2.yticklabelsvisible = false
    hm = heatmap!(ax2, hmeans'; hmapkw...)
    translate!(hm, 0, 0, -100)
    # heatmap of means along columns: averages for each diagnostic
    mean_per_diagnostic = mean(intermingledness; dims = 1)
    ax3 = Axis(fig[1,1];
        width = fixed_size ? unitlen * D : Auto(),
        height = fixed_size ? unitlen : Auto(), axbgkw...,
    )
    hm = heatmap!(ax3, mean_per_diagnostic'; hmapkw...)
    translate!(hm, 0, 0, -100)
    ax3.yticks = ([1], ["mean"])
    ax3.xticks = 1:D
    ax3.xticklabelsvisible = false
    # Layouting
    layout = fig isa Figure ? fig.layout : fig
    rowgap!(layout, 5)
    colgap!(layout, 5)
    if add_colorbar
        Colorbar(fig[:,3]; label = "intermingledness", width = 30, hmapkw...)
        colgap!(layout, 2, 15)
    end
    if hidexlabels
        axhm.xticklabelsvisible = false
        axhm.xlabelvisible = false
        ax3.xticklabelsvisible = false
        ax2.xticklabelsvisible = false
    end
    return fig
end

function plot_twice_intermingledness!(fig, imatrix1, imatrix2;
        title1 = "", title2 = "", cmap = to_color.([:white, :black]), ilims = (0, 1),
        names = string.(axes(intermingledness, 2)), names2 = names, kw...
    )

    fig1 = plot_intermingledness!(GridLayout(fig[1,1]), imatrix_features;
        names, kw..., cmap, ilims, hidexlabels = names2 == names,
        xlabel = "feature (extracted from diagnostics)", add_colorbar = false,
    )
    content(fig1[1,1]).title = title1
    fig2 = plot_intermingledness!(GridLayout(fig[2,1]), imatrix_basins;
        kw..., names = names2, cmap, ilims, add_colorbar = false,
    )
    content(fig2[1,1]).title = title2
    hmapkw = (colorrange = ilims, colormap = cgrad(cmap, 9, categorical = true), highclip = cmap[end])
    Colorbar(fig[:,2]; label = "intermingledness", width = 20, ticks = [0,1], labelpadding = -20, hmapkw...)
    rowgap!(fig, 1, 5)
    return
end
