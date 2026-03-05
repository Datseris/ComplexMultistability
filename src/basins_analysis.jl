function intermingledness_matrix(XIC, best_labels)
    if !(hasdim(XIC, InitCond) && length(dims(XIC)) == 2)
        error("need 2D input with init cond and diagnostic")
    end
    ukeys = unique(best_labels)
    # precompute indices of each cluster
    uidxs = [findall(isequal(i), best_labels) for i in ukeys]

    all_inter_cluster_distances = map(eachindex(dims(XIC, Diagnostic))) do j
        # get inter-cluster-distances for one dimension
        x = gnv(XIC[Diagnostic(j)]) # this is a `Vector`
        # separate into the clusters
        xs = [x[is] for is in uidxs]
        # now, for each cluster we need to calculate the mean distance
        inter_cluster_distances = map(eachindex(xs)) do i
            y = xs[i]
            # mean distance to all other clusters (and self)
            ds = mean_distance.(Ref(y), xs)
            # then normalize by intra-cluster distance
            ds = ds[i] ./ ds
            # drop the same cluster cluster entry
            deleteat!(ds, i)
            # and pick the max; arbitrary, could be median or whatever else
            maximum(ds)
        end

        inter_cluster_distances
    end

    matrix_inter_cluster_distances = transpose(hcat(all_inter_cluster_distances...))
    return matrix_inter_cluster_distances
end

function mean_distance(x::AbstractVector, y::AbstractVector)
    c = 0
    d = zero(eltype(x))
    for i in eachindex(x)
        for j in eachindex(y)
            d += @inbounds abs(x[i] - y[j])
            c += 1
        end
    end
    return d/c
end

function intermingledness(initcond::AbstractVector{<:AbstractVector}, basin_ids::AbstractVector{Int})

end