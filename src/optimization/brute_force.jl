using ProgressMeter
using Combinatorics: combinations
using StatsBase: sample

"""
    optimize_feature_selection(objective, allfeatures::StateSpaceSet; kw...)

Brute-force optimize which features selected by `allfeatures` maximize the given `objective` function,
based on the provided `optimizer`. `obj` is [`FeaturesClusteringQuality`](@ref).

Return `best_choice, best_result`, the first being a vector of indices of `allfeatures`
that maximize `obj`, the second being a container of results that includes the maximum,
the attractor labels, and more.

The brute force optimization generates all possible feature combinations by
selecting at least `minF` and at most `maxF`
features. Then it goes through all of them, computes the cost function,
and keeps track of the result with best quality of clustering.

At most `max_choices_per_dim` feature combinations are randomly sampled
for a given feature dimension to use to limit computational time.

If `reclustering > 0` the clustering process is applied iteratively again on each found
cluster `reclustring` amount of times.

## Keyword arguments
```
maxF::Int = 3
minF::Int = 1
max_choices_per_dim::Int = typemax(Int)
verbose::Bool = false
reclustering::Int = 0
```
"""
function optimize_feature_selection(
        features_clustering_quality, allfeatures::StateSpaceSet;
        maxF::Int = 3,
        minF::Int = 1,
        max_choices_per_dim::Int = typemax(Int),
        verbose::Bool = false,
        reclustering::Int = 0,
    )
    d = input_data_size(allfeatures)[1]
    diag_idxs = 1:d
    alloptions = generate_choices(diag_idxs, minF, maxF, max_choices_per_dim)
    best_choice = Int[]
    max_quality = -Inf
    best_result = FeaturesClusteringResults()
    ProgressMeter.@showprogress desc="Optimizing clustering..." for choice in alloptions
        features = allfeatures[:, choice]
        result = features_clustering_quality(features)
        quality = result.quality
        if verbose
            println("choice = ", choice, ", quality = ", quality, ", v_optimal = $(result.v_optimal)", ", L = $(result.outliers)")
        end
        if quality > max_quality
            max_quality = quality
            best_choice = choice
            best_result = result
        end
    end
    best_labels = best_result.labels
    # a final iterative reclustering step for the case of nested groups
    recluster!(best_labels, allfeatures, features_clustering_quality, reclustering; verbose = verbose)
    return best_choice, best_labels
end

function recluster!(best_labels, allfeatures, fcq, N::Int; verbose = false)
    for _ in 1:N
        ulabels = unique(best_labels)
        next_label = maximum(best_labels) + 1
        for existing_label in ulabels
            label_idxs = findall(isequal(existing_label), best_labels)
            existing_label_features = allfeatures[label_idxs]
            length(existing_label_features) ≤ dimension(existing_label_features) && continue
            cr = cluster(fcq.ca, existing_label_features)
            reclustered_labels = cr.labels
            L = length(unique(reclustered_labels))
            L == 1 && continue
            # we found a case where a found cluster further differentiates more
            # so the new label `1` gets assugned the `existing_label`, while the
            # rest obtain labels from the `next_label` counter
            if verbose
                println("Reclustering of original label $(existing_label) found total of $(L) labels.")
            end
            new_labels_mapping = [1 => existing_label]
            for jj in 2:L
                push!(new_labels_mapping, jj => next_label)
                next_label += 1
            end
            replace!(reclustered_labels, new_labels_mapping...)
            # and finally update main container
            best_labels[label_idxs] .= reclustered_labels
        end
    end
    return
end


# all combinatorially possible choices with minimum and maximum feature dimension
function generate_choices(diag_idxs, minF, maxF, max_choices_per_dim)
    alloptions = Vector{Int}[]
    for n in minF:maxF
        combs = SVector{n, Int}.(collect(combinations(diag_idxs, n)))
        if length(combs) ≤ max_choices_per_dim
            combs_reduced = combs
        else
            # combinations are too many, we reduce them by random sampling
            combs_reduced = sample(combs, max_choices_per_dim; replace = false)
        end
        append!(alloptions, combs_reduced)
    end
    return alloptions
end
