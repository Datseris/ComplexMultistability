# uses clustering api
using Statistics: mean

"""
    FeaturesClusteringQuality(;
        attractor_weight = 1.0
        feature_weight = 0.5
        punish_weight = 0.1
        gvc = GroupViaClustering(),
        verbose = false,
    ) → fcp

A type that can be used as a function to estimate the quality of a clustering of a
state space set of feature-vectors. `fcq(features)` returns
the clustering quality, given by the equation:

```math
q = v*(w_A A - w_F F - w_P P)
```
where ``v`` is the original clustering quality metric as returned by
ADBSCAN, ``A`` the number of attractors, ``F`` the number of features,
``P`` is the number of feature vectors that couldn't
be clustered and ``w_x`` the associated weights.

Typically an instance `fcp` is given directly to [`optimize_feature_selection`](@ref).
"""
@kwdef struct FeaturesClusteringQuality{C<:ClusteringAlgorithm}
    attractor_weight::Float64 = 1.0 # scaling of found attractors
    feature_weight::Float64 = 0.5 # scaling of used feature dimensions
    punish_weight::Float64 = 0.1 # failure punishment
    # Options for Attractors.jl and clustering
    ca::C = ADBSCAN()
end

@kwdef struct FeaturesClusteringResults
    # These two are for any type of grouping into attractors
    quality::Float64 = 0.0
    labels::Vector{Int} = Int[]
    # The rest are specific to the particular type
    v_optimal::Float64 = 0.0
    outliers::Int = 0
    attractor_weight::Float64 = 0.0
    feature_weight::Float64 = 0.0
    punish_weight::Float64 = 0.0
end

function (fcq::FeaturesClusteringQuality)(features)
    cr = cluster(fcq.ca, features)
    v_optimal = cr.v_optimal
    labels = cr.labels
    ulabels = unique(labels)
    if -1 ∈ ulabels
        A = length(ulabels)-1
        L = count(<(0), labels)
    else
        A = length(ulabels)
        L = 0
    end
    feature_dim = input_data_size(features)[1]
    quality = v_optimal*(1 + fcq.attractor_weight*A - fcq.feature_weight*feature_dim - fcq.punish_weight*L)
    return FeaturesClusteringResults(
        quality, labels, v_optimal, L,
        fcq.attractor_weight, fcq.feature_weight, fcq.punish_weight
    )
    return quality
end
