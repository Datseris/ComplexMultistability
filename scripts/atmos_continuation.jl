using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
include("atmos_shared.jl")

# %% computation across params

epsilons = [0.8, 0.85, 0.9, 0.95]

# we have defined in the shared file the `extract_features` so now we loop

labels_across_param = []
features_across_param = []
best_features_across_param = []

ca = ADBSCAN(; min_neighbors = 10)
opt = BruteForce(verbose = false, minF = 2, reclustering = 0, )
fcq = FeaturesClusteringQuality(; attractor_weight = 1.0, feature_weight = 0.5, ca)

for epsilon in epsilons
    filename = "MAOOAM_ep_$(epsilon).nc"
    ncfile = datadir("oisin", filename)
    X = ncread(ncfile, "value")
    t = gnv(dims(X, Tim))
    sampler = 10
    X = X[Tim(1:sampler:length(t))]
    t = t[1:sampler:length(t)]
    diagnostics = gnv(dims(X, Diagnostic))
    ics = gnv(dims(X, InitCond)) .+ 1
    allfeatures = extract_features(X)
    best_choice, best_labels = optimize_feature_selection(fcq, allfeatures;
        verbose = false, minF = 2, reclustering = 0,
    )
    push!(labels_across_param, best_labels)
    push!(best_features_across_param, best_choice)
    push!(features_across_param, allfeatures)
end

# %% matching (i.e., continuation)
# okay now we do the matching step, i.e., continuation, choosing the best features
# as the state space sets, and their centroids as the matching distance


using Attractors
using Statistics: mean

# cast data into "attractors": Dictionary mapping IDs to a State Space Set
# whose rows are the different features belonging in the same cluster.
function cast_to_attractor_dict(features::StateSpaceSet, labels)
    ukeys = unique(labels)
    gs = [findall(isequal(clu_idx), labels) for clu_idx in ukeys]
    return Dict(i => features[gs[i]] for i in ukeys)
end

atts_across_params = [unique(v) for v in labels_across_param]
attractors_cont = [cast_to_attractor_dict(f, l) for (f, l) in zip(features_across_param, labels_across_param)]

fidx = 27 # feature 27 is a common best feature
function mydistance(A, B)
    # cast to feature dimension
    a = A[:, fidx]; b = B[:, fidx]
    return abs(mean(a) - mean(b))
    return d
end

matcher = MatchBySSSetDistance(distance = mydistance)

rmaps = match_sequentially!(attractors_cont, matcher)

# matching worked fine but we were unlucky and got different labels from the labels in the
# previous figure. So we just adjust the `rmap` and rematch with the adjusted one.
consistency_rmap = Dict(1 => 2, 2 => 1)
for a in attractors_cont
    swap_dict_keys!(a, consistency_rmap)
end

# %% plot

# NOTE FOR READER: In this code I did not properly apply the tracking to the
# continuation of intermingledness and assummed (correctly) that the order of the
# attractors coincides with their labels. This allowed me to use the
# row number also as the attractor key. Normally you want to cast everything
# into a dictionary form mapping attractor ID to intermingledness value.
# This is done in the Attractors.jl implementation of intermingledness for example.

# continue intermingledness
intermingledness_cont = intermingledness.(attractors_cont)

fig, axs = axesgrid(2, 1; sharex = true, size = (figwidth÷2, figheight), xlabels = "ε", ylabels = [diagnostics[chosen], "mean intermingl."])
# plot attractors cont
chosen = 27 # which diagnostic to project to
a2r = A -> mean(A[:, chosen])
plot_attractors_curves!(axs[1], attractors_cont, a2r, epsilons; series_kwargs = (markersize = 20, linewidth = 1, ), add_legend = false)

# plot intermingledness cont, which will be average across all dims
inter_info = map(1:length(attractors_cont)) do i
    Dict(k => mean(intermingledness_cont[i][k, :]) for k in keys(attractors_cont[i]))
end
plot_continuation_curves!(axs[2], inter_info, epsilons; add_legend = false, series_kwargs = (markersize = 20, linewidth = 1, ))
# hideydecorations!.(axs; label = false)
xlims!(axs[end], nothing, nothing)
hideydecorations!(axs[1], label = false)
fig

# %%

wsave(papersdir("figures", "atmos_continuation"), fig)
