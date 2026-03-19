# This script unrolls the main loop of `optimize_feature_selection`
# to store the results so that we can see how all things behave
# # Load data and setup
using DrWatson
@quickactivate "ComplexMultistability"
include(srcdir("data_loading.jl"))
include(srcdir("clustering", "clustering_api.jl"))
include(srcdir("optimization", "optimization_api.jl"))

file = datadir("lohmann_veros_amoc", "Veros_Data_F8_1.nc")
X = ncread(file, "values")
diagnostics = gnv(dims(X, Diagnostic))
ics = gnv(dims(X, InitCond))
ttotal = gnv(dims(X, Ti))
Xvalid = select_in_common_time(X)
t = gnv(dims(Xvalid, Ti))
tstart = 19000
Y = Xvalid[Tim(Between(tstart, t[end]))]
using Statistics: mean
function allfeaturizer(A::StateSpaceSet, t) # each initial condition gives one `A`
    cols = columns(A)
    feats = map(c -> mean(c), cols)
    return SVector(feats)
end
allfeatures = cast_to_features(Y, allfeaturizer)
ca = ADBSCAN(; rescale_features = true, min_neighbors = 4)
fcq = FeaturesClusteringQuality(; attractor_weight = 1.0, feature_weight = 0.5, ca)

# Optimization loop
maxF::Int = 4
minF::Int = 1
max_choices_per_dim::Int = 2000 # keep things bounded for now

d = input_data_size(allfeatures)[1]
diag_idxs = 1:d
alloptions = generate_choices(diag_idxs, minF, maxF, max_choices_per_dim)
best_choice = Int[]
stored_choices = []
stored_optimals = []
stored_qualities = []
stored_attractors = []

ProgressMeter.@showprogress desc="BruteForce Clustering..." for choice in alloptions
    features = allfeatures[:, choice]
    result = fcq(features)
    s = result.v_optimal
    q = result.quality
    push!(stored_choices, choice)
    push!(stored_optimals, s)
    push!(stored_qualities, q)
    push!(stored_attractors, count(>(0), unique(result.labels)))
end

# %% analyze the results
include(srcdir("papertheme.jl"))

fig = Figure(size = (figwidth÷2, figheight))
ax = Axis(fig[1,1]; ylabel = "silhouette mean (s)", xlabel = "# groups created")
markers = [MARKERS[k] for k in length.(stored_choices)]

x_locs = stored_attractors .+ ((length.(stored_choices) .- 3) ./ 10)

scatter!(ax, x_locs, stored_optimals; marker = markers, color = (:black, 0.1))

Legend(fig[1,2],
    [MarkerElement(marker = MARKERS[k], color = :black, markersize = 20) for k in 1:4],
    string.(1:4),
    L"F",
)

display(fig)

# %%

wsave(papersdir("figures", "dbscan_test"), fig)
