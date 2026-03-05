include(srcdir("data_loading.jl"))
include(srcdir("clustering", "clustering_api.jl"))
include(srcdir("optimization", "optimization_api.jl"))
include(srcdir("intermingledness.jl"))

file = datadir("lohmann_veros_amoc", input_file_name)

X = ncread(file, "values")
diagnostics = gnv(dims(X, Diagnostic))
ics = gnv(dims(X, InitCond))
ttotal = gnv(dims(X, Ti))

# some data are missing, not all initial conditions are the same length.
# Identify the time that all data are nonmissing
Xvalid = select_in_common_time(X)
t = gnv(dims(Xvalid, Ti))

# We do the clustering of Veros again anyways as we need it to demonstrate the algorithm
# First create features
tstart = 19000
Y = Xvalid[Tim(Between(tstart, t[end]))]
using Statistics: mean
function allfeaturizer(A::StateSpaceSet, t) # each initial condition gives one `A`
    cols = columns(A)
    feats = map(c -> mean(c), cols)
    return SVector(feats)
end
allfeatures = StateSpaceSet(cast_to_features(Y, allfeaturizer))
# Then do the clustering
ca = ADBSCAN(; rescale_features = true, min_neighbors = 4)
fcq = FeaturesClusteringQuality(; attractor_weight = 1.0, feature_weight = 0.5, ca)
best_choice, best_labels = optimize_feature_selection(fcq, allfeatures;
    verbose = false, max_choices_per_dim = typemax(Int)
)

# Re-cluster only features with index 1
label1idxs = findall(isequal(1), best_labels)
allfeatures1 = allfeatures[label1idxs]
cluster1_choice, cluster1_labels = optimize_feature_selection(opt, fcq, allfeatures1)

# This indeed leads to 2 attractors, so we normalize the existing labels
cluster1_labels[findall(isequal(2), cluster1_labels)] .= 5
best_labels[label1idxs] .= cluster1_labels
attractors = unique(best_labels)
natt = length(attractors)