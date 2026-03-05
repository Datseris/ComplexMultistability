include(srcdir("data_loading.jl"))
include(srcdir("clustering", "clustering_api.jl"))
include(srcdir("optimization", "optimization_api.jl"))
include(srcdir("intermingledness.jl"))
using Statistics

# Change theme to paper
include(srcdir("papertheme.jl"))

# %% Load data
filename = "MAOOAM_ep_0.9.nc"
# filename = "epsilon=0.9.nc"
ncfile = datadir("oisin", filename)
ncdetails(ncfile)

X = ncread(ncfile, "value")
t = gnv(dims(X, Tim))
# right, this is literally insane amount of density in time sampling, so we subsample
sampler = 10
X = X[Tim(1:sampler:length(t))]
t = t[1:sampler:length(t)]
# Furthermore the diagnostics are so large that are not possible to plot!
# I reduce them (after finding out the best ones to make sure they are kept)
diagnostics = gnv(dims(X, Diagnostic))
ics = gnv(dims(X, InitCond)) .+ 1

# generate initial conditions

function initcondgener(A::StateSpaceSet, t)
    cols = columns(A)
    L = length(A)
    Ls = round(Int, 0.01L)
    u0s = map(c -> mean(c[1:Ls]), cols)
    return SVector(u0s)
end

u0s = StateSpaceSet(map(ics) do ic
    A = gnv(X[InitCond(ic)])
    A = StateSpaceSet(A)
    f = initcondgener(A, nothing)
end)

# %% Cluster
using Statistics
function featurizer(A::StateSpaceSet, t)
    cols = columns(A)
    L = length(A)
    Ls = round(Int, 0.75L)
    stds = map(c -> std(c[Ls:L]), cols)
    # We included the mean but no mean was chosen as optimal feature!
    # means = map(c -> mean(c[Ls:L]), cols)
    # return SVector(stds..., means...)
    return SVector(stds)
end

function extract_features(X)
    ics = gnv(dims(X, InitCond)) .+ 1
    allfeatures = map(ics) do ic
        A = gnv(X[InitCond(ic)]) # careful: assumes time axis is first axis!
        A = StateSpaceSet(A)
        f = featurizer(A, nothing)
    end
    allfeatures = rescale_to_01(allfeatures) # this makes it an SSSet as well
    return allfeatures
end

allfeatures = extract_features(X)

ca = ADBSCAN(; min_neighbors = 10)
fcq = FeaturesClusteringQuality(; attractor_weight = 1.0, feature_weight = 0.5, ca)

best_choice, best_labels = optimize_feature_selection(fcq, allfeatures;
    verbose = false, minF = 2, reclustering = 0,
)

best_choice
attractors = unique(best_labels)
natt = length(attractors)

feature_names = diagnostics


# %%
