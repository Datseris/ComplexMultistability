# In this script we will run the Veros data analysis, using only at most 2 features
# (for performance). Then we will vary the two key weights w_A, w_F to see how the
# sensitivity of the algorihm behaves.

using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
# Change theme to paper
include(srcdir("papertheme.jl"))
input_file_name = "Veros_Data_F8_1.nc"
include("veros_shared_code.jl")

weights = 0.1:0.3:1.0

result = map(Iterators.product(weights, weights)) do (w_A, w_F)
    fcq = FeaturesClusteringQuality(; attractor_weight = w_A, feature_weight = w_F, ca)
    best_choice, best_labels = optimize_feature_selection(fcq, allfeatures;
        verbose = false, max_choices_per_dim = typemax(Int), minF = 1, maxF = 3,
    )
    s = "$(sort(best_choice)): $(length(unique(best_labels)))"
end

# %%

fig = Figure(size = (figwidth÷2, figheight))
ax = Axis(fig[1,1]; xlabel = L"w_A", ylabel = L"w_F")

for (j, (w_A, w_F)) in enumerate(Iterators.product(weights, weights))
    text!(ax, w_A, w_F; text = string(result[j]), align = (:center, :center))
end

ax.xgridvisible = ax.ygridvisible = false
ax.limits = ((0, 1.1), (0, 1.1))
ax.xticks = weights
ax.yticks = weights
fig

# %%

wsave(papersdir("figures", "sensitivity_of_weights"), fig)