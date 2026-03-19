using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
include(srcdir("data_loading.jl"))
include(srcdir("clustering", "clustering_api.jl"))
include(srcdir("optimization", "optimization_api.jl"))
include(srcdir("intermingledness.jl"))

# Change theme to paper
include(srcdir("papertheme.jl"))

# load and process data (see notebooks for step by step)
file = datadir("samosa", "exoplasim_samosa460.nc")
X = ncread(file, "values")
ics = gnv(dims(X, InitCond)) .+ 1 # make 1-based
diagnostics = gnv(dims(X, Diagnostic))
correct_idxs = 1:12
habit_idxs = setdiff(eachindex(diagnostics), correct_idxs)
allfeatures_original = StateSpaceSet(X[:, correct_idxs])
diagnostics = diagnostics[correct_idxs]
allfeatures = StateSpaceSets.standardize(allfeatures_original)
habit_features = StateSpaceSet(X[:, habit_idxs])
pressure = ncread(file, "Pressure")
instellation = ncread(file, "Instellation")
parameters = StateSpaceSet(hcat(vec(log.(pressure)), vec(instellation)))

# create "basins of attraction" (habitability composite)
attractor_id(habit1, habit2) = (Bool(habit1) << 1) + Bool(habit2) + 1
best_labels = @. attractor_id(habit_features[:, 1], habit_features[:, 4]) .- 1

# %%
# Make overarching figure now

fig = Figure(size = (figwidth, 1.5figheight))

figinter = GridLayout(fig[1,2])
figbasins = GridLayout(fig[1,1])
figdiagnostics = GridLayout(fig[2, :])

# plot basins of parameters
ids = sort!(unique(best_labels))
axbasins = Axis(figbasins[1, 1], title = "a: habitability versus parameters", xlabel =  "log(pressure)", ylabel = "instellation")
colormap = cgrad(COLORS[ids], length(ids); categorical = true)
hmap = heatmap!(axbasins, columns(parameters)..., best_labels; colormap, colorrange = (0.5, 3.5))
cb = Colorbar(figbasins[1,2], hmap; ticks = (Int[1,2,3], string.([1,2,3])))

# calculate and plot intermingledness
imatrix = intermingledness(allfeatures, best_labels; summarizer = mean)
plot_intermingledness!(figinter, imatrix;
    names = diagnostics
)
figuretitle!(figinter, "b: intermingledness of diagnostics"; halign = :left)

# plot basins of diagnostics
pairs = [(6, 12), (1, 4), (5, 11)]
attractor_color = [COLORS[j] for j in best_labels]
attractor_marker = [MARKERS[j] for j in best_labels]
cols = columns(allfeatures_original)
axs = axesgrid!(figdiagnostics, 1, length(pairs))

for (k, ax) in enumerate(axs)
    ax = axs[k]
    i, j = pairs[k]
    scatter!(ax, cols[i], cols[j]; alpha = 0.75, color = attractor_color, marker = attractor_marker)
    ax.xlabel = diagnostics[i]
    ax.ylabel = diagnostics[j]
end
figuretitle!(figdiagnostics, "c: diagnostics \"basins\""; halign = :left)
# and add a legend
Legend(figdiagnostics[1,2],
    [MarkerElement(color = COLORS[i], marker = MARKERS[i], markersize = 20) for i in ids],
    string.(1:3), "ID"; halign = :left, valign = :bottom, nbanks = 1, tellwidth = false,
    margin = (6, 6, 6, 6), backgroundcolor = (:white, 0.9), patchlabelgap = -5,
)

# make sure there is enough space
rowsize!(fig.layout, 1, Relative(0.5))
colgap!(fig.layout, 30)
rowgap!(fig.layout, 00)
rowgap!(figdiagnostics, 5)
resize_to_layout!(fig)

display(fig)

# %% save
wsave(papersdir("figures", "samosa_analysis"), fig)
