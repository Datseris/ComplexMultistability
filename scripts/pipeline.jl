using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
# Change theme to paper
include(srcdir("papertheme.jl"))
input_file_name = "Veros_Data_F8_1.nc"
include("veros_shared_code.jl")

# With this in mind, we can now plot everything.
# %%
fig = Figure(size = (1.33figwidth, 2figheight))
top = GridLayout(fig[1,1])
bottom = GridLayout(fig[2,1])

figuretitle!(top[1,1], "1. complex multi-state system\n(many different simulations)")
ax0 = Axis(top[1,1][1,1])
hidedecorations!(ax0)

# Second plot: diagnostics timeseries
figuretitle!(top[1,2], "2. diagnostic variables timeseries\n(collapse space dimensions)")
chosen_ics = [24, 1, 90]
chosen_diagnostics = [5, 1]
axs = axesgrid!(top[1,2], length(chosen_diagnostics), 1; xlabels = "time (years)", sharex = true, ylabels = diagnostics[chosen_diagnostics])

j = 3401 # first nonmissing entries
timevec = ttotal[j:end] .- ttotal[j]
for ic in chosen_ics
    slice = X[InitCond(ic)]
    for (d, diag) in enumerate(chosen_diagnostics)
        x = gnv(slice[Diagnostic(diag)])
        id = best_labels[ic]
        lines!(axs[d], timevec, x[j:end]; color = Cycled(id), linestyle = Cycled(id), label = "")
    end
end
axislegend(axs[end], "initial condition"; backgroundcolor = (:white, 0.8), halign = :right, valign = 0.3, nbanks = 3)
xlims!(axs[end], 0, 3000)

figuretitle!(top[1,3], "3. cast into feature space\n(collapse time dimension)")
axf = Axis(top[1,3][1,1]; xlabel = "mean-last-1000-yr: salt_tot", ylabel = "std-total: temp_sub_SA")
axf.spinewidth = 3

using Statistics: mean, std
function featurizer(A, t)
    m = mean(A[end-200:end, chosen_diagnostics[1]])
    s = std(A[:, chosen_diagnostics[2]])
    return SVector(m, s)
end

features = cast_to_features(Xvalid, featurizer)

scatter!(axf, features; color = :gray, markersize = 10, label = "all i.c.")
for ic in chosen_ics
    id = best_labels[ic]
    scatter!(axf, features[ic]; markersize = 30,
        marker = Cycled(id), color = (COLORS[id], 0.75), strokecolor = COLORS[id],
        strokewidth = 3,
    )
end
axislegend(axf)
hidedecorations!(axf; label = false)

figuretitle!(top[1,4], "4. group features\ninto different groups")
axff = Axis(top[1,4][1,1]; xlabel = "mean-last-1000-yr: salt_tot", ylabel = "std-total: temp_sub_SA")
axff.spinewidth = 3

ca = ADBSCAN(; rescale_features = true, min_neighbors = 4, optimal_radius_method = "silhouettes")
fcq = FeaturesClusteringQuality(; attractor_weight = 5.0, feature_weight = 0.1, ca)
odd_choice, odd_labels = optimize_feature_selection(fcq, features;
    verbose = false, max_choices_per_dim = typemax(Int)
)
scatter!(axff, features;
    alpha = 0.75,
    color = [COLORS[l] for l in odd_labels], marker = [MARKERS[l] for l in odd_labels]
)
hidedecorations!(axff; label = false)

figuretitle!(bottom[1,4], "5. optimize choise of features\n(which features distinguish the most)")

figuretitle!(bottom[1,3], "6. unique groups identified")
axa = Axis3(bottom[1,3][1,1])

i, j, k = best_choice
axa.xlabel = "mean:\n"*diagnostics[i]
axa.ylabel = "mean:\n"*diagnostics[j]
axa.zlabel = "mean: "*diagnostics[k]

scatter!(axa, allfeatures[:, best_choice];
    alpha = 0.75,
    color = [COLORS[l] for l in best_labels], marker = [MARKERS[l] for l in best_labels]
)
hidedecorations!(axa; label = false)
axa.protrusions = (0, 0, 0, 0)
axa.viewmode = :fitzoom
axa.width = Relative(1)
axa.azimuth = 5.05530633326986
axa.elevation = 0.3
axa.xlabeloffset = 0
axa.zlabeloffset = 0
axa.ylabeloffset = 0
axa.xspinewidth = 3
axa.yspinewidth = 3
axa.zspinewidth = 3

figuretitle!(bottom[1,2], "7. approximate basins\n(if applicable)")
axb = Axis(bottom[1,2][1,1])

allu0s = zeros(length(ics), length(diagnostics))
for i in eachindex(ics)
    slice = X[InitCond(i)]
    # To get actual initial condition is not trivial due to different evolution time
    # so we find the first non-missing entry
    x = gnv(slice[Diagnostic(1)])
    j = findfirst(!ismissing, x)
    allu0s[i, :] .= gnv(slice[Time(j)])
end
i, j = best_choice
scatter!(axb, allu0s[:, i], allu0s[:, j];
    alpha = 0.75,
    color = [COLORS[l] for l in best_labels], marker = [MARKERS[l] for l in best_labels]
)
axb.xlabel = diagnostics[i]*"(t = 0)"
axb.ylabel = diagnostics[j]*"(t = 0)"

axislegend(axb,
    [MarkerElement(color = COLORS[l], marker = MARKERS[l], markersize = 15) for l in 1:5],
    string.(1:5), "group ID"; position = :lt,
)

figuretitle!(bottom[1,1], "8. intermingledness")

imatrix = intermingledness(StateSpaceSet(allu0s), best_labels)
imatrix = imatrix[:, best_choice]

cmap = to_color.([:white, :black])
hmapkw = (colorrange = (0.0,1.), colormap = cgrad(cmap, 9, categorical = true), highclip = :black,)
axbgkw = (xgridvisible = false, ygridvisible = false, xminorgridvisible = true,
    yminorgridvisible = true, xminorgridcolor = :grey, yminorgridcolor = :grey,
    xminorticks = IntervalsBetween(2), yminorticks = IntervalsBetween(2),
)
axhm = Axis(bottom[1,1][1,1]; ylabel = "group ID", xlabel = "feature dimension", axbgkw...)
hm = heatmap!(axhm, imatrix'; hmapkw...)
axhm.xticks = 1:3
translate!(hm, 0, 0, -100)
cb = Colorbar(bottom[1,1][1,2]; width = 20, hmapkw...)
cb.ticks = ([0.0, 1.0], ["0", "1"])
# cb.ticklabelrotation = -π/2
# cb.labelrotation = -π/2
cb.labelpadding = -20

colsize!(bottom, 4, Relative(0.4))

colgap!(top, 30)

display(fig)

# %%
wsave(papersdir("figures", "pipeline_original"), fig)
