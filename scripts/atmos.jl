using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
include("atmos_shared.jl")

# %% What we need

best_labels
best_choice
diagnostics
allfeatures
u0s
natt
featcols = columns(allfeatures)
u0cols = columns(u0s)

# from which we can generate the rest.
imatrix_features = intermingledness(allfeatures, best_labels)
imatrix_basins = intermingledness(u0s, best_labels)
choices = [best_choice[1:2], [23, 32], [21, length(diagnostics)]] # two plots from each

# We have too many diagnostics so we pick less
reduced = vcat(choices...)
append!(reduced, [1, 2, 3, 5, 10, 7, 8, 12, 18, 30, 24, 27, 34, 14, 38])
sort!(unique!(reduced))

imatrix_features = imatrix_features[:, reduced]
imatrix_basins = imatrix_basins[:, reduced]

COLORDict = Dict{Int, Any}(i => COLORS[i] for i in 1:12)
COLORDict[-1] = :red
MARKERDict = Dict{Int, Any}(i => MARKERS[i] for i in 1:12)
MARKERDict[-1] = 'F'

# %%

fig = Figure(size = (figwidth, 1.5figheight))
glleft = GridLayout(fig[:, 1])
figfeatures = GridLayout(glleft[1,1])
title1 = "features (std of last 4th)"
title2 = "basins (mean at first 10th)"
figuretitle!(figfeatures, "a: $(title1)"; halign = :left)
figbasins = GridLayout(glleft[2,1])
figuretitle!(figbasins, "b: $(title2)"; halign = :left)
figinter = GridLayout(fig[:,2][1,1]; alignmode = Outside())
figuretitle!(figinter, "c: intermingledness"; halign = :left)

# Features & basins plot (same style)
axs = axesgrid!(figfeatures, 1, 3; sharex = false)
hidedecorations!.(axs; label = false, grid = true )
axsba = axesgrid!(figbasins, 1, 3; sharex = false)
hidedecorations!.(axsba; label = false, grid = true )

for (k, (i, j)) in enumerate(choices)
    if k == 1
        labels = best_labels
    else
        cr = cluster(ca, vec(allfeatures[:, [i, j]]))
        labels = cr.labels
    end
    for (axs, cols) in zip((axs, axsba), (featcols, u0cols))
        ax = axs[k]
        if cols === u0cols
            j = mod1(j, length(diagnostics)) # for when we use both mean and std as features
        end
        # cluster, get labels, etc.
        scatter!(ax, cols[i], cols[j];
            marker = getindex.(Ref(MARKERDict), labels),
            color = getindex.(Ref(COLORDict), labels),
            markersize = 15,
            strokewidth = 0.25, strokecolor = "black",
        )
        ax.xlabel = feature_names[i]
        ax.ylabel = feature_names[j]
    end
end
text!(axs[1], 0.05, 0.4; text = "best\nfeatures", space = :relative, )
text!(axs[2], 0.5, 0.1; text = "expert\nchosen", space = :relative, )
text!(axs[3], 0.05, 0.3; text = "random\nfeatures", space = :relative, )

# # intermingledness

plot_twice_intermingledness!(figinter, imatrix_features, imatrix_basins,;
    title1, title2,
    names = feature_names[reduced], unitlen = 20, fixed_size = true,
)

# add timeseries
gl = GridLayout(glleft[3, 1])
ax = Axis(gl[1,1]; xlabel = "time (days)")
u0_idxs = [findfirst(isequal(k), best_labels) for k in 1:natt]
dim = best_choice[1]
for (j, u0idx) in enumerate(u0_idxs)
    x = gnv(X[Diagnostic(dim), InitCond(u0idx)])
    lines!(ax, t, x; color = Cycled(j), linewidth = 1, linestyle = :solid,)
end
ax.ylabel = diagnostics[dim]
hideydecorations!(ax; label = false)
xlims!(ax, 0, t[end])

Legend(gl[1,2],
    [MarkerElement(color = COLORS[l], marker = MARKERS[l], markersize = 15) for l in 1:natt],
    string.(1:natt), "group ID"; position = :lt, nbanks = 1, tellheight = false, tellwidth = true,
)

# colsize!(fig.layout, 1, Relative(0.5))
figinter.valign = :top
colgap!(fig.layout, 30)
# colsize!(fig.layout, 1, Relative(0.5))
rowgap!(figfeatures, 5)
colgap!(figfeatures, 10)
rowgap!(figbasins, 5)
colgap!(figbasins, 10)
resize_to_layout!(fig)


fig

# %%
wsave(papersdir("figures", "atmos_analysis"), fig)
