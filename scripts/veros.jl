using DrWatson
@quickactivate "ComplexMultistability"
using ClimateBase, StateSpaceSets
# Change theme to paper
include(srcdir("papertheme.jl"))
input_file_name = "Veros_Data_F8_1.nc"
include("veros_shared_code.jl")

# %% This is what I need for the figures
using Combinatorics
best_labels
best_choice
diagnostics
allfeatures
natt
u0s = StateSpaceSet(Xvalid[Tim(1)])
featcols = columns(allfeatures)
u0cols = columns(u0s)
# colormap = cgrad(COLORS[1:5], 5; categorical = true)
colormap = Makie.Categorical(COLORS[1:5])

# from which we can generate the rest.
imatrix_features = intermingledness(allfeatures, best_labels)
imatrix_basins = intermingledness(u0s, best_labels)
choices = [best_choice[1:2], [3, 11], [8, length(diagnostics)]] # two plots from each

# %%

fig = Figure(size = (figwidth, 1.3figheight))
glleft = GridLayout(fig[:, 1])
figfeatures = GridLayout(glleft[1,1])
figuretitle!(figfeatures, "a: features (last 1000yr mean)"; halign = :left)
figbasins = GridLayout(glleft[2,1])
figuretitle!(figbasins, "b: basins (value at time = 0)"; halign = :left)
Legend(glleft[3,1],
    [MarkerElement(color = COLORS[l], marker = MARKERS[l], markersize = 15) for l in 1:5],
    string.(1:5), "group ID"; position = :lt, nbanks = 5, tellheight = true, tellwidth = false,
)
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
        replace!(labels, -1 => 0) # plot differently bad solutions
    end
    for (axs, cols) in zip((axs, axsba), (featcols, u0cols))
        ax = axs[k]
        # cluster, get labels, etc.
        scatter!(ax, cols[i], cols[j];
            color = labels, alpha = 0.75,  colormap,  markersize = 10,
            strokewidth = 0.25, strokecolor = "black",
            marker = getindex.(Ref(MARKERS), labels),

        )
        ax.xlabel = diagnostics[i]
        ax.ylabel = diagnostics[j]
    end
end
text!(axs[1], 0.1, 0.3; text = "best\nfeatures", space = :relative, )
text!(axs[2], 0.1, 0.3; text = "random\nfeatures", space = :relative, )
text!(axs[3], 0.1, 0.3; text = "worse\nfeatures", space = :relative, )

# intermingledness
plot_twice_intermingledness!(figinter, imatrix_features, imatrix_basins,;
    title1 = "features (last 1000yr mean)", title2 = "basins (value at time 0)",
    names = diagnostics, unitlen = 20, fixed_size = true,
)

# # Layouting
figinter.valign = :top
colgap!(fig.layout, 30)
colsize!(fig.layout, 1, Relative(0.5))
rowgap!(figfeatures, 5)
colgap!(figfeatures, 10)
rowgap!(figbasins, 5)
colgap!(figbasins, 10)
resize_to_layout!(fig)

fig

# %%

wsave(papersdir("figures", "veros_analysis"), fig)
