# These are the default values
ENV["COLORSCHEME"] = "JuliaDynamicsLight" # or others, see docs
ENV["BGCOLOR"] = :transparent        # anything for `backgroundcolor` of Makie
ENV["AXISCOLOR"] = :white            # color of all axis elements (labels, spines, ticks)

using MakieForProjects, CairoMakie # this has now set the theme already!

# you may further edit the set theme by using
Makie.update_theme!(;
    # size = (figwidth, figheight),
)


cmap = to_color.([:blue, :white, :red])
intermingledness_heatmap_kw = (
    hmapkw = (colorrange = (0.0,1.0), colormap = cgrad(cmap, 9, categorical = true))
)
