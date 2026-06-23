# These are the default values
ENV["COLORSCHEME"] = "JuliaDynamicsLight" # or others, see docs
ENV["BGCOLOR"] = :transparent        # anything for `backgroundcolor` of Makie
ENV["AXISCOLOR"] = :white            # color of all axis elements (labels, spines, ticks)

using MakieForProjects, CairoMakie # this has now set the theme already!

# you may further edit the set theme by using
Makie.update_theme!(;
    # size = (figwidth, figheight),
)
