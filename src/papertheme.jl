# These are the default values
ENV["COLORSCHEME"] = "JuliaDynamics" # or others, see docs
ENV["BGCOLOR"] = :white        # anything for `backgroundcolor` of Makie
ENV["AXISCOLOR"] = :black            # color of all axis elements (labels, spines, ticks)

using MakieForProjects, CairoMakie # this has now set the theme already!

theme = MakieForProjects.make_theme()
Makie.set_theme!(theme)

# you may further edit the set theme by using
figwidth, figheight = (1200, 400)
Makie.update_theme!(;
    size = (figwidth, figheight),
)
