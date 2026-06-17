# This script compares analyses intermingledness
# the first component of the script compares it with basin entropy
using DrWatson
@quickactivate "ComplexMultistability"
using PredefinedDynamicalSystems, Attractors
include(srcdir("papertheme.jl"))
# we don't include the intermingledness src file, the function in Attractors.jl is better

ds = PredefinedDynamicalSystems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
xg = yg = range(-4, 4; length = 200)
grid = (xg, yg)
mapper = AttractorsViaRecurrences(psys, grid; Δt = 1)
basins, attractors = basins_of_attraction(mapper, grid; show_progress = false)

set_parameter!(psys, :γs, [1.0, 1.0, 1.0])
set_parameter!(psys, :α, 1.0)
reset_mapper!(mapper)
# mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 1)
basins2, attractors2 = basins_of_attraction(mapper, grid; show_progress = false)

set_parameter!(psys, :α, 2.0)
set_parameter!(psys, :d, 0.4)
reset_mapper!(mapper)
# mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 1)
basins3, attractors3 = basins_of_attraction(mapper, grid; show_progress = false)

points = ics_from_grid(grid)

# %%

# finally, make a basin with a smaller basin inside
basins4 = ones(Int, 200, 200)
L = 20
for i in 100-L:100+L
    for j in 100-L:100+L
        basins4[i, j] = 2
    end
end
basins5 = rand(1:2, 200, 200)

fig = Figure(size = (figwidth, 2figheight))

using Statistics: mean
for (i, b) in enumerate((basins, basins2, basins3))
    ax = Axis(fig[1,i])

    if length(unique(b)) > 3
        replace!(b, 4 => 3) # it's because it's sticking to the saddle point in the middle
    end

    heatmap!(ax, xg, yg, b; colormap = to_color.([COLORS[1], COLORS[2], COLORS[3]]))
    inter = Attractors.intermingledness(points, vec(b)) # points and labels must have same layout
    v = mean(values(inter))
    be = Attractors.basin_entropy(b, 5)[1]
    ax.title = "intermingledness = $(round(v; digits = 1))\nbasin entropy = $(round(be; digits = 1))"
    hidedecorations!(ax)

    if i == 1
        r = 0.5
        center = Point2f(-2.5, 2.5)
        pl = poly!(ax, Circle(center, r), color = (:red, 0.5), strokecolor = :red, strokewidth = 2.5)
        translate!(pl, 0, 0, +999)
    end

end

inter4 = Attractors.intermingledness(points, vec(basins4)) # points and labels must have same layout
ax4, = heatmap(fig[2,1], basins4; colormap = to_color.([COLORS[1], COLORS[2]]))
ax4.title = "intermingledness (purple) = $(round(inter4[1]; digits = 1))\nintermingledness (black) = $(round(inter4[2]; digits = 1))"
hidedecorations!(ax4)

inter5 = Attractors.intermingledness(points, vec(basins5))
ax5, = heatmap(fig[2,2], basins5; colormap = to_color.([COLORS[1], COLORS[2]]))
ax5.title = "intermingledness = $(round(mean(values(inter5)); digits = 1))"
hidedecorations!(ax5)

# finally, completely separated basins without touching
randompoint() = [rand([-1, 1]) + 0.1randn(), 0.1randn()]
points2 = [randompoint() for _ in 1:1000]
labels6 = [x[1] < 0 ? 1 : 2 for x in points2]
inter6 = Attractors.intermingledness(points2, labels6)
ax6, = scatter(fig[2,3], StateSpaceSet(points2); color = labels6, colormap = [COLORS[1], COLORS[2]])
ax6.title = "intermingledness = $(round(mean(values(inter6)); digits = 1))"

fig

# %%
wsave(papersdir("figures", "compare_with_be"), fig)
negate_remove_bg(papersdir("figures", "compare_with_be.png"))