# This script analyses intermingledness, boundary intermingledness, and basin entropy
# It uses the definitions from Attractors.jl
using DrWatson
@quickactivate "ComplexMultistability"
using PredefinedDynamicalSystems, Attractors
include(srcdir("papertheme.jl"))
# we don't include the intermingledness src file, the function in Attractors.jl is better
N = 100
ds = PredefinedDynamicalSystems.magnetic_pendulum(d=0.2, α=0.1, ω=0.8, N=3)
psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
xg = yg = range(-4, 4; length = N)
grid = (xg, yg)
mapper = AttractorsViaRecurrences(psys, grid; Δt = 1)
basins, attractors = basins_of_attraction(mapper, grid; show_progress = false)

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

# also prepare veros basin data
input_file_name = "Veros_Data_F8_1.nc"
include("veros_shared_code.jl")
# also need initial conditions
allu0s = zeros(length(ics), length(diagnostics))
for i in eachindex(ics)
    slice = X[InitCond(i)]
    # To get actual initial condition is not trivial due to different evolution time
    # so we find the first non-missing entry
    x = gnv(slice[Diagnostic(1)])
    j = findfirst(!ismissing, x)
    allu0s[i, :] .= gnv(slice[Time(j)])
end

veros_basins = best_labels
veros_points = StateSpaceSet(allu0s[:, best_choice[[1, 2]]])
veros_points = rescale_to_01(veros_points)

# %%

# finally, make a basin with a smaller basin inside
basins4 = ones(Int, N, N)
L = N÷10
for i in (N÷2)-L:(N÷2)+L
    for j in (N÷2)-L:(N÷2)+L
        basins4[i, j] = 2
    end
end

fig = Figure(size = (figwidth, 2figheight))

using Statistics: mean
for (i, b) in enumerate((basins, basins2, basins3))
    ax = Axis(fig[1,i])

    if length(unique(b)) > 3
        replace!(b, 4 => 3) # it's because it's sticking to the saddle point in the middle
    end

    heatmap!(ax, xg, yg, b; colormap = to_color.([COLORS[1], COLORS[2], COLORS[3]]))
    be = Attractors.basin_entropy(b, L)[1]
    inter = Attractors.intermingledness(points, vec(b)) # points and labels must have same layout
    inter2 = Attractors.boundary_intermingledness(points, vec(b)) # points and labels must have same layout
    ax.title = "intermingledness = $(round(mean(values(inter)); digits = 2))\nb-intermingledness = $(round(mean(values(inter2)); digits = 2))\nbasin entropy = $(round(be; digits = 2))"
    hidedecorations!(ax)

    if i == 1
        r = 0.5
        center = Point2f(-2.5, 2.5)
        pl = poly!(ax, Circle(center, r), color = (:red, 0.5), strokecolor = :red, strokewidth = 2.5)
        translate!(pl, 0, 0, +999)
    end

end

ax4, = heatmap(fig[2,1], basins4; colormap = to_color.([COLORS[1], COLORS[2]]))
inter4 = Attractors.intermingledness(points, vec(basins4)) # points and labels must have same layout
inter4new = Attractors.boundary_intermingledness(points, vec(basins4)) # points and labels must have same layout
ax4.title = "intermingledness (purple) = $(round(maximum(values(inter4)); digits = 2))\n b-intermingledness (purple) = $(round(maximum(values(inter4new)); digits = 2))"
hidedecorations!(ax4)

basins6 = ones(Int, N, N)
for i in 1:(N÷5):N
    basins6[i:(i + N÷10), :] .= 2
end
ax6, = heatmap(fig[2,2], basins6; colormap = to_color.([COLORS[1], COLORS[2]]))
inter6 = Attractors.intermingledness(points, vec(basins6))
inter6new = Attractors.boundary_intermingledness(points, vec(basins6))
ax6.title = "intermingledness = $(round(mean(values(inter6)); digits = 2))\nb-intermingledness = $(round(mean(values(inter6new)); digits = 2))"
hidedecorations!(ax6)

# do veros basins here  something better here:
ax5, = scatter(fig[2,3], veros_points; color = [COLORS[l] for l in best_labels])
inter5 = Attractors.intermingledness(veros_points, veros_basins)
inter5new = Attractors.boundary_intermingledness(veros_points, veros_basins)
ax5.title = "intermingledness (yellow) = $(round(maximum(values(inter5)); digits = 2))\nb-intermingledness (yellow) = $(round(maximum(values(inter5new)); digits = 2))"
hidedecorations!(ax5)

fig

# %%
wsave(papersdir("figures", "intermingledness"), fig)
negate_remove_bg(papersdir("figures", "intermingledness.png"))