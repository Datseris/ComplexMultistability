# ComplexMultistability

This reproducible code base accompanies the research article:
> Multistability and intermingledness in complex high-dimensional data
by Datseris et al., 2026.
The code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ComplexMultistability

# Applying this to your own data

To apply the same analysis of multistability to intermingledness as in the paper, you need to cast your data into diagnostic variables (Step 2. of the main figure of the paper).
The format must be a NetCDF file, with three dimensions (named exactly like this): `time`, `diagnostic`, `ic`. The data themselves should be a single field that covers all three aforementioned dimensions. Typical name is `values`, although you can easily alter this in the scripts. Once you have the NetCDF file, follow the instructions in the `notebooks/workflow.ipynb` file, which is a Jupyter notebook providing an exemplary application of the workflow.

# Paper reproducibility

The remainder of the files in the `scripts` folder reproduce the figures of the paper.
To achieve this you need to first download the model data used (available on request!).

Then, to reproduce the project, do the following:

0. Download this code base.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ComplexMultistability"
```
which auto-activate the project and enable local path handling from DrWatson.
