# optimization using JuMP
# I have asked about this on discourse in two places:
# https://discourse.julialang.org/t/how-to-prioritize-choosing-some-parameters-over-others-in-hyperopt-jl-or-similar-packages/102594
# https://discourse.julialang.org/t/using-a-jump-bin-vector-to-index-an-array-checkbounds-fails-with-type-variableref/102674
# but the answers there were disheartening/unhelpful eluding that my problem isn't solvable with JuMP.
# The problem also isn't solvable with black box optimizers, as the input variables are
# binary, not real.

# Now I found out: https://github.com/hdavid16/DisjunctiveProgramming.jl
# Then I asked the developer, who confirmed that it is not possible!
