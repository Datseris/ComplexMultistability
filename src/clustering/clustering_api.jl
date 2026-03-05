using ClusteringAPI
using ClusteringAPI: cluster

# These I am not sure yet whether they work:
# include("clustering_quickshift.jl")
# include("clustering_hclust_auto.jl")

# two helper functions for agnostic input data type
"""
    input_data_size(data) → (d, m)

Return the data point dimension and number of data points.
"""
input_data_size(A::AbstractMatrix) = size(A)
input_data_size(A::AbstractVector{<:AbstractVector}) = (length(first(A)), length(A))

"""
    each_data_point(data)

Return an indexable iterator over each data point in `data`, that can be
indexed with indices `1:m`.
"""
each_data_point(A::AbstractMatrix) = eachcol(A)
each_data_point(A::AbstractVector{<:AbstractVector}) = A

include("clustering_adbscan.jl")
