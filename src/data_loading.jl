using StateSpaceSets

const Diagnostic = Dim{:diagnostic}
const InitCond = Dim{:ic}

function select_in_common_time(X)
    if !any(ismissing, X)
        return X
    end
    t = gnv(dims(X, Ti))
    tvalid = t[1]
    Xred = X[Diagnostic(1)] # reduce the amount of data necessary to process
    for i in dims(X, InitCond)
        tvalid_index = findfirst(!ismissing, Xred[InitCond(At(i))])
        if t[tvalid_index] > tvalid
            tvalid = t[tvalid_index]
        end
    end
    Xvalid = X[Ti(Between(tvalid, t[end]))]
    return float.(Xvalid)
end

function cast_to_features(X::ClimArray, featurizer)
    # ensure X has dimensions ordered in the way we want
    X = permutedims(X, (Ti, Diagnostic, InitCond))
    t = gnv(dims(X, Ti))
    allfeatures = map(eachindex(dims(X, InitCond))) do i
        A = gnv(X[InitCond(i)])
        A = StateSpaceSet(A)
        f = featurizer(A, t)
    end
    return StateSpaceSet(allfeatures)
end

function initial_conditions(X::ClimArray)
    ics = gnv(dims(X, InitCond))
    XIC = ClimArray(zeros(length.((ics, diagnostics))), dims(X, InitCond, Diagnostic))
    for i in eachindex(ics)
        slice = X[InitCond(i)]
        # To get actual initial condition is not trivial due to different evolution time
        # so we find the first non-missing entry
        x = gnv(slice[Diagnostic(1)]) # the diagnostic doesn't matter here
        j = findfirst(!ismissing, x)
        XIC[InitCond(i)] .= gnv(slice[Time(j)])
    end
    return XIC
end


rescale_to_01(features::Vector{<:AbstractVector}) = rescale_to_01(StateSpaceSet(features))
function rescale_to_01(features::AbstractStateSpaceSet)
    mini, maxi = minmaxima(features)
    rescaled = map(f -> (f .- mini) ./ (maxi .- mini), features)
    return StateSpaceSet(rescaled)
end
