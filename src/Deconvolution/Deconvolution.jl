module Deconvolution

include("../NNLS/NNLS.jl")
import NNLS: nnls
export nnls

export ls_fit

function ls_fit(xs, ys)
    if length(xs) ≠ length(ys)
        throw(ArgumentError("vectors differ in size"))
    end
    if length(xs) ≤ 1 
        throw(ArgumentError("each vector must consist of at least two values"))
    end
    x̅, y̅ = (sum(xs), sum(ys)) ./ length(xs)
    numerator = zero(eltype(xs))
    denominator = zero(eltype(xs)) * zero(eltype(xs))
    for i in eachindex(xs)
        numerator += (xs[i] - x̅) * (ys[i] - y̅)
        denominator += (xs[i] - x̅)^2
    end
    slope = numerator / denominator
    intercept = y̅ - slope * x̅
    intercept, slope
end

end  # module