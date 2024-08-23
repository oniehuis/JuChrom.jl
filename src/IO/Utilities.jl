function buildxic(pointcounts, ionvec, intsvec)
    scancount = length(pointcounts)
    ions = sort(unique(ionvec))
    ionidx = Dict(ion => i for (i, ion) in enumerate(ions))
    xic = zeros(eltype(intsvec), scancount, length(ions))
    start = 1
    for iₛ in 1:scancount
        stop = start + pointcounts[iₛ] - 1
        for (i₂, ion) in enumerate(@view ionvec[start:stop])
            xic[iₛ, ionidx[ion]] += (@view intsvec[start:stop])[i₂]
        end
        start += pointcounts[iₛ]
    end
    ions, xic
end

abstract type Source end
struct File <: Source end
struct Path <: Source end

struct IOError <: Exception
    var::String
end

Base.showerror(io::IO, e::IOError) = print(io, e.var)

struct FileExistsError <: Exception
    var::String
end

Base.showerror(io::IO, e::FileExistsError) = print(io, e.var)
