"""
    JuChrom.InputOutput.buildxic(pointcounts::AbstractVector{<:Integer}, 
    ionvec::AbstractVector{<:Real}, intsvec::AbstractVector{<:Real})

Return a sorted list of ions and an intensity matrix with the number of columns equal to 
the number of ions in the sorted list and the number of rows equal to the number of scans. 
This is based on a vector that stores the number of ion-intensity pairs per scan 
(`pointcounts`), a vector of continuously recorded ions throughout the run (`ionvec`), and 
a vector of intensity values associated with those ions (`intsvec`).

# Examples
```jldoctest
julia> pointcounts = [2, 3, 2];

julia> ionvec = [85.1, 100.2, 85.2, 99.9, 112.1, 84.9, 100.6];

julia> intsvec = Int64[12, 234, 23, 324, 45422, 21, 523];

julia> mzs, xic = JuChrom.InputOutput.buildxic(pointcounts, ionvec, intsvec);

julia> mzs
7-element Vector{Float64}:
  84.9
  85.1
  85.2
  99.9
 100.2
 100.6
 112.1

julia> xic
3×7 Matrix{Int64}:
  0  12   0    0  234    0      0
  0   0  23  324    0    0  45422
 21   0   0    0    0  523      0
```
"""
function buildxic(pointcounts::AbstractVector{<:Integer}, ionvec::AbstractVector{<:Real}, 
    intsvec::AbstractVector{<:Real})
    length(ionvec) == length(intsvec) || throw(
        ArgumentError("ion and intensity values differ in their counts"))
    sum(pointcounts) == length(ionvec) == length(intsvec) || throw(
        ArgumentError("point counts incompatible with ion and intensity value counts"))
    ions = sort(unique(ionvec))
    ionidx = Dict(ion => i for (i, ion) in enumerate(ions))
    scancount = length(pointcounts)
    xic = zeros(eltype(intsvec), scancount, length(ions))
    start = 1
    for iₛ in 1:scancount
        stop = start + pointcounts[iₛ] - 1
        for (iᵢ, ion) in enumerate(@view ionvec[start:stop])
            xic[iₛ, ionidx[ion]] += (@view intsvec[start:stop])[iᵢ]
        end
        start += pointcounts[iₛ]
    end
    ions, xic
end


"""
    JuChrom.InputOutput.Source

Supertype for all data sources.

See also [`File`](@ref), [`Path`](@ref).
"""
abstract type Source end


"""
    JuChrom.InputOutput.File() <: Source

Type that indicates a source is a file.

See also [`Source`](@ref), [`Path`](@ref).
"""
struct File <: Source end


"""
    JuChrom.InputOutput.Path() <: Source

Type that indicates a source is a path.

See also [`Source`](@ref), [`File`](@ref).
"""
struct Path <: Source end


"""
    JuChrom.InputOutput.IOError([msg::AbstractString]) <: Exception

There was a problem reading or writing a file. `msg` is an optional descriptive error 
message.

See also [`FileExistsError`](@ref).
"""
struct IOError{T<:AbstractString} <: Exception
    msg::T
    function IOError{T}(msg::T) where T<:AbstractString
        new(string("IOError: ", msg))
    end
end

IOError(msg::T="") where T<:AbstractString = IOError{T}(msg)

Base.showerror(io::IO, e::IOError) = print(io, e.msg)


"""
    JuChrom.InputOutput.FileExistsError([msg::AbstractString]) <: Exception

File already exists. `msg` is an optional descriptive error message.

See also [`IOError`](@ref).
"""
struct FileExistsError{T<:AbstractString} <: Exception
    msg::String
    function FileExistsError{T}(msg::T) where T<:AbstractString
        new(string("FileExistsError: ", msg))
    end
end


FileExistsError(msg::T="") where T<:AbstractString = FileExistsError{T}(msg)


Base.showerror(io::IO, e::FileExistsError) = print(io, e.msg)
