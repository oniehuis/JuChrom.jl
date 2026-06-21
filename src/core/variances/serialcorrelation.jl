"""
    vif(ρ::Real, n::Integer; ρmax=0.8, nonnegative=true, nmin=1)

Return the variance inflation factor associated with lag correlation `ρ` and sample
size `n`:

    VIF = 1 + 2 * ρ * (n_eff - 1) / n_eff

where `n_eff = max(n, nmin)`. Non-finite correlations are treated as zero. When
`nonnegative=true`, negative correlations are clipped to zero; otherwise correlations are
clamped symmetrically to `[-ρmax, ρmax]`. The returned factor is always at least one.

See also
[`binretentions`](@ref)
"""
vif(ρ::Real, n::Integer; kwargs...) = vif(float(ρ), n; kwargs...)

function vif(
    ρ::T1,
    n::Integer;
    ρmax::Real=0.8,
    nonnegative::Bool=true,
    nmin::Integer=1,
) where {T1<:AbstractFloat}
    0 ≤ ρmax ≤ 1 || throw(ArgumentError("ρmax must be in [0, 1]"))
    nmin ≥ 1 || throw(ArgumentError("nmin must be ≥ 1"))

    ρcap = T1(ρmax)
    neff = T1(max(n, nmin))
    ρ₀ = isfinite(ρ) ? ρ : zero(T1)
    ρ₁ = nonnegative ? max(ρ₀, zero(T1)) : ρ₀
    ρ₂ = nonnegative ? min(ρ₁, ρcap) : clamp(ρ₁, -ρcap, ρcap)

    max(one(T1) + T1(2) * ρ₂ * (neff - one(T1)) / neff, one(T1))
end
