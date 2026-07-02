# Mapper containers

## JLD2 support

JuChrom also ships a `JLD2` extension so [`RetentionMapper`](@ref) objects can be stored 
and restored with [`JLD2`](https://github.com/JuliaIO/JLD2.jl). The extension loads 
automatically once `JLD2` is available.

Typical usage looks like:

```julia
using JLD2

save_object("retention_mapper.jld2", mapper)
mapper_loaded = load_object("retention_mapper.jld2")

jldsave("retention_mappers.jld2"; mapper, mapper_smooth)
mapper_reloaded = JLD2.load("retention_mappers.jld2", "mapper")
mapper_smooth_loaded = JLD2.load("retention_mappers.jld2", "mapper_smooth")
```

For a runnable persistence example, see the
[retention mapping workflow](mapping_workflow.md).

## Types

```@docs
JuChrom.AbstractRetentionMapper
JuChrom.RetentionMapper
```

## Metadata

```@docs
JuChrom.extras(::JuChrom.AbstractRetentionMapper)
```

## Calibration Anchors

```@docs
JuChrom.retentionunit_A
JuChrom.retentionunit_B
JuChrom.retentions_A
JuChrom.retentions_B
JuChrom.rawretentions_A
JuChrom.rawretentions_B
```

## Mapping Bounds

```@docs
JuChrom.mapmin
JuChrom.mapmax
JuChrom.invmapmin
JuChrom.invmapmax
JuChrom.rawmapmin
JuChrom.rawmapmax
JuChrom.rawinvmapmin
JuChrom.rawinvmapmax
```
