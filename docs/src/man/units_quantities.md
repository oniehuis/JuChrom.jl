# Units and Quantities

JuChrom supports optional units via
[Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl). You can work with unitless values
or Unitful quantities interchangeably, as long as each operation uses compatible units. This
page gives a short, package‑focused primer on the syntax you need for JuChrom.

## Quick primer

- Attach a unit to a numeric value with the `u""` string macro:
```@example 1
1.0u"s"        # (seconds)
```@example 2
12.5u"minute"
```
```@example 3
100.0u"Th"     # (Thomson)
```
```@example 4
5.3u"pA"       # (picoampere)
```
- Attach units to arrays in the same way:
```@example 5
ret = [1.0, 2.0, 3.0]u"s"
```
- Convert between compatible units with `uconvert`:
```@example 6
uconvert(u"minute", 120.0u"s")` → `2.0u"minute"
```
- Convert units for every element of a vector with broadcasting:
  - `uconvert.(u"minute", [30.0, 60.0, 90.0]u"s")` → `[0.5, 1.0, 1.5]u"minute"`
- Remove units when you need raw numbers:
  - `ustrip(12.5u"s")` → `12.5`
- Remove units from every element of a vector with broadcasting:
  - `ustrip.([30.0, 60.0, 90.0]u"s")` → `[30.0, 60.0, 90.0]`
- Unitful checks dimensional consistency for you:
  - `1.0u"s" + 500u"ms"` is valid
  - `1.0u"s" + 1.0u"Th"` raises an error

JuChrom accepts either unitless or unitful inputs for retentions, intensities, and m/z
values. If you pass unitful data, JuChrom preserves and propagates units where appropriate.
If you pass unitless data, the corresponding unit is stored as `nothing`.

## Common units in JuChrom

Below are the units you will see most often in this package:

- Time: `u"s"` (second), `u"ms"` (millisecond), `u"minute"`
- Current (intensity units for some detectors, e.g., FID): `u"A"` (ampere), `u"mA"`,
  `u"μA"`, `u"nA"`, `u"pA"`
- Mass‑to‑charge: `u"Th"` (Thomson; optional), or unitless m/z

### Metric prefixes

Unitful follows standard SI prefixes. The most common are:

- milli (`m`): `1e-3` (e.g., `u"ms"`, `u"mA"`)
- micro (`μ`): `1e-6` (e.g., `u"μs"`, `u"μA"`)
- nano (`n`): `1e-9` (e.g., `u"ns"`, `u"nA"`)
- pico (`p`): `1e-12` (e.g., `u"ps"`, `u"pA"`)

## The Thomson (`u"Th"`) in JuChrom

Mass‑to‑charge (m/z) is conventionally treated as unitless in many GC–MS workflows.
However, some instruments and data formats label m/z using the **Thomson** (symbol `Th`).
In JuChrom, `u"Th"` is supported as a convenience unit for m/z values:

- `mz = [57.0, 73.0, 91.0]u"Th"`

The Thomson is not an officially recognized derived unit in the International System of
Units (SI). It is used informally in mass spectrometry as shorthand for “m/z in daltons per
elementary charge (Da/e).” Use of `u"Th"` is optional; you can treat m/z as unitless if you
prefer. Use of `u"Th"` in JuChrom example code is primarily to illustrate support for
unitful m/z values.

## Further reading

For a full introduction to units and quantities in Julia, see the
[Unitful.jl documentation](https://painterqubits.github.io/Unitful.jl/stable/).
