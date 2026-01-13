# JuChrom.jl

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] [![][codecov-img]][codecov-url] |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://oniehuis.github.io/JuChrom.jl/dev

[GHA-img]: https://github.com/oniehuis/JuChrom.jl/actions/workflows/CI.yml/badge.svg?branch=main
[GHA-url]: https://github.com/oniehuis/JuChrom.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/oniehuis/JuChrom.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/oniehuis/JuChrom.jl

## Testing

Run tests with the dedicated test environment to match CI:

```bash
julia --project=test -e 'import Pkg; Pkg.test()'
```
