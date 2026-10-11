# JuChromExplorer compatibility

`workflows/downstream.yml` tests the JuChrom revision checked out for each pull
request or push to `main` against JuChromExplorer's public `main` branch on Julia
1.12 and 1.13. For pull requests, GitHub's checkout defaults to the proposed merge
revision. The logs record both repositories' commits.

The test script creates a temporary environment from the explorer's test project,
overrides its JuChrom source with the checkout under review, and verifies the
loaded package paths before running the explorer's full test suite. It leaves
both repositories' project files and dependency pins unchanged. Linux uses Xvfb
and software OpenGL. No cross-repository dispatch token is needed for public
repositories.

## Enable the required GitHub check

First publish JuChromExplorer (including its tests) and this workflow, then let
the workflow pass once. In **JuChrom.jl → Settings → Branches**, edit the existing
`main` rule or add a classic branch protection rule for `main`:

1. Enable **Require a pull request before merging**. Approving reviews need not be
   required for a sole maintainer.
2. Enable **Require status checks to pass before merging** and add
   **JuChromExplorer compatibility**. Retain existing required checks.
3. Enable **Require branches to be up to date before merging**.
4. Enable **Do not allow bypassing the above settings** to apply the rules to
   administrators too, then save.

The required check aggregates both Julia jobs and fails if either fails, is
cancelled, or is skipped. Its name stays constant if the test matrix changes.
Repository files alone do not enable branch protection; this setting must be
saved on GitHub. Avoid path filters on this required workflow.

## Run locally

With a graphical desktop available, from the JuChrom checkout:

```sh
julia --startup-file=no .github/scripts/test_explorer.jl ../JuChromExplorer
```

JuChromExplorer's normal CI and documentation use its tested JuChrom commit.
These downstream tests deliberately override that commit to detect upcoming
incompatibilities. Coordinate breaking changes before merging and update the
explorer's package, test, and documentation source pins together after testing.
