# Numerics Local Development Rules

## Technical Authority and Algorithm Change Control (NON-NEGOTIABLE)

Haden Smith is the final technical and numerical authority for this repository. AI agents,
including Codex and Claude, must never change an established algorithm, formula,
probability-combination rule, numerical method, convergence rule, default tolerance,
random-seed behavior, clipping or normalization policy, or reference-result contract without
Haden Smith's explicit prior approval. A failing test authorizes diagnosis and a proposed fix
only; it does not authorize tuning or replacing the algorithm. When an agent's reasoning
conflicts with Haden Smith's direction or domain judgment, the agent must stop, present the
evidence, and request a decision. AI confidence is not technical authority.

## Quality Gate

- No compiler errors.
- No compiler warnings.
- Do not suppress warnings to make a build pass; fix the code, XML documentation, or invalid reference.
- Empty catch blocks are not allowed. Preserve exception context when rethrowing or wrapping.

## XML Documentation

- XML documentation is mandatory for all public API additions and changes.
- Public classes, structs, interfaces, enums, constructors, properties, methods, events, and fields require accurate XML documentation.
- Private or internal methods with non-obvious algorithms should also be documented.
- Include `<summary>`, `<param>`, `<returns>`, `<exception>`, and `<remarks>` whenever they apply.
- Use `<remarks>` for numerical method notes, assumptions, references, and domain-specific constraints.
- XML documentation warnings `CS1570`, `CS1571`, `CS1572`, `CS1573`, `CS1574`, `CS1584`, `CS1587`, `CS1589`, and `CS1591` are build errors when `EnforceXmlDocumentation=true`.

## Unit Tests

- Every behavior change must include unit tests or an explicit reason no test is useful.
- Every new public class requires corresponding coverage in `Test_Numerics`.
- Statistical and numerical methods must cover known values, edge cases, invalid inputs, and regression cases.
- Use small inline fixtures in tests unless a shared fixture already exists for the exact scenario.

## Required Validation

Run validation before committing:

```powershell
dotnet build -c Release
dotnet test -c Release --no-build
```

Both commands must complete with zero errors, zero warnings, and zero failed tests. If validation fails, fix the issue before committing.

## Git Workflow

- Commit after successful build and test validation for the completed logical change.
- Stage only files that belong to the validated change.
- Do not include unrelated tracked changes or untracked local files.
- Use concise commit messages that describe the change.
- Do not push unless explicitly asked.

## PR and Release Workflow

Use this workflow for release-preparation pull requests and public package releases.

### Release PR Planning

- Base release PRs on `bug-fixes-and-enhancements` targeting `main` unless the user gives a different branch plan.
- Identify the last merged PR from the release branch and list every commit since that merge.
- Merge `origin/main` into the release branch before editing so current README, JOSS, citation, and metadata updates are preserved.
- Define the release outcome up front, including the exact `RMC.Numerics` version that NuGet.org must show as the latest public package after publishing.

### Version and Metadata Checklist

Update and cross-check every version-bearing file for each release:

- `Numerics/Numerics.csproj`: set `Version`, `AssemblyVersion`, and `PackageReleaseNotes`.
- `CITATION.cff`: set `version` and `date-released`.
- `codemeta.json`: set `version` and `dateModified`.
- `.github/workflows/Snapshot.yml`: advance the snapshot version to the next patch after the release.
- Keep `PackageReleaseNotes` concise because NuGet.org displays this field on the package page.

### PR Message Checklist

Use a title like `Prepare vX.Y.Z release`.

Include these sections in the PR body:

- `Summary`: package, citation, CodeMeta, and snapshot metadata updates.
- `Summary`: user-facing fixes, reliability changes, API changes, and test coverage added since the previous release branch merge.
- `Validation`: required commands and package inspection results.

Use this validation checklist:

```markdown
- [ ] `dotnet restore`
- [ ] `dotnet build Numerics/Numerics.csproj -c Release /p:Version=X.Y.Z`
- [ ] `dotnet test -c Release` with `VSTEST_CONNECTION_TIMEOUT=600`
- [ ] `dotnet pack Numerics/Numerics.csproj -c Release /p:Version=X.Y.Z --no-build -o ./packages`
- [ ] Inspect the `.nupkg` metadata and confirm version `X.Y.Z` and release notes are correct.
```

### Release Message Checklist

Create a GitHub Release titled `Numerics vX.Y.Z`.

Include these sections in the release body:

- `Highlights`: short release overview.
- `Reliability and Fixes`: behavior fixes and reliability improvements.
- Domain-specific sections for API, statistics, time-series, numerical methods, or other notable changes.
- `Install`: package install command.

Use this install command format:

```bash
dotnet add package RMC.Numerics --version X.Y.Z
```

### Publish Workflow

- Merge the validated PR to `main`.
- Tag the merge commit as `vX.Y.Z`.
- Push the tag so `.github/workflows/Release.yml` publishes to the internal Nexus feed.
- Create and publish a non-draft, non-prerelease GitHub Release for `vX.Y.Z`.
- Confirm `.github/workflows/NuGetPublish.yml` publishes `RMC.Numerics.X.Y.Z.nupkg` to NuGet.org.
- Verify NuGet.org shows `RMC.Numerics X.Y.Z` as the latest public package.
