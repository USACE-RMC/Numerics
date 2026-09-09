<#
.SYNOPSIS
Resumes the approved BestFit allowlist one exact method at a time, stopping for any failure.
.DESCRIPTION
Only the latest exact, passing, unchanged-source evidence for the requested checkout pair is
reused. Failed, absent, or stale methods run through the guarded single-method runner. The
inventory retains excluded and locally absent methods without expanding execution to them.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$BestFitRoot,
    [Parameter(Mandatory = $true)][string]$NumericsRoot,
    [string[]]$Methods,
    [switch]$ForceRerun
)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'bestfit-regression-evidence.ps1')
$evidenceRoot = Split-Path $PSScriptRoot -Parent
$inventory = Get-Content -LiteralPath (Join-Path $evidenceRoot 'docs\distributions\bestfit-regression-inventory.json') -Raw | ConvertFrom-Json
if ($inventory.methods.Count -ne 293) { throw 'The approved inventory must contain exactly 293 methods.' }
$BestFitRoot = (Resolve-Path -LiteralPath $BestFitRoot).Path
$NumericsRoot = (Resolve-Path -LiteralPath $NumericsRoot).Path
$numericsVersion = Get-BestFitRegressionSourceVersion $NumericsRoot 'Numerics'
$bestFitVersion = Get-BestFitRegressionSourceVersion $BestFitRoot @('src/RMC.BestFit','src/RMC.BestFit.Verification','src/TestCommon','verification/data')
$ledgerPath = Join-Path $evidenceRoot 'docs\distributions\bestfit-regression-runs.jsonl'
$latest = @{}
if (Test-Path -LiteralPath $ledgerPath) {
    foreach ($line in Get-Content -LiteralPath $ledgerPath) {
        $record = $line | ConvertFrom-Json
        if ($null -eq $record.numerics.PSObject.Properties['sourceSha256'] -or $null -eq $record.bestFit.PSObject.Properties['sourceSha256']) { continue }
        if ($record.numerics.root -eq $NumericsRoot -and $record.bestFit.root -eq $BestFitRoot -and
            $record.numerics.sourceSha256 -eq $numericsVersion.sourceSha256 -and $record.bestFit.sourceSha256 -eq $bestFitVersion.sourceSha256) {
            $latest[$record.method] = $record
        }
    }
}
if ($Methods) {
    if (@($Methods | Sort-Object -Unique).Count -ne $Methods.Count) { throw 'Duplicate requested methods are not allowed.' }
    foreach ($method in $Methods) { if ($method -notin $inventory.methods.method) { throw "Method is outside the approved inventory: $method" } }
    $ordered = @($Methods)
} else {
    $ordered = @($inventory.methods | Sort-Object @{Expression={
        if ($_.method -match 'FittingAnalysisRecoveryTests\.(PearsonTypeIII|LogPearsonTypeIII|LnNormal)_') { 0 }
        elseif ($_.method -match '\.ArrFlikeTests\.') { 1 }
        elseif ($_.method -match 'Pearson|LP3|B17CCovarianceTests|NonstationaryParentTrendCoverageTests') { 2 }
        elseif ($_.method -match 'CompetingRiskRecoveryTests|MixtureRecoveryTests|CompositePredictiveRecoveryTests') { 3 }
        elseif ($_.area -eq 'DistributionFitting') { 4 }
        elseif ($_.method -match '\.ViglioneEtAlTests\.') { 5 }
        elseif ($_.area -eq 'ModelEstimation') { 6 }
        elseif ($_.area -eq 'Univariate') { 7 }
        elseif ($_.area -eq 'Bivariate') { 8 }
        elseif ($_.area -eq 'RatingCurve') { 9 }
        elseif ($_.area -eq 'TimeSeries') { 10 }
        else { 11 }
    }},method | Select-Object -ExpandProperty method)
}
$pending = @($ordered | Where-Object {
    $record = $latest[$_]
    $ForceRerun -or $null -eq $record -or $record.result -ne 'Passed' -or !$record.exactIdentityVerified -or $record.sourceChangedDuringRun
})
$alreadyPassed = $ordered.Count - $pending.Count
Write-Output "$alreadyPassed current passes retained; $($pending.Count) exact methods pending out of $($ordered.Count) selected."
$logDirectory = Join-Path $evidenceRoot 'artifacts\regression-repair\logs'
New-Item -ItemType Directory -Path $logDirectory -Force | Out-Null
$batchLog = Join-Path $logDirectory ((Get-Date -Format 'yyyyMMdd-HHmmss-fff') + '-allowlist.log')
$index = 0
foreach ($method in $pending) {
    $index++
    Write-Output "[$index/$($pending.Count)] $method"
    try {
        & (Join-Path $PSScriptRoot 'run-bestfit-regression-method.ps1') -Test $method -BestFitRoot $BestFitRoot -NumericsRoot $NumericsRoot -Phase 'stabilized-pass' -CaseNote 'Resumable approved allowlist; no scientific settings altered.' *>&1 | Out-File -LiteralPath $batchLog -Append -Encoding utf8
        $record = Get-Content -LiteralPath $ledgerPath -Tail 1 | ConvertFrom-Json
        if ($record.method -ne $method) { throw 'The latest ledger identity does not match the completed method.' }
        if ($record.numerics.sourceSha256 -ne $numericsVersion.sourceSha256 -or $record.bestFit.sourceSha256 -ne $bestFitVersion.sourceSha256) { throw 'Source changed between methods; restart to reconcile stale passes.' }
        Write-Output "$($record.result): test $($record.testSeconds)s; build $($record.buildSeconds)s; $($alreadyPassed + $index)/$($ordered.Count) current."
    } catch {
        if (Test-Path -LiteralPath $batchLog) { Get-Content -LiteralPath $batchLog -Tail 28 }
        throw
    }
}
Write-Output "All $($ordered.Count) selected exact methods have passing evidence for this source fingerprint pair."
