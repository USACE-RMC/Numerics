<#
.SYNOPSIS
Runs one approved BestFit regression method and appends independently checked TRX evidence.
.DESCRIPTION
Uses BestFit's guarded runner without changing its filter or scientific settings. A failed method
is recorded and then reported as an error so an outer serial loop stops for immediate diagnosis.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Test,
    [Parameter(Mandatory = $true)][string]$BestFitRoot,
    [Parameter(Mandatory = $true)][string]$NumericsRoot,
    [string]$Phase = 'iteration',
    [string]$CaseNote = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$evidenceRoot = Split-Path $PSScriptRoot -Parent
$pausePath = Join-Path $evidenceRoot 'artifacts\regression-repair\pause-before-next-method.txt'
if (Test-Path -LiteralPath $pausePath) {
    throw ('Paused before executing the next method: ' + (Get-Content -LiteralPath $pausePath -Raw).Trim())
}
$inventoryPath = Join-Path $evidenceRoot 'docs\distributions\bestfit-regression-inventory.json'
$ledgerPath = Join-Path $evidenceRoot 'docs\distributions\bestfit-regression-runs.jsonl'
$inventory = Get-Content -LiteralPath $inventoryPath -Raw | ConvertFrom-Json
$entry = @($inventory.methods | Where-Object method -eq $Test)
if ($entry.Count -ne 1) { throw "Not one exact approved method: $Test" }
$BestFitRoot = (Resolve-Path -LiteralPath $BestFitRoot).Path
$NumericsRoot = (Resolve-Path -LiteralPath $NumericsRoot).Path

. (Join-Path $PSScriptRoot 'bestfit-regression-evidence.ps1')

$numericsBefore = Get-BestFitRegressionSourceVersion $NumericsRoot 'Numerics'
$bestFitBefore = Get-BestFitRegressionSourceVersion $BestFitRoot @('src/RMC.BestFit','src/RMC.BestFit.Verification','src/TestCommon','verification/data')
$logDirectory = Join-Path $evidenceRoot 'artifacts\regression-repair\logs'
New-Item -ItemType Directory -Path $logDirectory -Force | Out-Null
$runId = (Get-Date -Format 'yyyyMMdd-HHmmss-fff') + '-' + $Test.Replace('.','_')
$logPath = Join-Path $logDirectory ($runId + '.log')
$started = [DateTimeOffset]::UtcNow
$stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
$failure = $null
$oldLocal = $env:UseLocalRmcNumerics
$oldProject = $env:RmcNumericsProjectPath
Push-Location $BestFitRoot
try {
    $env:UseLocalRmcNumerics = 'true'
    $env:RmcNumericsProjectPath = Join-Path $NumericsRoot 'Numerics\Numerics.csproj'
    & (Join-Path $BestFitRoot 'scripts\run-verification-test.ps1') -Test $Test *>&1 | Tee-Object -FilePath $logPath
} catch {
    $failure = $_.Exception.ToString()
    $failure | Add-Content -LiteralPath $logPath
} finally {
    Pop-Location
    $env:UseLocalRmcNumerics = $oldLocal
    $env:RmcNumericsProjectPath = $oldProject
    $stopwatch.Stop()
}

$logs = Get-Content -LiteralPath $logPath -Raw
$trxLine = [regex]::Match($logs, '(?m)^TRX: (.+)\s*$')
$record = [ordered]@{
    runId=$runId; method=$Test; phase=$Phase; startedUtc=$started.ToString('o');
    wallSeconds=$stopwatch.Elapsed.TotalSeconds; testSeconds=$null; buildSeconds=$null;
    result='InfrastructureFailure'; trx=$null; log=$logPath; exactIdentityVerified=$false;
    numerics=$numericsBefore; bestFit=$bestFitBefore; numericsAssemblySha256=$null;
    bestFitAssemblySha256=$null; verificationAssemblySha256=$null;
    dependencies=$entry[0].declaredDependencies; sourceChangedDuringRun=$false;
    failure=$failure; note=$CaseNote
}
$buildTime = [regex]::Match($logs, 'Time Elapsed (\d+:\d+:\d+(?:\.\d+)?)')
if ($buildTime.Success) { $record.buildSeconds = [TimeSpan]::Parse($buildTime.Groups[1].Value,[cultureinfo]::InvariantCulture).TotalSeconds }
if ($trxLine.Success) {
    $trxPath = $trxLine.Groups[1].Value.Trim()
    $trxFiles = @(Get-ChildItem -LiteralPath (Split-Path $trxPath) -Filter *.trx -Recurse -File)
    if ($trxFiles.Count -ne 1) { throw 'Expected exactly one TRX in the guarded run directory.' }
    [xml]$trx = Get-Content -LiteralPath $trxPath -Raw
    $results = @($trx.SelectNodes("//*[local-name()='UnitTestResult']"))
    $definitions = @($trx.SelectNodes("//*[local-name()='TestMethod']"))
    if ($results.Count -ne 1 -or $definitions.Count -ne 1) { throw 'Expected one TRX result and one test definition.' }
    $actual = $definitions[0].className + '.' + $definitions[0].name
    if ($actual -ne $Test) { throw "TRX identity mismatch: $actual" }
    $record.trx = $trxPath
    $record.exactIdentityVerified = $true
    $record.result = [string]$results[0].outcome
    $record.testSeconds = [TimeSpan]::Parse($results[0].duration,[cultureinfo]::InvariantCulture).TotalSeconds
    $errorNode = $results[0].SelectSingleNode(".//*[local-name()='Message']")
    if ($null -ne $errorNode) { $record.failure=$errorNode.InnerText }
    $assemblyRoot = Join-Path $BestFitRoot 'src\RMC.BestFit.Verification\bin\Debug\net10.0'
    $numericsAssembly = Join-Path $NumericsRoot 'Numerics\bin\Debug\net10.0\Numerics.dll'
    $record.numericsAssemblySha256=(Get-FileHash -LiteralPath (Join-Path $assemblyRoot 'Numerics.dll') -Algorithm SHA256).Hash
    if ($record.numericsAssemblySha256 -ne (Get-FileHash -LiteralPath $numericsAssembly -Algorithm SHA256).Hash) { throw 'Verification did not load the requested Numerics build.' }
    $record.bestFitAssemblySha256=(Get-FileHash -LiteralPath (Join-Path $assemblyRoot 'RMC.BestFit.dll') -Algorithm SHA256).Hash
    $record.verificationAssemblySha256=(Get-FileHash -LiteralPath (Join-Path $assemblyRoot 'RMC.BestFit.Verification.dll') -Algorithm SHA256).Hash
}
$numericsAfter = Get-BestFitRegressionSourceVersion $NumericsRoot 'Numerics'
$bestFitAfter = Get-BestFitRegressionSourceVersion $BestFitRoot @('src/RMC.BestFit','src/RMC.BestFit.Verification','src/TestCommon','verification/data')
$record.sourceChangedDuringRun = $numericsBefore.sourceSha256 -ne $numericsAfter.sourceSha256 -or $bestFitBefore.sourceSha256 -ne $bestFitAfter.sourceSha256
if ($null -ne $failure -and $record.result -eq 'Passed') {
    # A passing TRX cannot override a nonzero guarded-runner exit or infrastructure exception.
    $record.result = 'InfrastructureFailure'
}
$record | ConvertTo-Json -Depth 10 -Compress | Add-Content -LiteralPath $ledgerPath -Encoding utf8
Write-Host "Ledger: $($record.result); test=$($record.testSeconds)s; build=$($record.buildSeconds)s; source changed=$($record.sourceChangedDuringRun)"
if ($record.result -ne 'Passed' -or !$record.exactIdentityVerified -or $record.sourceChangedDuringRun) { throw "Method requires investigation: $Test" }
