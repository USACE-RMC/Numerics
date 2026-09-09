# Shared source fingerprinting for the exact-method runner and resumable allowlist.
function Get-BestFitRegressionSourceVersion([string]$Root, [string[]]$SourceDirectory) {
    $head = & git -C $Root rev-parse HEAD
    if ($LASTEXITCODE -ne 0) { throw "Cannot resolve source version: $Root" }
    $files = @(& git -C $Root ls-files --cached --others --exclude-standard -- $SourceDirectory 'Directory.Build.*' 'Directory.Packages.props' 'global.json' 'NuGet.Config') | Sort-Object -Unique
    if ($LASTEXITCODE -ne 0) { throw "Cannot resolve source files: $Root" }
    $parts = foreach ($relative in $files) {
        $absolute = Join-Path $Root $relative
        if (Test-Path -LiteralPath $absolute -PathType Leaf) {
            $relative + ':' + (Get-FileHash -LiteralPath $absolute -Algorithm SHA256).Hash
        } else { $relative + ':absent' }
    }
    $bytes = [System.Text.Encoding]::UTF8.GetBytes(($parts -join "`n"))
    $sha = [System.Security.Cryptography.SHA256]::Create()
    try { $fingerprint = [Convert]::ToHexString($sha.ComputeHash($bytes)) } finally { $sha.Dispose() }
    return [ordered]@{ root=$Root; commit=[string]$head; sourceSha256=$fingerprint }
}
