param(
    [int]$WaitForPid = 0,
    [int]$Workers = 6
)

$ErrorActionPreference = "Stop"

$repositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
$resultsRoot = Join-Path $repositoryRoot "invariant_region_mlp\results\mechanism_suite_r3"
$logPath = Join-Path $resultsRoot "pipeline.log"
$statusPath = Join-Path $resultsRoot "pipeline_status.json"

New-Item -ItemType Directory -Force -Path $resultsRoot | Out-Null
Set-Location -LiteralPath $repositoryRoot

function Write-PipelineLog {
    param([string]$Message)
    $timestamp = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
    "$timestamp $Message" | Tee-Object -FilePath $logPath -Append
}

function Invoke-PythonStage {
    param([string[]]$Arguments)
    Write-PipelineLog ("START python " + ($Arguments -join " "))
    & python @Arguments 2>&1 | Tee-Object -FilePath $logPath -Append
    if ($LASTEXITCODE -ne 0) {
        throw "python exited with code ${LASTEXITCODE}: $($Arguments -join ' ')"
    }
    Write-PipelineLog ("DONE python " + ($Arguments -join " "))
}

try {
    if ($WaitForPid -gt 0 -and (Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue)) {
        Write-PipelineLog "Waiting for active S1 process $WaitForPid"
        Wait-Process -Id $WaitForPid
        Write-PipelineLog "Active S1 process $WaitForPid exited"
    }

    $module = "invariant_region_mlp.experiments.mechanism_suite_r3"
    Invoke-PythonStage @("-m", "$module.run", "--stage", "s1", "--workers", "$Workers")
    Invoke-PythonStage @("-m", "$module.run", "--stage", "s3", "--workers", "$Workers")
    Invoke-PythonStage @("-m", "$module.run", "--stage", "s4", "--workers", "$Workers")
    Invoke-PythonStage @("-m", "$module.analyze")
    Invoke-PythonStage @("-m", "unittest", "$module.test_round3")

    $status = [ordered]@{
        status = "complete"
        completed_utc = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
        log = $logPath
    }
    $status | ConvertTo-Json | Set-Content -LiteralPath $statusPath -Encoding utf8
    Write-PipelineLog "PIPELINE COMPLETE"
}
catch {
    $status = [ordered]@{
        status = "failed"
        failed_utc = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
        error = $_.Exception.Message
        log = $logPath
    }
    $status | ConvertTo-Json | Set-Content -LiteralPath $statusPath -Encoding utf8
    Write-PipelineLog ("PIPELINE FAILED: " + $_.Exception.Message)
    exit 1
}
