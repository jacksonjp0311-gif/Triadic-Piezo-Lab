$ErrorActionPreference = 'SilentlyContinue'
Set-Location $PSScriptRoot
$python = Get-Command py -ErrorAction SilentlyContinue
if (-not $python) { $python = Get-Command python -ErrorAction SilentlyContinue }
if ($python) {
    Start-Process -WindowStyle Minimized -FilePath $python.Source -ArgumentList @('-m','http.server','8765','--bind','127.0.0.1')
    Start-Sleep -Milliseconds 800
    Start-Process 'http://127.0.0.1:8765/OPEN_LAB.html'
} else {
    Start-Process (Join-Path $PSScriptRoot 'OPEN_LAB.html')
}
