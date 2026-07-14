@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
  start "TRIAD Local Server" /min py -m http.server 8765 --bind 127.0.0.1
  timeout /t 1 /nobreak >nul
  start "" "http://127.0.0.1:8765/OPEN_LAB.html"
  exit /b 0
)
where python >nul 2>nul
if %errorlevel%==0 (
  start "TRIAD Local Server" /min python -m http.server 8765 --bind 127.0.0.1
  timeout /t 1 /nobreak >nul
  start "" "http://127.0.0.1:8765/OPEN_LAB.html"
  exit /b 0
)
start "" "%~dp0OPEN_LAB.html"
