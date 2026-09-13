@echo off
setlocal
cd /d "%~dp0"
if exist venv\Scripts\python.exe (
    venv\Scripts\python.exe python\plot_two_models.py
) else (
    py python\plot_two_models.py
)
endlocal
