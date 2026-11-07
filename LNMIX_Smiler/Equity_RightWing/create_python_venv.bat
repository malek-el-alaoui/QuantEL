@echo off
setlocal
cd /d "%~dp0"
py -m venv venv
if errorlevel 1 exit /b 1
call venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r python\requirements.txt
echo.
echo Python environment ready. Run: venv\Scripts\python.exe python\plot_two_models.py
endlocal
