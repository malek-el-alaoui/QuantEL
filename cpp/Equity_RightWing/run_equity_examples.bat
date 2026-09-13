@echo off
REM Run from the Equity_RightWing project directory after building Release x64.
REM The executable accepts: <dataset> <eta>.
set EXE=x64\Release\VolatilitySmileTwoModels.exe
if not exist "%EXE%" set EXE=Release\VolatilitySmileTwoModels.exe

%EXE% all 0.10
if exist output_eta10 rmdir /s /q output_eta10
xcopy /e /i /y output output_eta10 >nul

%EXE% all 0.02
if exist output_eta02 rmdir /s /q output_eta02
xcopy /e /i /y output output_eta02 >nul

echo Done.
