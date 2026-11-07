@echo off
REM Run from the FX_FlatMixture project directory after building Release x64.
REM The executable accepts: <dataset> <eta>.  FX uses NO equity wing control.
set EXE=x64\Release\VolatilitySmileTwoModels.exe
if not exist "%EXE%" set EXE=Release\VolatilitySmileTwoModels.exe

%EXE% eurusd 0.10
if exist output_eta10_eurusd rmdir /s /q output_eta10_eurusd
xcopy /e /i /y output output_eta10_eurusd >nul

%EXE% eurusd 0.02
if exist output_eta02_eurusd rmdir /s /q output_eta02_eurusd
xcopy /e /i /y output output_eta02_eurusd >nul

%EXE% usdhkd 0.10
if exist output_eta10_usdhkd rmdir /s /q output_eta10_usdhkd
xcopy /e /i /y output output_eta10_usdhkd >nul

%EXE% usdhkd 0.02
if exist output_eta02_usdhkd rmdir /s /q output_eta02_usdhkd
xcopy /e /i /y output output_eta02_usdhkd >nul

echo Done.
