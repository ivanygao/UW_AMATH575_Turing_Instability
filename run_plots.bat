@echo off
:menu
cls
echo =============================
echo Select the file to run:
echo -----------------------------
echo 1 - plot_figure5.py
echo 2 - plot_figure6.py
echo 3 - plot_figure7.py
echo 4 - plot_figure8.py
echo 5 - lbm.py
echo 6 - Run ALL
echo 0 - Exit
echo =============================
set /p choice=Enter your choice (0-5): 

if "%choice%"=="1" (
    echo Running Figure 5...
    python plot_figure5.py
    pause
    goto menu
)

if "%choice%"=="2" (
    echo Running Figure 6...
    python plot_figure6.py
    pause
    goto menu
)

if "%choice%"=="3" (
    echo Running Figure 7...
    python plot_figure7.py
    pause
    goto menu
)

if "%choice%"=="4" (
    echo Running Figure 8...
    python plot_figure8.py
    pause
    goto menu
)

if "%choice%"=="5" (
    echo Running lbm...
    python lbm.py
    pause
    goto menu
)

if "%choice%"=="6" (
    echo Running ALL figures...
    python plot_figure5.py
    python plot_figure6.py
    python plot_figure7.py
    python plot_figure8.py
    pause
    goto menu
)

if "%choice%"=="0" (
    echo Exiting...
    exit /b
)

echo Invalid input. Please try again.
pause
goto menu
