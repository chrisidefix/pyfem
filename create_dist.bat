@echo off

@echo "Creating PyFem distribution"

RMDIR /S /Q buid

python setup.py bdist_wininst
pause
