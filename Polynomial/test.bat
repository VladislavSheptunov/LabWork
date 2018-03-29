@echo off

set PARAM=%1

cd Polynomial

IF "%PARAM%" NEQ "" (
	IF "%PARAM%" NEQ "-p2.7" (
		IF "%PARAM%" NEQ "-p3.6" (
			echo Invalid parameters
			)
		)
) ELSE (
	echo No arguments
) 

 
if %PARAM% == -p2.7 (
	C:\Python27\python ./utest.py
	)

if %PARAM% == -p3.6 (
	C:\ProgramData\Anaconda3\python ./utest.py
	) 

