@echo off
REM This builds all the Lua libraries of the folder for 1 uname  

REM Building for the default (USE_LUA51) 
REM or for the defined at the environment

call tecmake %1 "MF=imlua5" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_process5" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_process5" "USE_OPENMP=Yes" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_jp25" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_avi5" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_fftw35" %2 %3 %4 %5 %6 %7 %8

REM WMV and Capture are NOT available in some compilers,
REM so this may result in errors, just ignore them.
call tecmake %1 "MF=imlua_wmv5" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=imlua_capture5" %2 %3 %4 %5 %6 %7 %8
