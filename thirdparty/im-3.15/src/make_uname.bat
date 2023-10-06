@echo off
REM This builds all the libraries of the folder for 1 uname  

call tecmake %1 %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_process" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_process" "USE_OPENMP=Yes" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_jp2" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_avi" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_fftw3" %2 %3 %4 %5 %6 %7 %8

REM WMV and Capture are NOT available in some compilers,
REM so this may result in errors, just ignore them.
call tecmake %1 "MF=im_wmv" %2 %3 %4 %5 %6 %7 %8
call tecmake %1 "MF=im_capture" %2 %3 %4 %5 %6 %7 %8
  
REM Building for the default (USE_LUA51) 
REM or for the defined at the environment

REM call make_uname_lua %1 %2 %3 %4 %5 %6 %7 %8
call make_uname_lua %1 "USE_LUA51=Yes" %2 %3 %4 %5 %6 %7 %8
call make_uname_lua %1 "USE_LUA52=Yes" %2 %3 %4 %5 %6 %7 %8
call make_uname_lua %1 "USE_LUA53=Yes" %2 %3 %4 %5 %6 %7 %8
call make_uname_lua %1 "USE_LUA54=Yes" %2 %3 %4 %5 %6 %7 %8
