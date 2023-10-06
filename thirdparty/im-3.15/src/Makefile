ifeq ($(OS), Windows_NT)
  TECMAKE_CMD = $(MAKE) --no-print-directory -f ../tecmakewin.mak TEC_UNAME=$(TEC_UNAME)
  WINLIBS = im_avi im_wmv imlua_avi5 imlua_wmv5 im_capture imlua_capture5
else
  TECMAKE_CMD = $(MAKE) --no-print-directory -f ../tecmake.mak
endif

.PHONY: do_all im im_jp2 im_process im_fftw3 imlua5 imlua_jp25 imlua_process5 imlua_fftw35 $(WINLIBS)
do_all: im im_jp2 im_process im_process_omp im_fftw3 imlua5 imlua_jp25 imlua_process5 imlua_process_omp5 imlua_fftw35 $(WINLIBS)

im:
	@$(TECMAKE_CMD)
im_jp2:
	@$(TECMAKE_CMD) MF=im_jp2
im_avi:
	@$(TECMAKE_CMD) MF=im_avi
im_process:
	@$(TECMAKE_CMD) MF=im_process
im_process_omp:
	@$(TECMAKE_CMD) MF=im_process USE_OPENMP=Yes
im_fftw3:
	@$(TECMAKE_CMD) MF=im_fftw3
  
# Building for the default (USE_LUA51) 
# or for the defined at the environment

imlua5:
	@$(TECMAKE_CMD) MF=imlua5
imlua_jp25:
	@$(TECMAKE_CMD) MF=imlua_jp25
imlua_avi5:
	@$(TECMAKE_CMD) MF=imlua_avi5
imlua_process5:
	@$(TECMAKE_CMD) MF=imlua_process5
imlua_process_omp5:
	@$(TECMAKE_CMD) MF=imlua_process5 USE_OPENMP=Yes
imlua_fftw35:
	@$(TECMAKE_CMD) MF=imlua_fftw35
  
# WMV and Capture are NOT available in some compilers,
# so this may result in errors, just ignore them.
im_wmv:
	@$(TECMAKE_CMD) MF=im_wmv
imlua_wmv5:
	@$(TECMAKE_CMD) MF=imlua_wmv5 
im_capture:
	@$(TECMAKE_CMD) MF=im_capture
imlua_capture5:
	@$(TECMAKE_CMD) MF=imlua_capture5 
