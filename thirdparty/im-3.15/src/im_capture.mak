PROJNAME = im
LIBNAME = im_capture
OPT = YES
             
INCLUDES = ../include
DEPENDDIR = dep

ifneq ($(findstring Win, $(TEC_SYSNAME)), )
  # We use ISampleGrabberCB from qedit.h which is available only in an old Windows SDK (6.1 - Vista)
  # and also need an old DirectX SDK (9.15).
  DXSDK = d:/lng/dxsdk
  WINSDK = d:/lng/winsdk
  
  INCLUDES += $(WINSDK)/include $(DXSDK)/include 
  SRC = im_capture_dx.cpp
endif

#ifneq ($(findstring Linux, $(TEC_UNAME)), )
#  SRC = im_capture_v4l.cpp
#endif
             
LIBS = strmiids

#mingw3-dll:                    
#	@echo Importing MingW stub library
#	@cd ../lib/dll
#	@dlltool -d im_capture.def -D im_capture.dll -l ../lib/mingw3/libim_capture.a
#	@cd ../src

#VC6 DLL must be available
mingw4-dll:                    
	@echo Importing MingW stub library
	@cd ../lib/dll
	@dlltool -d im_capture.def -D im_capture.dll -l ../lib/mingw4/libim_capture.a
	@cd ../src
	@cp -f ../lib/dll/im_capture.dll ../lib/mingw4/

#VC6 DLL must be available
dllw4-dll:                    
	@echo Importing MingW stub library
	@cd ../lib/dll
	@dlltool -d im_capture.def -D im_capture.dll -l ../lib/dllw4/libim_capture.a
	@cd ../src
	@cp -f ../lib/dll/im_capture.dll ../lib/dllw4/

#bc56-dll:                    
#	@echo Importing Bcc stub library
#	@d:/lng/cbuilderx/bin/implib -a ../lib/bc56/im_capture.lib ../lib/dll/im_capture.dll

#owc1-dll:                    
#	@wlib -b -c -n -q -fo -io ../lib/owc1/im_capture.lib @im_capture.wlib
# TEST	@wlib -b -c -n -q -fo -io ../lib/owc1/im_capture.lib +../lib/dll/im_capture.dll


ifneq ($(findstring gcc, $(TEC_UNAME)), )
  $(error No support for DirectX in Cygwin)
endif
ifneq ($(findstring mingw, $(TEC_UNAME)), )
  $(error No support for DirectX in MingW)
endif
ifneq ($(findstring cygw, $(TEC_UNAME)), )
  $(error No support for DirectX in Cygwin)
endif
ifneq ($(findstring dllw, $(TEC_UNAME)), )
  $(error No support for DirectX in MingW)
endif
ifneq ($(findstring dllg, $(TEC_UNAME)), )
  $(error No support for DirectX in Cygwin)
endif
ifneq ($(findstring owc, $(TEC_UNAME)), )
  $(error No support for DirectX in OpenWatcom)
endif
ifneq ($(findstring bc, $(TEC_UNAME)), )
  $(error No support for DirectX in BorlandC)
endif
