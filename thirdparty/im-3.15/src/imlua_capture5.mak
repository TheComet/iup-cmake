PROJNAME = im
LIBNAME = imlua_capture
DEF_FILE = imlua_capture.def

OPT = YES
DEPENDDIR = dep

SRCDIR = lua5

SRC = imlua_capture.c

LIBS = im_capture

INCLUDES = lua5

ifdef USE_LUA_VERSION
  USE_LUA51:=
  USE_LUA52:=
  USE_LUA53:=
  USE_LUA54:=
  ifeq ($(USE_LUA_VERSION), 54)
    USE_LUA54:=Yes
  endif
  ifeq ($(USE_LUA_VERSION), 53)
    USE_LUA53:=Yes
  endif
  ifeq ($(USE_LUA_VERSION), 52)
    USE_LUA52:=Yes
  endif
  ifeq ($(USE_LUA_VERSION), 51)
    USE_LUA51:=Yes
  endif
endif

ifdef USE_LUA54
  LUASFX = 54
else
ifdef USE_LUA53
  LUASFX = 53
else
ifdef USE_LUA52
  LUASFX = 52
else
  USE_LUA51 = Yes
  LUASFX = 51
endif
endif
endif

LIBNAME := $(LIBNAME)$(LUASFX)

USE_IMLUA = Yes
# To not link with the Lua dynamic library in UNIX
NO_LUALINK = Yes
# To use a subfolder with the Lua version for binaries
LUAMOD_DIR = Yes
IM = ..

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
