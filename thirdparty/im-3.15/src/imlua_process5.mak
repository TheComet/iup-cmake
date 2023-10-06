PROJNAME = im
LIBNAME = imlua_process

IM = ..

OPT = YES
# To not link with the Lua dynamic library in UNIX
NO_LUALINK = Yes
# To use a subfolder with the Lua version for binaries
LUAMOD_DIR = Yes
USE_BIN2C_LUA = Yes
NO_LUAOBJECT = Yes
DEPENDDIR = dep

USE_IMLUA = YES

SRC = lua5/imlua_process.c lua5/imlua_kernel.c lua5/imlua_convert.c
DEF_FILE = lua5/imlua_process.def

SRCLUA = lua5/im_process.lua lua5/im_processconvert.lua
SRCLUADIR = lua5

DEFINES += IM_PROCESS
INCLUDES = lua5
LIBS = im_process

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

ifdef USE_OPENMP
  LIBNAME := $(LIBNAME)_omp
  LIBS := im_process_omp
endif
LIBNAME := $(LIBNAME)$(LUASFX)

ifdef NO_LUAOBJECT
  DEFINES += IMLUA_USELH
  USE_LH_SUBDIR = Yes
  LHDIR = lua5/lh
else
  DEFINES += IMLUA_USELOH
  USE_LOH_SUBDIR = Yes
  LOHDIR = lua5/loh$(LUASFX)
endif

ifneq ($(findstring MacOS, $(TEC_UNAME)), )
  USE_IMLUA:=
  INCLUDES += ../include
  LDIR = ../lib/$(TEC_UNAME)
endif
