PROJNAME = im
LIBNAME = imlua

IM = ..

OPT = YES
# To not link with the Lua dynamic library in UNIX
NO_LUALINK = Yes
# To use a subfolder with the Lua version for binaries
LUAMOD_DIR = Yes
USE_BIN2C_LUA = Yes
NO_LUAOBJECT = Yes
DEPENDDIR = dep

USE_IM = YES

SRC = lua5/imlua.c lua5/imlua_aux.c lua5/imlua_convert.c lua5/imlua_file.c \
      lua5/imlua_image.c lua5/imlua_palette.c lua5/imlua_util.c
DEF_FILE = lua5/imlua.def

SRCLUA = lua5/im_image.lua lua5/im_convert.lua
SRCLUADIR = lua5

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
ifdef NO_LUAOBJECT
  DEFINES += IMLUA_USELH
  USE_LH_SUBDIR = Yes
  LHDIR = lua5/lh
else
  DEFINES += IMLUA_USELOH
  USE_LOH_SUBDIR = Yes
  LOHDIR = lua5/loh$(LUASFX)
endif
