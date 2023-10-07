cmake_minimum_required (VERSION 3.21)

set (IUP_VERSION "3.30.0")
project ("IUP-GUI"
    LANGUAGES
    VERSION ${IUP_VERSION})

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

message (STATUS "System name      : ${CMAKE_SYSTEM_NAME}")
message (STATUS "System version   : ${CMAKE_SYSTEM_VERSION}")
message (STATUS "System processor : ${CMAKE_SYSTEM_PROCESSOR}")

if (CMAKE_SYSTEM_NAME MATCHES Linux)
    if (CMAKE_SYSTEM_VERSION VERSION_LESS_EQUAL 2.6)
        set (NO_GTK_DEFAULT ON)
    endif ()
endif ()

if (NOT IUP_TOOLKIT AND NOT NO_GTK_DEFAULT)
    if (CMAKE_SYSTEM_NAME MATCHES CYGWIN)
        set (GTK_DEFAULT ON)
    elseif (CMAKE_SYSTEM_NAME MATCHES CentOS)
        set (GTK_DEFAULT ON)
    elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
        message (STATUS "GTK_DEFAULT=ON")
        set (GTK_DEFAULT ON)
    elseif (CMAKE_SYSTEM_NAME MATCHES Darwin)
        set (GTK_DEFAULT ON)
    elseif (CMAKE_SYSTEM_NAME MATCHES FreeBSD)
        set (GTK_DEFAULT ON)
    elseif (CMAKE_SYSTEM_NAME MATCHES SunOS)
        if (CMAKE_SYSTEM_PROCESSOR MATCHES x86)
            set (GTK_DEFAULT ON)
        endif ()
    endif ()
endif ()

if (NOT IUP_TOOLKIT AND GTK_DEFAULT)
    if (CMAKE_SYSTEM_NAME MATCHES Linux)
        if (CMAKE_SYSTEM_VERSION VERSION_GREATER_EQUAL 3.1)
            set (IUP_TOOLKIT "GTK3")
        endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES CYGWIN)
        set (IUP_TOOLKIT "GTK3")
    # Homebrew
    elseif (CMAKE_SYSTEM_NAME MATCHES Darwin)
        if (CMAKE_SYSTEM_VERSION VERSION_GREATER_EQUAL 10)
            set (IUP_TOOLKIT "GTK3")
        endif ()
    else ()
        set (IUP_TOOLKIT "GTK2")
    endif ()
endif ()

if (NOT IUP_TOOLKIT)
    if (CMAKE_SYSTEM_NAME MATCHES Windows)
        set (IUP_TOOLKIT "Win32")
    endif ()
endif ()

include (CMakeDependentOption)

set (IUP_TOOLKIT "${IUP_TOOLKIT}" CACHE STRING "Select which GUI toolkit to use as a driver for IUP")
set_property (CACHE IUP_TOOLKIT PROPERTY STRINGS "OFF;GTK2;GTK3;Motif;Win32;Haiku")
set (IUP_LIB_TYPE "SHARED" CACHE STRING "Build IUP components as DLLs/shared libraries instead of static")
set_property (CACHE IUP_LIB_TYPE PROPERTY STRINGS "SHARED;STATIC")
cmake_dependent_option (IUP_GDK_NULL "Completely remove X11 calls (only works for GTK3)" OFF "IUP_TOOLKIT MATCHES GTK3" OFF)
cmake_dependent_option (IUP_NEW_DRAW "Draw driver with alpha and anti-aliasing in Windows and Linux enabled" ON "IUP_TOOLKIT" OFF)
cmake_dependent_option (IUP_IMGLIB "Build iupimglib: Library for loading and drawing images, icons and cursors" ON "IUP_TOOLKIT" OFF)
cmake_dependent_option (IUP_IMGLIB_LARGE "Enable for large icon support" ON "IUP_IMGLIB" OFF)
cmake_dependent_option (IUP_IM "Build image toolkit" ON "IUP_TOOLKIT" OFF)
cmake_dependent_option (IUP_CD "Build Canvas Draw library" ON "IUP_TOOLKIT" OFF)
cmake_dependent_option (IUP_VLED "Build VLED editor" ON "IUP_CD" OFF)
cmake_dependent_option (IUP_EXAMPLES "Build examples" ON "IUP_IMGLIB" OFF)

set (IUP_BUILD_BINDIR "bin" CACHE STRING "")
set (IUP_BUILD_LIBDIR "lib" CACHE STRING "")

set (IUP_INSTALL_BINDIR "bin" CACHE STRING "")
set (IUP_INSTALL_LIBDIR "lib" CACHE STRING "")
set (IUP_INSTALL_INCLUDEDIR "include" CACHE STRING "")
set (IUP_INSTALL_CMAKEDIR "share/iup/cmake" CACHE STRING "")

if (IUP_CD)
    if (NOT TARGET freetype)
        if (IUP_LIB_TYPE MATCHES "SHARED")
            set (FREETYPE_PIC ON CACHE INTERNAL "" FORCE)
        endif ()
        add_subdirectory ("thirdparty/freetype-2.10.2")
    endif ()
    if (NOT TARGET ftgl)
        if (IUP_LIB_TYPE MATCHES "SHARED")
            set (FTGL_PIC ON CACHE INTERNAL "" FORCE)
        endif ()
        add_subdirectory ("thirdparty/ftgl-2.1.5")
    endif ()
    if (NOT TARGET cd)
        if (IUP_LIB_TYPE MATCHES "SHARED")
            set (CD_PIC ON CACHE INTERNAL "" FORCE)
        endif ()
        add_subdirectory ("thirdparty/cd-5.14")
    endif ()
endif ()

if (IUP_CD OR IUP_IM)
    if (NOT TARGET ZLIB::ZLIB)
        if (IUP_LIB_TYPE MATCHES "SHARED")
            set (ZLIB_PIC ON CACHE INTERNAL "" FORCE)
        endif ()
        add_subdirectory ("thirdparty/zlib-1.2.11")
    endif ()
endif ()

if (IUP_TOOLKIT)
    add_subdirectory ("src")
endif ()
if (IUP_IMGLIB)
    add_subdirectory ("srcimglib")
endif ()
if (IUP_IM)
    if (IUP_LIB_TYPE MATCHES "SHARED")
        set (IM_PIC ON CACHE INTERNAL "" FORCE)
    endif ()
    add_subdirectory ("thirdparty/im-3.15")
    add_subdirectory ("srcim")
endif ()
if (IUP_CD)
    add_subdirectory ("srccd")
endif ()
if (IUP_VLED)
    #    add_subdirectory ("srcvled")
endif ()

if (IUP_EXAMPLES)
    add_subdirectory ("examples")
endif ()

# All subprojects have added targets to this export
# group, so make sure to install it
install (
    EXPORT IUPTargets
    NAMESPACE IUP::
    DESTINATION ${IUP_INSTALL_CMAKEDIR})
#export (
#    EXPORT IUPTargets
#    NAMESPACE IUP::
#    FILE ${IUP_BUILD_CMAKEDIR}/IUPTargets.cmake)

install (
    DIRECTORY "include/"
    DESTINATION ${IUP_INSTALL_INCLUDEDIR})
#install (
#    DIRECTORY "cmake/modules/"
#    DESTINATION ${IUP_INSTALL_CMAKEDIR}/modules)

#include(CMakePackageConfigHelpers)
#configure_package_config_file("cmake/templates/IUPConfig.cmake.in"
#    "${IUP_BUILD_CMAKEDIR}/IUPConfig.cmake"
#    INSTALL_DESTINATION "${IUP_INSTALL_CMAKEDIR}")
#install(
#    FILES "${IUP_BUILD_CMAKEDIR}/IUPConfig.cmake"
#    DESTINATION ${IUP_INSTALL_CMAKEDIR})
