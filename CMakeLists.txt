cmake_minimum_required (VERSION 3.21)

set (IUP_VERSION "3.30.0")
project ("IUP-GUI"
    LANGUAGES
    VERSION ${IUP_VERSION})

option (USE_HAIKU "Add support for Haiku window manager (BeOS)" OFF)
option (USE_GTK "" OFF)
option (USE_MOTIF "" OFF)
option (USE_WIN32 "Add support for Win32 API" ON)
option (USE_NEW_DRAW "Draw driver with alpha and anti-aliasing in Windows and Linux enabled" ON)

add_subdirectory ("src")
add_subdirectory ("srcimglib")
add_subdirectory ("examples")
