# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindEXIF
--------

Find the Joint Photographic Experts Group (EXIF) library (``libexif``)

Imported targets
^^^^^^^^^^^^^^^^

.. versionadded:: 3.12

This module defines the following :prop_tgt:`IMPORTED` targets:

``EXIF::EXIF``
  The EXIF library, if found.

Result variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``EXIF_FOUND``
  If false, do not try to use EXIF.
``EXIF_INCLUDE_DIRS``
  where to find libexif/exif-data.h, etc.
``EXIF_LIBRARIES``
  the libraries needed to use EXIF.

Cache variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``EXIF_INCLUDE_DIRS``
  where to find exiflib.h, etc.

#]=======================================================================]

find_path(EXIF_INCLUDE_DIR exif-data.h
    PATH_SUFFIXES libexif)

set(exif_names ${EXIF_NAMES} exif exif-static libexif libexif-static)
foreach(name ${exif_names})
  list(APPEND exif_names_debug "${name}d")
endforeach()

if(NOT EXIF_LIBRARY)
  find_library(EXIF_LIBRARY_RELEASE NAMES ${exif_names} NAMES_PER_DIR)
  find_library(EXIF_LIBRARY_DEBUG NAMES ${exif_names_debug} NAMES_PER_DIR)
  include(SelectLibraryConfigurations)
  select_library_configurations(EXIF)
  mark_as_advanced(EXIF_LIBRARY_RELEASE EXIF_LIBRARY_DEBUG)
endif()
unset(exif_names)
unset(exif_names_debug)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXIF
  REQUIRED_VARS EXIF_LIBRARY EXIF_INCLUDE_DIR
  VERSION_VAR EXIF_VERSION)

if(EXIF_FOUND)
  set(EXIF_LIBRARIES ${EXIF_LIBRARY})
  set(EXIF_INCLUDE_DIRS "${EXIF_INCLUDE_DIR}")

  if(NOT TARGET EXIF::EXIF)
    add_library(EXIF::EXIF UNKNOWN IMPORTED)
    if(EXIF_INCLUDE_DIRS)
      set_target_properties(EXIF::EXIF PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EXIF_INCLUDE_DIRS}")
    endif()
    if(EXISTS "${EXIF_LIBRARY}")
      set_target_properties(EXIF::EXIF PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${EXIF_LIBRARY}")
    endif()
    if(EXISTS "${EXIF_LIBRARY_RELEASE}")
      set_property(TARGET EXIF::EXIF APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(EXIF::EXIF PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
        IMPORTED_LOCATION_RELEASE "${EXIF_LIBRARY_RELEASE}")
    endif()
    if(EXISTS "${EXIF_LIBRARY_DEBUG}")
      set_property(TARGET EXIF::EXIF APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(EXIF::EXIF PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
        IMPORTED_LOCATION_DEBUG "${EXIF_LIBRARY_DEBUG}")
    endif()
  endif()
endif()

mark_as_advanced(EXIF_LIBRARY EXIF_INCLUDE_DIR)
