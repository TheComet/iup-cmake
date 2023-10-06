/* IMLIB - jas_config.h for IM */

/* Avoid problems due to multiple inclusion. */
#ifndef JAS_CONFIG_H
#define JAS_CONFIG_H

#include <jasper/jas_dll.h>

/* This preprocessor symbol identifies the version of JasPer. */
#define JAS_VERSION "2.0.14"

#define JAS_HAVE_FCNTL_H		1
#define JAS_HAVE_SYS_TYPES_H	1

#undef JAS_HAVE_IO_H
#undef JAS_HAVE_WINDOWS_H
#undef JAS_HAVE_UNISTD_H
#undef JAS_HAVE_SYS_TIME_H
#undef JAS_HAVE_GETTIMEOFDAY
#undef JAS_HAVE_GETRUSAGE
#undef JAS_HAVE_SNPRINTF

#ifndef INT_FAST32_MAX
#define INT_FAST32_MAX   2147483647
#endif
#ifndef INT_FAST32_MIN
#define INT_FAST32_MIN   (-2147483647 - 1)
#endif
#ifndef SIZE_MAX
#define SIZE_MAX 0xffffffff
#endif

#if 0
#ifndef __cplusplus
#undef inline
#define inline __inline
#endif
#endif

#if !defined(JAS_DEC_DEFAULT_MAX_SAMPLES)
#define JAS_DEC_DEFAULT_MAX_SAMPLES (64 * ((size_t) 1048576))
#endif

#if defined(__GNUC__) && !defined(__clang__)
#define JAS_ATTRIBUTE_DISABLE_USAN \
  __attribute__((no_sanitize_undefined))
#elif defined(__clang__)
#define JAS_ATTRIBUTE_DISABLE_USAN \
  __attribute__((no_sanitize("undefined")))
#else
#define JAS_ATTRIBUTE_DISABLE_USAN
#endif

#endif
