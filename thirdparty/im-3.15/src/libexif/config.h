/* essencial config.h for IM */

#define GETTEXT_PACKAGE "libexif"

#ifdef WIN32
#define snprintf _snprintf
#if _MSC_VER >= 1900 /* IMLIB vc14 */
#undef snprintf
#endif
#if __GNUC__ >= 6 /* IMLIB mingw6 */
#undef snprintf
#endif
#endif

/* Others:

- removed inline from code 
- replaced ssize_t with size_t
- removed , in last enum at *-tag.h

*/
