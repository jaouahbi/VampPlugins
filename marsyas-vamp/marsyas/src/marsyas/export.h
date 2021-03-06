#ifndef _marsyas_DLLDEFINES_H
#define _marsyas_DLLDEFINES_H

#define MARSYAS_STATIC 1

#if defined (_WIN32) && !(MARSYAS_STATIC)

#ifdef _MSC_VER
#pragma warning(disable: 4251)
#endif

#if defined (marsyas_EXPORTS)
#define marsyas_EXPORT __declspec(dllexport)
#else
#define marsyas_EXPORT __declspec(dllimport)
#endif

#else

#define marsyas_EXPORT

#endif

#endif
