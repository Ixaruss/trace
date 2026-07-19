#ifndef EXPORT_ADAPT_H
#define EXPORT_ADAPT_H

#ifdef _WIN32
    #ifdef ODLIB_EXPORTS
        #define ODLIB_API __declspec(dllexport)
    #else
        #define ODLIB_API __declspec(dllimport)
    #endif
#else
    #define ODLIB_API __attribute__ ((visibility ("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

ODLIB_API int LoadData(const int* flat_data, int rows, int cols);
ODLIB_API int DisplayData();

#ifdef __cplusplus
}
#endif

#endif // EXPORT_ADAPT_H
