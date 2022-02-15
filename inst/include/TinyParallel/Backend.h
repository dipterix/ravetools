
#ifndef __RAVETOOLS_PARALLEL_BACKEND__
#define __RAVETOOLS_PARALLEL_BACKEND__

#include <cstdlib>
#include <cstring>

extern "C" {
void REprintf(const char*, ...);
}

namespace TinyParallel {
namespace internal {

enum backend_type {
   BACKEND_TINYTHREAD
};

inline backend_type defaultBackend()
{
   return BACKEND_TINYTHREAD;
}

inline const char* backendToString(backend_type backend)
{
   return "tinythread";
}

inline backend_type backend()
{
   const char* requestedBackend = std::getenv("RAVETOOLS_BACKEND");
   if (requestedBackend == NULL)
   {
      return defaultBackend();
   }
   else if (std::strcmp(requestedBackend, "tbb") == 0)
   {
      const char* msg =
         "tbb backend is not available; using tinythread instead";

      REprintf("%s\n", msg);
      return BACKEND_TINYTHREAD;
   }
   else if (strcmp(requestedBackend, "tinythread") == 0)
   {
      return BACKEND_TINYTHREAD;
   }
   else
   {
      const char* fmt = "unknown parallel backend '%s'; using '%s' instead\n";
      REprintf(fmt, requestedBackend, backendToString(defaultBackend()));
      return defaultBackend();
   }
}

} // namespace internal
} // namespace TinyParallel

#endif /* __RAVETOOLS_PARALLEL_BACKEND__ */
