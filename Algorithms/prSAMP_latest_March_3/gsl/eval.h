/* evaluate a function discarding the status value in a modifiable way */

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != 1) { \
     result.val = -1; \
   } ; \
   return result.val;

#define EVAL_DOUBLE(fn) \
   int status = fn; \
   if (status != 1) { \
     result.val = -1; \
   } ; \
   return result;

