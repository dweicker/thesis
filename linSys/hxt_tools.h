#ifndef HXT_TOOLS_H
#define HXT_TOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "hxt_message.h" // already include hxt_api (for HXT_status_t)


/*************************************************************
  Check for supported standard for memory allocation alignment
 **************************************************************/
#if (defined (__STD_VERSION__) && (__STD_VERSION__ >= 200112L))  || defined (_ISOC11_SOURCE)
#define HXT_HAVE_ALIGNED_MALLOC
#define HXT_HAVE_NORMAL_ALIGNED_FREE
#define HXT_C11 // in C11, function aligned_alloc(alignment,size) is supported
#elif (defined (__POSIX_VERSION) && (_POSIX_VERSION >= 200112L)) || (defined (_POSIX_C_SOURCE) && (_POSIX_C_SOURCE >= 200112L))
#define HXT_HAVE_ALIGNED_MALLOC
#define HXT_HAVE_NORMAL_ALIGNED_FREE
#define HXT_POSIX // in posix, function posix_memalign(ptr,alignment,size) is supported
#elif defined (_MSC_VER)
#include <malloc.h>
#define HXT_HAVE_ALIGNED_MALLOC
#define HXT_MICROSOFT // in microsoft, _aligned_malloc(size, alignment) is supported (it use _aligned_free to free and has an _aligned_realloc function)
#elif defined (__INTEL_COMPILER)
#include <malloc.h>
#define HXT_HAVE_ALIGNED_MALLOC
#define HXT_INTEL // intel support _mm_malloc and _mm_free()
#else
#warning "no supported aligned allocation"
#define HXT_HAVE_NORMAL_ALIGNED_FREE
#endif

/* define SIMD ALIGNMENT */
#ifndef SIMD_ALIGN
#ifdef HXT_HAVE_ALIGNED_MALLOC
#ifdef __AVX512F__
#define SIMD_ALIGN 64
#elif defined(__AVX2__)
#define SIMD_ALIGN 32
#else
#define SIMD_ALIGN 16 /* we align to 16 anyway even if no sse */
#endif
#else
#define SIMD_ALIGN (sizeof(size_t))
#endif
#endif

// SIMD_REST(x) = what to add to x to get a multiple of SIMD_ALIGN
#define SIMD_REST(x) ((SIMD_ALIGN -1) - (((x)-1) & (SIMD_ALIGN-1)))
#define SIMD_FLOOR(x) ((x) & ~(SIMD_ALIGN-1))

// declare alignement of pointer allocated on the stack or in a struct
#ifdef HXT_MICROSOFT
#define hxt_declare_aligned  __declspec(align(SIMD_ALIGN))
#ifndef __restrict__
#define __restrict__ __restrict
#endif
// no attribute in MSVC
//#ifndef __attribute__
//#define __attribute__(x)
//#endif
#else
#define hxt_declare_aligned  __attribute__((aligned(SIMD_ALIGN)))
#endif


#if !defined(__INTEL_COMPILER) || (__INTEL_COMPILER < 1300)
#define __assume(x)
#define __assume_aligned(x,y)
#endif


/*********************************************************
 * Hextreme malloc implementation
 *********************************************************/

/************************************************************************
  WARNING: the pointer passed (p) should be a double pointer !!

  use it like this:
___
|  int array_length = ...;
|  int *array;
|  HXT_CHECK(
|    HXT_malloc((void**) &array, sizeof(int)*array_length) );
|  
|  array[0] = ...;
|  [...]
|  array[array_length-1] = ...;
|
|  [...]
|
|  HXT_CHECK(
|    HXT_free(&array);
|___
 ************************************************************************/

static inline HXT_status_t HXT_malloc(void* ptr_to_ptr, size_t size)
{
  void** p = (void**)ptr_to_ptr;
  *p = malloc(size);
  if (*p==NULL)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  return HXT_STATUS_OK;
}


// allocate num element of size size and zero the memory
static inline HXT_status_t HXT_calloc(void* ptr_to_ptr, size_t num, size_t size)
{
  void** p = (void**)ptr_to_ptr;
  *p = calloc(num, size);
  if (*p==NULL)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  return HXT_STATUS_OK;
}


static inline HXT_status_t HXT_free(void* ptr_to_ptr)
{
  void** p = (void**)ptr_to_ptr;
  free(*p);
  *p=NULL;
  return HXT_STATUS_OK;
}


static inline HXT_status_t HXT_realloc(void* ptr_to_ptr, size_t size)
{
  void** p = (void**)ptr_to_ptr;
  void* newptr = realloc(*p, size);
  if (newptr==NULL && *p!=NULL && size!=0)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  *p = newptr;
  return HXT_STATUS_OK;
}



/*********************************************************
 * Hextreme aligned malloc implementation
 *********************************************************/
 // stuff to have an aligned malloc
static inline HXT_status_t HXT_aligned_malloc(void* ptr_to_ptr, size_t size)
{
  void** p = (void**)ptr_to_ptr;
#ifdef HXT_C11
  *p = aligned_alloc(SIMD_ALIGN, size);
#elif defined (HXT_POSIX)
  if (posix_memalign(p, SIMD_ALIGN, size) != 0) {
    *p = NULL;
  }
#elif defined (HXT_MICROSOFT)
  *p = _aligned_malloc(size, SIMD_ALIGN);
#elif defined (HXT_INTEL)
  *p = _mm_malloc(size, SIMD_ALIGN);
#else
  *p = malloc(size);
#endif
  if (*p == NULL)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  return HXT_STATUS_OK;
}


static inline HXT_status_t HXT_aligned_free(void* ptr_to_ptr)
{
  void** p = (void**)ptr_to_ptr;
#ifdef HXT_HAVE_NORMAL_ALIGNED_FREE
  free(*p);
#elif defined (HXT_MICROSOFT)
  _aligned_free(*p);
#elif defined (HXT_INTEL)
  _mm_free(*p);
#else
#error "you should not arrive here"
#endif
  *p = NULL;
  return HXT_STATUS_OK;
}


static inline HXT_status_t HXT_aligned_realloc(void* ptr_to_ptr, size_t size, size_t size_to_copy)
{
  void** p = (void**)ptr_to_ptr;
#ifdef HXT_MICROSOFT
  void* newptr = _aligned_realloc(*p, size, SIMD_ALIGN);
  if (newptr == NULL && *p != NULL && size != 0) {
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  }
#elif !defined(HXT_HAVE_ALIGNED_MALLOC)
  void* newptr = realloc(*p, size);
  if (newptr == NULL && *p != NULL && size != 0) {
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  }
#else
  if (size == 0) {
    HXT_CHECK(HXT_aligned_free(p));
    return HXT_STATUS_OK;
  }

  if (*p == NULL) {
    HXT_CHECK(HXT_aligned_malloc(p, size));
    return HXT_STATUS_OK;
  }

  void* newptr = NULL;
  HXT_CHECK(
    HXT_aligned_malloc(&newptr, size));

  memcpy(newptr, *p, size_to_copy);

  HXT_CHECK( HXT_aligned_free(p) );
#endif

  *p = newptr;
  return HXT_STATUS_OK;
}


HXT_status_t HXT_det3x3(double mat[3][3], double *det);
HXT_status_t HXT_inv3x3(double mat[3][3], double inv[3][3], double *det);

HXT_status_t HXT_inv4x4ColumnMajor(double mat[16], double inv[16], double *det);



#endif
