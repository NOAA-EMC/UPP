#ifndef MS_SUA
# include <stdio.h>
# include <stdlib.h>
#endif
#include <string.h>
#include "streams.h"

#ifndef CRAY
# ifdef NOUNDERSCORE
#      define INT_GET_TI_HEADER_C  int_get_ti_header_c
#      define INT_GEN_TI_HEADER_C  int_gen_ti_header_c
# else
#   ifdef F2CSTYLE
#      define INT_GET_TI_HEADER_C  int_get_ti_header_c__
#      define INT_GEN_TI_HEADER_C  int_gen_ti_header_c__
#   else
#      define INT_GET_TI_HEADER_C  int_get_ti_header_c_
#      define INT_GEN_TI_HEADER_C  int_gen_ti_header_c_
#   endif
# endif
#endif

#ifdef MEMCPY_FOR_BCOPY
# define bcopy(A,B,C) memcpy((B),(A),(C))
#endif

int
INT_GEN_TI_HEADER_C ( char * hdrbuf, int * hdrbufsize,           /* hdrbufsize is in bytes */
                    int * itypesize, int * typesize,
                    int * DataHandle, char * Data,
                    int * Count, int * code )
{
  int i ;
  char * p ;
  p = hdrbuf ;
  p += sizeof(int) ;
  bcopy( code, p, sizeof(int) ) ; p += sizeof(int) ;       /* 2 */
  bcopy( DataHandle, p, sizeof(int) ) ; p += sizeof(int) ; /* 3 */
  bcopy( typesize, p, sizeof(int) ) ; p += sizeof(int) ;   /* 4 */
  bcopy( Count, p, sizeof(int) ) ; p += sizeof(int) ;      /* 5 */
  bcopy( Data, p, *Count * *typesize ) ; p += *Count * *typesize ; /* 6++ */
  *hdrbufsize = (int) (p - hdrbuf) ;
  bcopy( hdrbufsize, hdrbuf, sizeof(int) ) ;
  return(0) ;
}

int
INT_GET_TI_HEADER_C ( char * hdrbuf, int * hdrbufsize, int * n,  /* hdrbufsize and n are in bytes */
                    int * itypesize, int * typesize,
                    int * DataHandle, char * Data,
                    int * Count, int * code )
{
  int i ;
  char * p ;
  p = hdrbuf ;
  bcopy( p, hdrbufsize, sizeof(int) ) ;     p += sizeof(int) ;        /* 1 */
  bcopy( p, code, sizeof(int) ) ;           p += sizeof(int) ;        /* 2 */
  bcopy( p, DataHandle, sizeof(int) ) ;     p += sizeof(int) ;        /* 3 */
  bcopy( p, typesize, sizeof(int) ) ;       p += sizeof(int) ;        /* 4 */
  bcopy( p, Count, sizeof(int) ) ;          p += sizeof(int) ;        /* 5 */
  if ( *Count * *typesize > 0 ) {
  bcopy( p, Data, *Count * *typesize ) ;  p += *Count * *typesize ; /* 6++ */
  }
  *n = (int)( p - hdrbuf ) ;
  return(0) ;
}

