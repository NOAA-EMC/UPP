#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*
**  Define the maximum length of a mnemonic in GRIB2 Code Table 4.2
*/
#define MXG2MNEM 16
#define MXG2MNEMP4 ( MXG2MNEM + 4 )

/*
**  Allocate internal memory for table entries in chunks of 1000
*/
#define NUMALLOC 1000

/*
**  On certain operating systems, the FORTRAN compiler appends an underscore
**  to subprogram names in its object namespace.  Therefore, on such systems,
**  a matching underscore must be appended to any C language references to the
**  same subprogram names so that the linker can correctly resolve such
**  references across the C <-> FORTRAN interface at link time.
*/
#ifdef LINUX
#define open_and_read_4dot2 open_and_read_4dot2_ 
#define search_for_4dot2_entry search_for_4dot2_entry_
#define close_4dot2 close_4dot2_
#endif

/*
**  In order to ensure that the C <-> FORTRAN interface works properly (and
**  portably!), the default size of an "INTEGER" declared in FORTRAN must be
**  identical to that of an "int" declared in C.  If this is not the case (e.g.
**  some FORTRAN compilers, most notably AIX via the -qintsize= option, allow the
**  sizes of INTEGERs to be definitively prescribed outside of the source code
**  itself!), then the following conditional directive (or a variant of it) can
**  be used to ensure that the size of an "int" in C remains identical to that
**  of an "INTEGER" in FORTRAN.
*/
#ifdef F77_INTSIZE_8
    typedef long f77int;
#else
    typedef int f77int;
#endif

/*
**  Define global variables
*/
struct TableEntry {
    char mnemonic[MXG2MNEMP4];
    int discipline;
    int category;
    int parameter;
};

struct TableEntry *pe0 = NULL;

size_t nentry = 0;

/*
**  Declare function prototypes for ANSI C compatibility
*/
void open_and_read_4dot2( char *, f77int * );
void search_for_4dot2_entry( char *, f77int *, f77int *, f77int *, f77int *, f77int * );
int compar( const struct TableEntry *, const struct TableEntry * );
void close_4dot2( f77int * );

/*
**  Define the internal comparison function for use with qsort and bsearch.
*/
int compar( const struct TableEntry *pte1, const struct TableEntry *pte2 )
{
    return strcmp( pte1->mnemonic, pte2->mnemonic );
}

/*
**  Define the remaining subprograms which will be called by user applications.
*/

/********************************************************************************
 * open_and_read_4dot2								*
 *										*
 * Opens and reads the GRIB2 Code Table 4.2 into an internal memory structure.	*
 *										*
 * Usage:									*
 *      call open_and_read_4dot2 ( filename, iret )				*
 *										*
 * Input parameters:								*
 *	filename	character*(*)	Location of table file on filesystem;	*
 *					directory prefixes or other local	*
 *					filesystem notation is allowed up to	*
 *					120 total characters.			*
 *										*
 * Output parameters:								*
 *	iret		integer		Return code:				*
 *					 0 = normal return			*
 *					-1 = input filename was more than	*
 *					     120 characters			*
 *					-2 = table file could not be opened	*
 *					-3 = memory allocation error		*
 **										*
 *  Log:				                                	*
 *  J. Ator/NCEP	02/10							*
 *******************************************************************************/
void open_and_read_4dot2( char *filename, f77int *iret )
{
#define MXFNLEN 120

    struct TableEntry *pra;

    char lfn[MXFNLEN+1], str[81], cflag, cub = '_';

    size_t i;

    FILE *pfn;

/*
**  Copy the input filename into a local variable and check it for validity.
**  This is especially important in case the filename was passed in as a
**  string literal by the calling program or else doesn't have a trailing
**  NULL character.
*/
    for ( i = 0; ( ! isspace( filename[i] ) && ! iscntrl( filename[i] ) ); i++ ) {
	if ( i == MXFNLEN ) {
	    *iret = ( f77int ) -1;
	    return;
	}
	lfn[i] = filename[i];
    }
    lfn[i] = '\0';

/*
**  Open the file.
*/
    if ( ( pfn = fopen( lfn, "r" ) ) == NULL ) { 
	*iret = ( f77int ) -2;
        printf("cant open file,%s\n",lfn);
	return;
    }

/*
**  Read the file contents into an internal memory structure.  Memory will be
**  allocated as needed for NUMALLOC entries at a time.
*/
    while ( fgets( str, 80, pfn ) != NULL ) {
	if ( str[0] != '!' ) {    /* ignore comment lines */
	    if ( ( nentry % NUMALLOC ) == 0 ) {
/*
**		Allocate additional memory.
*/
		pra = realloc( pe0, ( NUMALLOC * sizeof( struct TableEntry ) ) );
		if ( pra == NULL ) {
		    *iret = ( f77int ) -3;
		    return;
		}
		pe0 = pra;
	    }
	    sscanf( str, "%d%d%d%*3c%c%*c%s",
		    &pe0[nentry].discipline, &pe0[nentry].category,
		    &pe0[nentry].parameter, &cflag,
		    pe0[nentry].mnemonic );
	    strncat( pe0[nentry].mnemonic, &cub, 1 );
	    strncat( pe0[nentry].mnemonic, &cflag, 1 );
	    nentry++;
	}
    }

/*
**  Sort the entries within the internal memory structure.
*/
    qsort( pe0, nentry, sizeof( struct TableEntry ),
		( int (*) ( const void *, const void * ) ) compar );

/*
**  Close the file.
*/
    fclose ( pfn );   

    *iret = ( f77int ) 0;
}

/********************************************************************************
 * search_for_4dot2_entry							*
 *										*
 * Searches for a specified mnemonic within the previously-opened GRIB2 Code	*
 * Table 4.2 and returns the corresponding product discipline, parameter	*
 * category and parameter number.  A binary search algorithm is used.		*
 *										*
 * Usage:									*
 *      call search_for_4dot2_entry ( nemo, locflg, disc, catg, parm, iret )	*
 *										*
 * Input parameters:								*
 *	nemo		character*(*)	Mnemonic (of up to 16 characters in	*
 *					length) to search for within table.	*
 *	locflg		integer		Version of mnemonic to be returned, in	*
 *					case of duplication within table:	*
 *					 0 = international version (default)	*
 *					 1 = local version			*
 *										*
 * Output parameters:								*
 *	disc		integer		Product discipline			*
 *	catg		integer		Parameter category			*
 *	parm		integer		Parameter number			*
 *	iret		integer		Return code:				*
 *					 0 = normal return			*
 *					-1 = nemo not found within table for	*
 *					     specified locflg version		*
 **										*
 *  Log:				                                	*
 *  J. Ator/NCEP	02/10							*
 *******************************************************************************/
void search_for_4dot2_entry( char nemo[MXG2MNEM], f77int *locflg,
			     f77int *disc, f77int *catg, f77int *parm,
			     f77int *iret )
{
  /* chg from short to long because of a PGI compiler error */
    unsigned long n = 0, n2 = 0;

    long llf;

    struct TableEntry key, *pbs;

    size_t ipt;

/*
**  Make a local copy of nemo.  Mnemonics may consist of any combination of
**  alphanumeric, underscore and dash characters.
*/
    while (  ( n < MXG2MNEM ) &&
		( ( isalnum( ( int ) nemo[n] ) ||
		  ( nemo[n] == '_' ) || ( nemo[n] == '-' ) ) )  ) {
	    key.mnemonic[n2++] = nemo[n++];
    }

/*
**  Append an underscore followed by the locflg in order to generate the
**  mnemonic to search for. 
*/
    key.mnemonic[n2++] = '_';
    llf = ( long ) *locflg;
    if ( llf != 1 ) llf = 0;   /* default to using international entry unless
				  local is specified */
    sprintf( &(key.mnemonic[n2]), "%ld", llf );  /* trailing null will be automatically
						  appended by sprintf */

/*
**  Search for the mnemonic in the Code Table and return appropriate output values.
*/
    pbs = bsearch( &key, pe0, nentry, sizeof( struct TableEntry ),
			( int (*) ( const void *, const void * ) ) compar );
    if ( pbs == NULL ) {
	*iret = ( f77int ) -1;
    }
    else {
	*iret = ( f77int ) 0;
	ipt = pbs - pe0;
	*disc = ( f77int ) pe0[ipt].discipline;
	*catg = ( f77int ) pe0[ipt].category;
	*parm = ( f77int ) pe0[ipt].parameter;
    }
	
}

/********************************************************************************
 * close_4dot2									*
 *										*
 * This subroutine should be called one time at the end of the application	*
 * program in order to free all allocated memory.				*
 *										*
 * Usage:									*
 *      call close_4dot2 ( iret )						*
 *										*
 * Output parameters:								*
 *	iret		integer		Return code:				*
 *					 0 = normal return			*
 **										*
 *  Log:				                                	*
 *  J. Ator/NCEP	02/10							*
 *******************************************************************************/
void close_4dot2( f77int *iret )
{
    free ( pe0 );

    *iret = ( f77int ) 0;
}
