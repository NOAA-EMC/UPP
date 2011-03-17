#!/usr/bin/perl
#
# Build configuration file used during UPP compile command
# BE SURE TO RUN AS ./configure (to avoid getting a system configure command by mistake)
# 
# This file will ask the user which type of compile they would like to
# configure for - based on machine and available compilers
#
# After receiving valid user input the preamble file is placed in the configuration
# file, followed by the machine dependent compiler/linker/archive setting, followed
# by the postamble.  This will be the configure.upp file which is used to compile
# all of UPP or any subdirectory
#

#
# Initialize variables
$sw_netcdf_path = "" ;
$sw_usenetcdff = "" ;    # for 3.6.2 and greater, the fortran bindings might be in a separate lib file
$sw_wrf_path = "" ;
$sw_os = "ARCH" ;             # ARCH will match any
$sw_mach = "ARCH" ;           # ARCH will match any
$sw_fc = "\$(SFC)" ;
$sw_cc = "\$(SCC)" ;
$sw_f90 = "\$(SF90)" ;
$sw_dmparallel = "" ;
$sw_ompparallel = "" ;        # Not supported
$sw_comms_obj = "" ;
$sw_comms_objst = "" ;
$sw_comms_lib = "" ;
$sw_serial_mpi_stub = "" ;    # Assume parallel build
$sw_serial_mpi_lib = "" ;
$sw_bindir = "" ;   # bin directory 
$sw_incmod = "" ;   # include directory 
$sw_libdir = "" ;   # library directory


#
# Read in command line arguments :: set local variables
while ( substr( $ARGV[0], 0, 1 ) eq "-" )
 {
  if ( substr( $ARGV[0], 1, 7 ) eq "netcdf=" )
  {
    $sw_netcdf_path = substr( $ARGV[0], 8 ) ;
  }
  if ( substr( $ARGV[0], 1, 8 ) eq "wrfpath=" )
  {
    $sw_wrf_path = substr( $ARGV[0], 9 ) ;
  }
  if ( substr( $ARGV[0], 1, 3 ) eq "os=" )
  {
    $sw_os = substr( $ARGV[0], 4 ) ;
  }
  if ( substr( $ARGV[0], 1, 5 ) eq "mach=" )
  {
    $sw_mach = substr( $ARGV[0], 6 ) ;
  }
  if ( substr( $ARGV[0], 1, 11 ) eq "dmparallel=" )
  {
    $sw_dmparallel=substr( $ARGV[0], 12 ) ;
  }
  if ( substr( $ARGV[0], 1, 12 ) eq "ompparallel=" )
  {
    $sw_ompparallel=substr( $ARGV[0], 13 ) ;
  }
  if ( substr( $ARGV[0], 1, 11 ) eq "USENETCDFF=" )
  {
    $sw_usenetcdff = substr( $ARGV[0], 12 ) ;
  }
  if ( substr( $ARGV[0], 1, 7 ) eq "bindir=" )
  {
    $sw_bindir = substr( $ARGV[0], 8 ) ;
  }
  if ( substr( $ARGV[0], 1, 7 ) eq "incmod=" )
  {
    $sw_incmod = substr( $ARGV[0], 8 ) ;
  }
  if ( substr( $ARGV[0], 1, 7 ) eq "libdir=" )
  {
    $sw_libdir = substr( $ARGV[0], 8 ) ;
  }
  shift @ARGV ;
 }

#
#  Check for compression libraries needed in support of GRIB2 format
$JASPTERLIB = "";
$JASPERLIB  =  $ENV{JASPERLIB};
if ($JASPERLIB) {
  $sw_grib2_libs="-L$sw_wrf_path/external/io_grib2 -lio_grib2 -L$JASPERLIB -ljasper";
 print "grib2lib = $sw_grib2_libs";
}

#
# Check for HWRF environment variable set if applicable
if ( $ENV{HWRF}  ) {
  $sw_hwrf_libs="-L$sw_wrf_path/external/atm_ocn -latm_ocn";
} else {
  $sw_hwrf_libs="";
}

#
# Display the choices to the user and get selection
$validresponse = 0 ;

## UPP only supports serial @platforms = qw ( serial dmpar ) ;
@platforms = qw ( serial dmpar ) ;

until ( $validresponse ) {
  printf "------------------------------------------------------------------------\n" ;
  printf "Please select from among the following supported platforms.\n\n" ;

  open CONFIGURE_DEFAULTS, "< ./arch/configure.defaults" 
      or die "Cannot open ./arch/configure.defaults for reading" ;

#
# Read configure.defaults :: display all records which contain the ARCH directive,
# a matching OS, and a matching machine
  $opt = 1 ;
  while ( <CONFIGURE_DEFAULTS> )
  {
    for $paropt ( @platforms ) 
    {
      if ( substr( $_, 0, 5 ) eq "#ARCH" && ( index( $_, $sw_os ) >= 0 ) && ( index( $_, $sw_mach ) >= 0 )
	   && ( index($_, $paropt) >= 0 )  )
      {
        $optstr[$opt] = substr($_,6) ;
        $optstr[$opt] =~ s/^[ 	]*// ;
        $optstr[$opt] =~ s/#.*$//g ;
        chomp($optstr[$opt]) ;
        $optstr[$opt] = $optstr[$opt]." (".$paropt.")" ;
        if ( substr( $optstr[$opt], 0,4 ) ne "NULL" )
        {
          printf "  %2d.  %s\n",$opt,$optstr[$opt] ;
          $opt++ ;
        }
      }
    }
  }
  close CONFIGURE_DEFAULTS ;

#
# Get to end of our array
  $opt -- ;

#
# Get response - ask again and again - unless -1 -> exit
  printf "\nEnter selection [%d-%d] : ",1,$opt ;
  $response = <STDIN> ;

  if ( $response == -1 ) { exit ; }

  if ( $response >= 1 && $response <= $opt ) 
  { $validresponse = 1 ; }
  else
  { printf("\nInvalid response (%d)\n",$response);}
}
printf "------------------------------------------------------------------------\n" ;

$optchoice = $response ;

# 
# Open configure.defaults again to read record configuration settings
open CONFIGURE_DEFAULTS, "< ./arch/configure.defaults" 
      or die "Cannot open ./arch/configure.defaults for reading" ;
$latchon = 0 ;
while ( <CONFIGURE_DEFAULTS> )
{
  if ( substr( $_, 0, 5 ) eq "#ARCH" && $latchon == 1 )
  {
    $latchon = 0 ;
  }

#
# Got our record make substitutions with local variables set above
  if ( $latchon == 1 )
  {
    $_ =~ s/CONFIGURE_NETCDF_PATH/$sw_netcdf_path/g ;
    $_ =~ s/CONFIGURE_NETCDF_LIBS/$sw_usenetcdff -lnetcdf/g ;
    $_ =~ s/CONFIGURE_WRF_PATH/$sw_wrf_path/g ;
    $_ =~ s/CONFIGURE_FC/$sw_fc/g ;
    $_ =~ s/CONFIGURE_F90/$sw_f90/g ;
    $_ =~ s/CONFIGURE_CC/$sw_cc/g ;

    @machopts = ( @machopts, $_ ) ;
  }

#
# Loop through records to find the match based on ARCH directive and latchon
# matching OS, matching machine, matching processor selection
  for $paropt ( @platforms )
  {
    if ( substr( $_, 0, 5 ) eq "#ARCH" && $latchon == 0
          && ( index( $_, $sw_os ) >= 0 ) && ( index( $_, $sw_mach ) >= 0 )
          && ( index($_, $paropt) >= 0 ) )
    {
      $x=substr($_,6) ;
      $x=~s/^[     ]*// ;
      $x =~ s/#.*$//g ;
      chomp($x) ;
      $x = $x." (".$paropt.")" ;
      if ( $x eq $optstr[$optchoice] )
      {
        $latchon = 1 ;
        $sw_ompparallel = "" ;
        $sw_dmparallel = "" ;
        $validresponse = 0 ;

# Serial compile uses a stub library for mpi calls
        if ( $paropt eq 'serial' )
        {
          $sw_serial_mpi_stub  = "wrfmpi_stubs" ;
          $sw_serial_mpi_lib   = "-lmpi" ;
        }
# DM parallel
        elsif ( $paropt eq 'dmpar' ) 
        {
          $sw_comms_lib = "" ;
          $sw_comms_obj = "" ;
          $sw_comms_objst = "";
          $sw_dmparallel = "" ;
          $sw_dmparallelflag = "-DDM_PARALLEL" ;
          $sw_fc = "\$(DM_FC)" ;
          $sw_f90 = "\$(DM_F90)" ;
          $sw_cc = "\$(DM_CC)" ;
        }
      }
    }
  }
}
close CONFIGURE_DEFAULTS ;

#
# Build configure.upp
open CONFIGURE_UPP, "> configure.upp" or die "cannot append configure.upp" ;

#
# preamble
open ARCH_PREAMBLE, "< arch/preamble" or die "cannot open arch/preamble" ;
my @preamble;
#
# apply substitutions to the preamble...
while ( <ARCH_PREAMBLE> )
  {
  @preamble = ( @preamble, $_ ) ;
  }
close ARCH_PREAMBLE ;

print CONFIGURE_UPP @preamble  ;
close ARCH_PREAMBLE ;

#
# machine/compiler configuration values
printf CONFIGURE_UPP "# Settings for %s", $optstr[$optchoice] ;
print CONFIGURE_UPP @machopts  ;

#
# postamble
open ARCH_POSTAMBLE, "< arch/postamble" or die "cannot open arch/postamble" ;
while ( <ARCH_POSTAMBLE> ) { 
    $_ =~ s/CONFIGURE_NETCDF_PATH/$sw_netcdf_path/g ;
    $_ =~ s/CONFIGURE_NETCDF_LIBS/$sw_usenetcdff -lnetcdf/g ;
    $_ =~ s/CONFIGURE_WRF_PATH/$sw_wrf_path/g ;
    $_ =~ s/CONFIGURE_COMMS_OBJST/$sw_comms_objst/g ;
    $_ =~ s/CONFIGURE_COMMS_OBJ/$sw_comms_obj/g ;
    $_ =~ s/CONFIGURE_COMMS_LIB/$sw_comms_lib/g ;
    $_ =~ s/CONFIGURE_GRIB2_LIBS/$sw_grib2_libs/g ;
    $_ =~ s/CONFIGURE_HWRF_LIBS/$sw_hwrf_libs/g ;
    $_ =~ s/CONFIGURE_SERIAL_MPI_STUB/$sw_serial_mpi_stub/g ;
    $_ =~ s/CONFIGURE_SERIAL_MPI_LIB/$sw_serial_mpi_lib/g ;
    $_ =~ s/CONFIGURE_BLD_BINDIR/$sw_bindir/g ;
    $_ =~ s/CONFIGURE_BLD_INCMOD/$sw_incmod/g ;
    $_ =~ s/CONFIGURE_BLD_LIBDIR/$sw_libdir/g ;
  print CONFIGURE_UPP;
 }
close ARCH_POSTAMBLE ;
close CONFIGURE_UPP ;

printf "Configuration successful. To build the UPP, type: compile \n" ;
printf "------------------------------------------------------------------------\n" ;


