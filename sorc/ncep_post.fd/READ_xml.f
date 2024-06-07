!> @file
!> @brief read_postxconfig() reads the post available field XML file and post control XML file. 
! Each set of output fields going to one output file will be saved and processed later. 
! In other words, post control file will be read in whole once. 
!> PROGRAM HISTORY LOG:
!>   01_27_2012  Jun Wang - INITIAL CODE
!>   03_10_2015  Lin Gan  - Replace XML file with flat file implementation with parameter marshalling
!>   07_08_2016 J. Carley - Clean up prints 

      SUBROUTINE READ_xml()

       use xml_perl_data,only: post_avblflds,paramset,read_postxconfig
       use grib2_module, only: num_pset
       use rqstfld_mod,only: num_post_afld,MXLVL,lvlsxml
       use CTLBLK_mod, only: me
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none

!>    READ post available field table

      call read_postxconfig()
      num_post_afld=size(paramset(1)%param)
      num_pset=size(paramset)

      RETURN
      end subroutine read_xml

