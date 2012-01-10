      subroutine TRANSPORTWIND(im,jm,lm,jsta,jend,jsta_2l,jend_2u &
      ,mixht,u10,v10,uh,vh,zmid,pint,utran,vtran,modelname) 
!      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,jpds,jgds,kpds,buf)
                                                                                
      implicit none
      INCLUDE "mpif.h"
!      INCLUDE "CTLBLK.comm"
      character(len=4) :: modelname
      integer im,jm,lm,jsta,jend,jsta_2l,jend_2u
      integer i,j,l,ie,iw
      real dp,mixhtv,pv1,pv2,zmidv
      real mixht(im,jsta_2l:jend_2u),u10(im,jsta_2l:jend_2u) &
           ,v10(im,jsta_2l:jend_2u)      
      real utran(im,jsta:jend),vtran(im,jsta:jend),psum(im,jsta:jend)
      real uh(im,jsta_2l:jend_2u,lm),vh(im,jsta_2l:jend_2u,lm) &
          ,zmid(im,jsta_2l:jend_2u,lm),pint(im,jsta_2l:jend_2u,lm+1)
                                                                                
      utran=0.
      vtran=0.
      psum=0.
      
      IF(MODELNAME == 'NMM')THEN
        do l=1,lm
          do j=jsta,jend
            do i=1,im
	      ie=i+mod(j,2)
	      iw=i+mod(j,2)-1
	      if((j==1 .or. j==jm) .and. i < im)then   !s & n bc
                mixhtv=0.5*(mixht(i,j)+mixht(i+1,j))
		zmidv=0.5*(zmid(i,j,l)+zmid(i+1,j,l))
              else if((i==1 .or. i==im) .and. mod(j,2) == 0) then   !w & e even bc
                mixhtv=0.5*(mixht(i,j-1)+mixht(i,j+1))
		zmidv=0.5*(zmid(i,j-1,l)+zmid(i,j+1,l))
              else
                mixhtv=0.25*(mixht(iw,j)+mixht(ie,j) &
                +mixht(i,j+1)+mixht(i,j-1))
                zmidv=0.25*(zmid(iw,j,l)+zmid(ie,j,l) &
                +zmid(i,j+1,l)+zmid(i,j-1,l))
              end if
	      if(zmidv<=mixhtv)then
	        if((j==1 .or. j==jm) .and. i < im)then   !s & n bc
                  pv1=0.5*(pint(i,j,l)+pint(i+1,j,l))
		  pv2=0.5*(pint(i,j,l+1)+pint(i+1,j,l+1))
                else if((i==1 .or. i==im) .and. mod(j,2) < 0) then   !w & e even bc
                  pv1=0.5*(pint(i,j-1,l)+pint(i,j+1,l))
		  pv2=0.5*(pint(i,j-1,l+1)+pint(i,j+1,l+1))
                else
                  pv1=0.25*(pint(iw,j,l)+pint(ie,j,l) &
                  +pint(i,j+1,l)+pint(i,j-1,l))
                  pv2=0.25*(pint(iw,j,l+1)+pint(ie,j,l+1) &
                  +pint(i,j+1,l+1)+pint(i,j-1,l+1))
                end if
	        dp=pv2-pv1
	        utran(i,j)=utran(i,j)+uh(i,j,l)*dp
	        vtran(i,j)=vtran(i,j)+vh(i,j,l)*dp
	        psum(i,j)=psum(i,j)+dp
	      end if
	    end do
	  end do
        end do 
      else    
        do l=1,lm
          do j=jsta,jend
            do i=1,im
	      if(zmid(i,j,l)<=mixht(i,j))then
	        dp=pint(i,j,l+1)-pint(i,j,l)
	        utran(i,j)=utran(i,j)+uh(i,j,l)*dp
	        vtran(i,j)=vtran(i,j)+vh(i,j,l)*dp
	        psum(i,j)=psum(i,j)+dp
	      end if
	    end do
	  end do
        end do	      
      end if	      
! divide by psum
      do j=jsta,jend
        do i=1,im
	  if(psum(i,j)>0.)then
	    utran(i,j)=utran(i,j)/psum(i,j)
	    vtran(i,j)=vtran(i,j)/psum(i,j)
	  else
	    utran(i,j)=u10(i,j) ! if no mix layer, specify 10 m wind, per DiMego,
	    vtran(i,j)=v10(i,j)
	  end if
	end do
      end do 	            	  
      print*, 'sample transport wind= ',utran(im/2,jsta),vtran(im/2,jsta)                                                                          
      RETURN
      END
