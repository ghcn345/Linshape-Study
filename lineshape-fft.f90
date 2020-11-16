module type_kinds
  implicit none

  Private

  !Integer Kinds
  Public :: byte,short,long

  !Floating point Kinds
  Public :: single,double

  !String Lengths
  Public :: title,path,line,comment 

  !> string length for a name/title
  integer, parameter:: title=50 
  !> string length for a directory path
  integer, parameter:: path=100
  !> string length for a line in a file
  integer, parameter:: line=100
  !> string length for a message/comment
  integer, parameter:: comment=500

  !> byte integer size 
  integer, parameter:: byte=selected_int_kind(1)
  !> short integer size 
  integer, parameter:: short=selected_int_kind(4)
  !> long integer size 
  integer, parameter:: long=selected_int_kind(8)

  !> single precision float size 
  integer, parameter:: single=selected_real_kind(6)
  !> double precision float size 
  integer, parameter:: double=selected_real_kind(15)

end module type_kinds
!===========================================================
!===========================================================
!> \brief
!! Timming Module
!! \details
!! Supplies the software with a standard method of computing
!! wall clock time in units of mils, sec, hours or even days.
!! \authors
!! Author: Daniel Montemayor
!<===========================================================
module wallclock
  use type_kinds
  implicit none
  
  private
  public:: date,time,zone,timearray
  
     character(len=8)::d
     character(len=10)::t
     character(len=5)::z
     integer(long) :: a(8)
  
contains
  !----------------------------
  function date()!<returns date in 8 char string format YYYYMMDD
    character(len=8)::date
    call date_and_time(d,t,z,a)
    date=d
  end function date
  !----------------------------
  function time()!<returns time in 10 char string format hhmmss.sss
    character(len=10)::time
    call date_and_time(d,t,z,a)
    time=t
  end function time
  !----------------------------
  function zone()!<returns zone relative to Coordinated Univeral Time (UTC) in 5 char string format (+-)hhmm
    character(len=5)::zone
    call date_and_time(d,t,z,a)
    zone=z
  end function zone
  !----------------------------
  function timearray()!<returns 8 element integer array representing time in format Yr,Mo,Day,UTCmin_dif,hr,min,sec,msec
    integer(long) :: timearray(8)
    call date_and_time(d,t,z,a)
    timearray=a
  end function timearray

end module wallclock
!===========================================================
!> \brief
!! Enables MPI
!! \details
!! Checks if MPI is available from preprocessor and sets up appropriate variables
!<
module MPIframework
  use type_kinds
  implicit none
  
!!$# ifdef MPI
!!$  include "mpif.h"
!!$  integer(long)::status(MPI_STATUS_SIZE)
!!$# endif
  
  integer(long), parameter :: master=0
  integer(long)::myid=0
  integer(long)::nproc=1
  
end module MPIFRAMEWORK
!===========================================================
!===========================================================
module rand
!> \brief
!! Random number generators
!<

  use type_kinds
  use wallclock
  use MPIframework
  implicit none

  private
  public :: seed,setseed,ran0,gran

  integer(long) :: idum=0_long

contains
  subroutine setseed(seedin)
    IMPLICIT NONE
    integer(long) :: seedin
    idum=seedin
    RETURN
  END subroutine setseed
  !----------------------------
  function seed()
    IMPLICIT NONE
    integer(long) :: seed
    seed=idum
    RETURN
  END function seed
  !----------------------------
  FUNCTION RAN0(seedin)
    IMPLICIT NONE
    REAL(double):: RAN0
    integer(long), optional:: seedin
    INTEGER(long),PARAMETER:: IA=16807_long
    INTEGER(long),PARAMETER:: IM=2147483647_long
    INTEGER(long),PARAMETER:: IQ=127773_long
    INTEGER(long),PARAMETER:: IR=2836_long
    INTEGER(long),PARAMETER:: MASK=20410876_long
    REAL(double),PARAMETER:: AM=1.0_double/IM
    INTEGER(long):: k

    integer(long)::array(8)

    if(present(seedin))idum=seedin
    if(idum.EQ.0)then
       do idum=0, myid*100000!E5
          array=timearray()
       end do
       idum=myid*1E4
       idum=idum+array(8)
       idum=idum+array(7)*1E3
       idum=idum+array(6)*6E4
       idum=idum+array(5)*36E5
       !write(*,*)"core ",myid,"generating new random seed ",idum)
    end if

    idum=IEOR(idum,MASK)
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF (idum.LT.0) idum=idum+IM
    RAN0=AM*idum
    idum=IEOR(idum,MASK)
    RETURN
  END FUNCTION RAN0
  !----------------------------
  FUNCTION GRAN(seed)
    !         FUNCTION TO GENERATE GAUSSIAN RANDOM NUMBERS
    !         USING BOX-MUELLER/WEINER METHOD.
    IMPLICIT NONE
    integer(long), optional:: seed
    REAL(double) :: GRAN,w1,w2
    real(double) , parameter :: PI=3.1415926535898_double

    if(present(seed))idum=seed

    w1=RAN0()
    w2=RAN0()
    gran=SQRT(-2.*LOG(w1))*COS(2.*PI*w2)

  END FUNCTION GRAN

end module rand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module gl




end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE four1(data,nn,isign)
implicit none

INTEGER isign,nn
REAL*8 data(2*nn)
INTEGER i,istep,j,m,mmax,n
REAL*8 tempi,tempr
real*8 theta,wi,wpi,wpr,wr,wtemp

n=2*nn
j=1
do i=1,n,2
   if(j.gt.i)then
      tempr=data(j)
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
   end if
   m=n/2
1  if ((m.ge.2).and.(j.gt.m)) then
      j=j-m
      m=m/2
      goto 1
   end if
   j=j+m
end do
mmax=2
2 if (n.gt.mmax) then
   istep=2*mmax
   theta=6.28318530717959d0/(isign*mmax)
   wpr=-2.d0*sin(0.5d0*theta)**2
   wpi=sin(theta)
   wr=1.d0
   wi=0.d0
   do m=1,mmax,2
      do i=m,n,istep
         j=i+mmax
         tempr=wr*data(j)-wi*data(j+1)
         tempi=wr*data(j+1)+wi*data(j)
         data(j)=data(i)-tempr
         data(j+1)=data(i+1)-tempi
         data(i)=data(i)+tempr
         data(i+1)=data(i+1)+tempi
      end do
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
   end do
   mmax=istep
   goto 2
end if
return 
END SUBROUTINE four1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swap(a,b)
implicit none

real*8::a,b,temp

temp=a
a=b
b=temp

end subroutine swap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program main
use rand
implicit none

integer::ki
real*8::j(0:1099)=0.0d0
real*8::w
integer::ti
integer::m

real*8::re
real*8::im
real*8::wi
real*8::in

real*8:: b=529.610525d0     ! kT=0.00094409
real*8:: wb=0.006835371d0
real*8:: n=0.29d0
real*8:: r,dt,t
integer::mi,mn,l,fi,max=8192
real*8::data(0:32768*2-1)=0.0d0

r=0.2d0*wb

write (*,*) 'Starting program'

 do m=1,600
         
         w=m*0.00005d0
         j(m-1)=n*w*wb*wb*wb*wb/((wb*wb-w*w)*(wb*wb-w*w)+4*w*w*r*r)
   	
      end do

      !  in=0.0d0
     ! dt=20.0d0
      dt=1.0d0
   do ti=-max/2,max/2-1 
      t=ti*dt
      re=0.d0
      im=0.d0
      
      do m=1,600         
         w=m*0.00005d0
                 
  	 re=j(m-1)*(1-cos(w*t))*0.00005d0/(w*w*tanh(b*w))+re
         im=j(m-1)*(sin(w*t)-w*t)*0.00005d0/(w*w)+im
         write(22,*)w,j(m-1),re,im
      end do
      stop
    !  write(33,*)t,in
    ! in=in+exp(-re)*cos(im+wi*t)
     if(ti<0)then
         mn=ti+max
         data(mn*2)=dt*exp(-re)*cos(im)
         data(mn*2+1)=-dt*exp(-re)*sin(im)
        
      else
         data(ti*2)=dt*exp(-re)*cos(im)
         data(ti*2+1)=-dt*exp(-re)*sin(im)
     
      end if
    end do
    call four1(data,max,-1)

    
      do l=0,max-1
         call swap(data(l),data(l+max))
      end do

 ! ki=-300,300
 ! wi=ki*0.0001d0         
    
    do ki=3000,5000
       wi=6.28318530717959*(ki-max/2)/(max*dt)
       
       write(11,*) wi,data(ki*2)
       end do
write(*,*)'Program Done!'	
	

end Program main
