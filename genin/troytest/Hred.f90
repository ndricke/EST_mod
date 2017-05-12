      Program Main
      Use FCIUtils
      Implicit None
!     Reads number of sites from the command line
!     This is a dummy routine to get around Fortran's restrictions
!     with dynamic memory allocation
      Integer*4 Nimp,No,NStr
      Character(128) :: Nimp_str,No_str

!     Set number of impurity sites from command line
      Call getarg(1, Nimp_str)
      Read(Nimp_str, *) Nimp
      Call getarg(2, No_str)
      Read(No_str, *) No

      NStr=NcR(Nimp,No)
!     Find out where the input file is
      !Call getarg(3, tdir)
      !inname  = tdir(:len_trim(tdir)) // '/in.bin'
      !outname = tdir(:len_trim(tdir)) // '/out.bin'
!      write(*,*) NStr
      
      Call HredWrap(Nimp,No,NStr)

      End Program

!I'm calling Nimp as N here, as this could be more general
      Subroutine HredWrap(N,No,NStr) 
      Use FCIUtils
!      Integer*4 Intent(In) :: N,No
      Integer*4 N,No
      Integer*4 NStr,N0,N2,Max1,Max2
      Real*8 h(N,N),V(N,N,N,N)
      Real*8 X1(NStr,NStr),X2(NStr,NStr),X3(NStr,NStr),X4(NStr,NStr)
      Real*8 H1s(4,4)
      Integer*4, Parameter :: N0Max=12
      Character(128) :: inname,outname

      Max1=NcR(N-1,No-1); Max2=NcR(N-2,No-2)
      N2=NcR(N,2); N0=Min(N0Max,NStr)

      inname = 'in.bin'
      outname = 'out.bin'
      open(unit=80, file=inname, status='old', action='read', form='unformatted')

      read(80) h
      read(80) V
      read(80) X1
      read(80) X2
      read(80) X3
      read(80) X4
      close(unit=80)

      Call Hred(N,No,N2,Max1,Max2,Nstr,X1,X2,X3,X4,h,V,H1s)

      open(unit=81, file=outname, status='replace', action='write', form='unformatted')
!      write(81) NStr
!      write(81) X1
      write(81) H1s
!      write(81) h
!      write(81) V
      close(unit=81)
      End Subroutine

