      Program Main
      Use FCIUtils
      Implicit None
!     Reads number of sites from the command line
!     This is a dummy routine to get around Fortran's restrictions
!     with dynamic memory allocation
      Integer*4 Nimp,No,Nstr,NS
!     Set number of excited states 
!     (1=ground only; 2=ground+one excited...)
      Real*8 t,U
      Character(128) :: Nimp_str,No_str,NS_str,p2_only_str
      Character(128) :: tdir,inname,outname,read_guess_str
      Logical :: read_guess,p2_only
!
!     Command line arguments:
!       1) Number of impurities
!       2) Number of states
!       3) Temporary folder path
!       4) Should we read the guess?

!     Set number of impurity sites from command line
      Call getarg(1, p2_only_str)
      Read(p2_only_str, *) p2_only
      Call getarg(2, Nimp_str)
      Read(Nimp_str, *) Nimp
      Call getarg(3, No_str)
      Read(No_str, *) No
!      Call getarg(4, NS_str)
!      Read(NS_str, *) NS

!     Find out where the input file is
      !Call getarg(3, tdir)
      !inname  = tdir(:len_trim(tdir)) // '/in.bin'
      !outname = tdir(:len_trim(tdir)) // '/out.bin'
      inname = 'in.bin'
      outname = 'out.bin'

!     Parse guess flag
!      Call getarg(5, read_guess_str)
!      Read(read_guess_str, *) read_guess

!     Nstr is the number of CI strings
      Nstr=NcR(Nimp,No)

!     Apparently file handles are global, so we'll open here to avoid
!     passing names. Yes, I know I'm the worst person.
      open(unit=80, file=inname, status='old', action='read', form='unformatted')
      open(unit=81, file=outname, status='replace', action='write', form='unformatted')
!      Call CIWrap(Nimp,No,NStr,NS,read_guess)

      If(p2_only) then
        Call P2Wrap(Nimp,No,Nstr)
      Else 
        Call getarg(4, NS_str)
        Read(NS_str, *) NS
!       Parse guess flag
        Call getarg(5, read_guess_str)
        Read(read_guess_str, *) read_guess
        Call CIWrap(Nimp,No,Nstr,NS,read_guess)
      End If

      End Program

      Subroutine P2Wrap(N,No,NStr)
      Use FCIUtils
      Implicit None
      Integer*4, Intent(In) :: N,No,NStr
      Integer*4 N2,N0,Max1,Max2
      Real*8 X(NStr,NStr),XH(NStr,NStr)
      Real*8 P2(N,N,N,N),T1(N,N),T2(N,N,N,N)
      Real*8 P1(N,N)
      Integer*4, Parameter :: N0Max=12
      read(80) X
      close(unit=80)

      N2=NcR(N,2); N0=Min(N0Max,NStr)
      Max1=NcR(N-1,No-1); Max2=NcR(N-2,No-2)
      Call Pdm2(N,No,N2,Max1,Max2,Nstr,X,T1,T2,XH,P1,P2)
      write(81) P1
      write(81) P2
      write(81) X
      close(unit=81)

      End Subroutine

      Subroutine CIWrap(N,No,NStr,NS,read_guess)
      Use FCIUtils
      Implicit None
      ! Read the big python data blob and fill the relevant arrays
      Integer*4, Intent(In) :: N,No,NStr,NS
      Logical, Intent(In) :: read_guess
      Integer*4 N2,N0,Max1,Max2,i,j,k,iS,jS,ij,c,l
      Real*8 h(N,N),V(N,N,N,N),hnat(N,N),Vnat(N,N,N,N),Xi(NStr,NStr,NS)
      Real*8 v1(No,No),V2(No,No,No,No)
      Real*8 P(N,N,NS,NS),P2(N,N,N,N),Ei(NS),Himp(NS,NS),Etmp
      Real*8 H1s(4,4),S1s(4,4)
!     N0Max controls size of preconditioner in FCI
!        Larger preconditioners make each step slower, but
!        fewer steps are (typically) required.
      Integer*4, Parameter :: N0Max=12


!     Initialize CI vectors
      Xi=0d0

!     Read input from python
      read(80) h
      read(80) V
      read(80) hnat
      read(80) Vnat
      if (read_guess) then
          read(80) Xi
      end if
      close(unit=80)

!     Set a Bunch of Info For FCI calculations
      N2=NcR(N,2); N0=Min(N0Max,NStr)
      Max1=NcR(N-1,No-1); Max2=NcR(N-2,No-2)

!H1s is the 1-site FCI hamlitonian that will be diagonalized in the python script
      P=0d0; P2=0d0; Ei=0d0; Himp=0d0; H1s=0d0; S1s = 0d0
      Call FCIman(N,No,N2,Nstr,Max1,Max2,N0,NS,h,V,hnat,Vnat,Xi,P,P2,Ei,H1s,S1s)

!!     Reconstruct Energy from P and OnTop (tests accuracy of P,OnTop)
!      Do iS=1,NS; Do jS=iS,NS
!         Etmp=0d0
!         Do i=1,N; Do j=1,N
!            Etmp=Etmp+P(i,j,iS,jS)*h(i,j)*2d0
!         End Do; End Do
!         Do i=1,N
!            Etmp=Etmp+OnTop(i,iS,jS)*V(i,i,i,i)
!         End Do
!         If(iS==jS) Write(6,*)'@#@#@Energy',iS,Etmp,Ei(iS)
!         If(iS.ne.jS) Write(6,*)'@#@#@Coupling',iS,Etmp
!      End Do; End Do

!     Construct Impurity Hamiltonian from P and OnTop
!     Note: Uses "central site" prescription for now 
!      ij=No/2; If((No/2)*2.ne.No)ij=ij+1
!      h(1:No,1:No)=h(1:No,1:No)-v1
!      V(1:No,1:No,1:No,1:No)=V(1:No,1:No,1:No,1:No)-V2
!      Himp=0d0
!      Do iS=1,NS; Do jS=iS,NS
!         Do i=1,N
!!     Note: factor of two for spin
!            Himp(iS,jS)=Himp(iS,jS)+P(i,ij,iS,jS)*h(i,ij)*2d0
!         End Do
!         Himp(iS,jS)=Himp(iS,jS)+OnTop(ij,iS,jS)*V(ij,ij,ij,ij)
!         Himp(jS,iS)=Himp(iS,jS)
!      End Do; End Do
!     Write results to python

      write(81) P
      write(81) P2
      write(81) Ei
      write(81) Xi
      write(81) H1s
!      write(81) S1s
      close(unit=81)
            
      End Subroutine
