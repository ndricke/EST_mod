      Subroutine FCIman(N,No,N2,Nstr,Max1,Max2,N0,NS,h,V,hnat,Vnat,
     $     Xi,P,P2,Ei,H1s,S1s)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     N     - Number of Basis Functions
!     No    - Number of Alpha electrons (=Number of Beta Electrons)
!     N2    - N*(N-1)/2
!     Nstr  - Number of Alpha CI Strings (=Number of Beta Strings)
!     Max1  - Max. # of states a 1 particle operator can connect to
!     Max2  - Max. # of states a 2 particle operator can connect to
!     N0    - Number of strings in preconditioner space
!     NS    - Number of eigenstates desired
!
!     h     - One electron Hamilontian
!     V     - Two electron Hamiltonian
!
!     Xi    - Input: NS Approximate guess eigenvectors
!             Output: Lowest NS Eigenvectors
!     P     - One particle alpha transition density matrix
!     P2    - 2 particle density matrix. Currently only for ground state calc
!     Ei    - Lowest NS eigenvalues
!     H1s   - Full hamiltonian for 1-site embedding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Nstr,Max1,Max2,NS,N0,MyMod,iSym,iX,iS,jS,iSpin
      Integer*4 i,j,k,l,m,a,b,ij,kl,ab,ii,jj,kk,ll,ierr,iter,Info
      Integer*4 Istr(N0)
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 Iocc(No),Isubst(No),Inda,Indb,Isigna,Isignb
      Real*8 Energy,q,dum,Uncertainty(NS)
      Real*8 Norm,DE,EX,EY,NormX,NormY,fac,DeltaE,scale
      Real*8 X(Nstr,Nstr),X1(Nstr,Nstr),XH(Nstr,Nstr),X1H(Nstr,Nstr)
      Real*8 Xi(Nstr,Nstr,NS),Ei(NS),Sym(NS)
      Real*8 Xn1(Nstr,Nstr),Xn2(Nstr,Nstr)
      Real*8 Xn2t(Nstr,Nstr),Xn2d(Nstr,Nstr)
      Real*8 H0(N0*N0,N0*N0),E0(N0*N0),U0(N0*N0,N0*N0)
      Real*8 Hd(Nstr,Nstr)
      Real*8 Hm(2,2),eig(2),U(2,2),H1s(4,4),S1s(4,4),eps
      Real*8 h(N,N),hnat(N,N),V(N,N,N,N),Vnat(N,N,N,N)
      Real*8 P(N,N,NS,NS),T1(N,N),T2(N,N,N,N),P2(N,N,N,N)
      Real*8 X1a(Nstr,Nstr),X2(Nstr,Nstr),X3(Nstr,Nstr)
      Real*8 X1Ha(Nstr,Nstr),X2H(Nstr,Nstr),X3H(Nstr,Nstr)
      Real*8 X4(Nstr,Nstr),X4H(Nstr,Nstr)
      Common /ZIndex/ Zptr
      Integer*4,Pointer :: Zptr(:,:)
      Integer*4,Target :: Zindex(N,N)
      Real*8 zero,one,two,three,four
      Real*8, Parameter :: Thresh=1.0d-9
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      Energy=zero; Zptr=>Zindex; Fac=One;
!     Build indexing array for future use
      ZIndex=0
      Do k=1,No; Do l=k,N-No+k 
         If(k==No) Then
            ZIndex(k,l)=l-k
         Else
            Do m=N-l+1,N-k
               ZIndex(k,l)=ZIndex(k,l)+NcR(m,No-k)-NcR(m-1,No-k-1)
            End Do
         End If
      End Do; End Do

!!----Determine which strings are connected by various operators.---!!
      Call IString(N,No,N2,Max1,Max2,Ex1,Ex2)
      
!     Build Diagonal part of H
      Hd=Zero
      Call GetHd(N,No,Nstr,h,V,Iocc,Isubst,Hd,1)

!     Get N0 lowest Strings
      Call GetIstr(N,No,Nstr,N0,Hd,Istr)
      Call Isort(N0,Istr,Info)

!     Build + Diagonalize H0
      Call GetH0(N,No,N0,N2,Max1,Max2,Nstr,Ex1,Ex2,Istr,h,V,H0)
      Call Diagonalize(N0*N0,H0,E0,U0,Info)

!     Big Loop over states
      Ei=Zero; iX=0
      Do iS=1,NS
!     Initial vector for this state (ensure it is a singlet for now)
!     This restarts from the input Xi vector if Xi is nonzero
      X=Xi(:,:,iS); iSpin=1; Norm=MDot(X,X) 
      If(Norm.lt.0.01d0) Then
         iSpin=-1
         Do While(iSpin==-1)
            iX=iX+1
            X=Zero; ij=0
            Do i=1,N0; Do j=1,N0; ij=ij+1
               X(Istr(j),Istr(i))=U0(ij,iX)
            End Do; End Do
            X1=X-Transpose(X)
            Norm=MDot(X1,X1)/MDot(X,X)
            If(Norm<1d-2) iSpin=1
         End Do
      End If
 
!     Check if the state has even (+1) or odd (-1) S
!     Even is singlet, quintet...  | Odd is triplet, sestet...
      X1=X-Transpose(X)
      If(MDot(X1,X1).le.1d-1) iSym=1
      If(MDot(X1,X1).gt.1d-1) iSym=-1
      Sym(iS)=iSym

!     Fix up broken spin symmetry (if it exists)
      X=(X+iSym*Transpose(X))/Two; X=X/Sqrt(MDot(X,X))

!     Make initial vector orthogonal to lower eigenvectors
      Call GS(X,Xi,iS)

!     Compute Initial guess energy
      Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)
      Energy=MDot(X,XH); DE=Energy/Two; iter=0
!      Write(6,*)'Guess:',E0(iS),Energy

!     Iteratively determine eigenvalue #iS of the Hamiltonian
      Fac=1.0d0
      Do While(Abs(fac)>Thresh.and.iter<5000)
         iter=iter+1; DE=Energy
!     Make the (orthogonal component of the) Davidson Update
         X1=-(XH-Energy*X)/(Hd-Energy)

!     Build (H0-Energy)^-1 (excluding eigenvalues that might blow up)
         H0=Zero
         Do k=1,N0*N0
            fac=0d0 
            If(abs(E0(k)-Energy).gt.1d-2) Then
               fac=1d0/(E0(k)-Energy)
            End If
            Do i=1,N0*N0; Do j=1,N0*N0
               H0(i,j)=H0(i,j)+U0(i,k)*U0(j,k)*fac
            End Do; End Do
         End Do

!     Build the Davidson Update using H0
         ij=0
         Do i=1,N0; Do j=1,N0; ij=ij+1; kl=0 
            X1(Istr(j),Istr(i))=Zero
            Do k=1,N0; Do l=1,N0; kl=kl+1
               X1(Istr(j),Istr(i))= X1(Istr(j),Istr(i))
     $           -(XH(Istr(l),Istr(k))-Energy*X(Istr(l),Istr(k)))
     $           *H0(ij,kl)
         End Do; End Do; End Do; End Do      

!     Apply Spin Symmetry and Gramm-Schmidt to Update
         X1=(X1+iSym*Transpose(X1))/Two; X1=X1/Sqrt(MDot(X1,X1))
         X1=X1-MDot(X1,X)*X; X1=X1/Sqrt(MDot(X1,X1))
         X1=X1-MDot(X1,X)*X; X1=X1/Sqrt(MDot(X1,X1))

!     Correct if Davidson has given us a bad vector
!        X1 should not be orthogonal to (H-E)X
!        If it is (nearly) orthogonal, add a bit of (H-E)X to X1
!        If we don't do this, it occasionally gives false convergence
!        This should happen rarely and eps controls how often 
!         this correction is invoked
!        A more elegant fix might be to just have a better
!         preconditioner.
         X1H=XH-Energy*X; Call GS(X1H,Xi,iS);
         X1H=X1H/Sqrt(MDot(X1H,X1H))
         eps=0.1d0; fac=abs(MDot(X1,X1H))
         If(fac.lt.eps) Then
            X1=X1+2d0*eps*X1H
            X1=X1-MDot(X1,X)*X; X1=X1/Sqrt(MDot(X1,X1))
            X1=X1-MDot(X1,X)*X; X1=X1/Sqrt(MDot(X1,X1))
         End If

!     Make X1 orthogonal to lower eigenvectors
         Call GS(X1,Xi,iS)

!     Act H on Davidson Vector
         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X1,h,V,X1H)

!     Build Hm
         Hm(1,1)=MDot(X,XH)
         Hm(2,2)=MDot(X1,X1H)
         Hm(1,2)=MDot(X,X1H); Hm(2,1)=Hm(1,2)

!     Diagonalize Hm
         Call GetEig2(Hm,Eig,U)

!     Keep Lowest Eigenvector
         If(iter<50.or.(iter/10)*10/=iter) Then
            fac=U(2,1)
            X=U(1,1)*X+U(2,1)*X1
            XH=U(1,1)*XH+U(2,1)*X1H
         Else
!     If convergence is slow, sometimes it is because we are
!     Swiching back and forth between two near-solutions. To
!     fix that, every once in a while, we take half the step.
!     A more elegant fix might be to use more than one Davidson 
!     vectors in the space.
            fac=.50d0*U(2,1); Norm=One/Pythag(One,fac)
            X=Norm*(X+fac*X1)
            XH=Norm*(XH+fac*X1H)
         End If

         Norm=Sqrt(MDot(X,X)); X=X/Norm; XH=XH/Norm

!     Avoid accumulating roundoff error
         If((iter/4)*4==iter)
     $        Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)

!     Print out Convergence Information
         Energy=MDot(X,XH)/MDot(X,X); DE=DE-Energy
!         If(iter==1)Write(6,*)'Step    Energy             Predicted     
!     $     Delta_E          Residual'
!         If((iter/10)*10==iter)
!     $        Write(6,114)iter,Energy,Eig(1),DE,fac

         XH=0d0
         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)

 114     Format('##',I4,5f18.13)
 113     Format(I4,3f18.13)
      End Do

!     Calculate Uncertainty in Energy
      Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)
      Uncertainty(iS)=MDot(XH,XH)/MDot(XH,X)/MDot(XH,X)-1
!      Write(6,*)'-----Fractional Uncertainty in Energy:'
!     $     ,Uncertainty(iS)

!     Store this Eigenvecor and Eigenvalue
      Xi(:,:,iS)=X; Ei(iS)=Energy

!     End Big Loop Over States
      End Do


!      Write(6,*) '-------------Summary of Results-------------'
!      Write(6,*) '       State        Energy         ',
!     $     '       Spin Symmetry         Uncertainty'
!      Do i=1,NS
!         Write(6,*)i,Ei(i),Sym(i),Uncertainty(i)
!      End Do
!      Write(6,*) '-------------Summary of Results-------------'

! Get One Particle (Transition) Density Matrix 
      T2=0d0; P=0d0
      Do iS=1,NS
         X=Xi(:,:,iS)
         Do i=1,N; Do j=i,N
            T1=Zero
            T1(i,j)=T1(i,j)+0.25d0
            T1(j,i)=T1(j,i)+0.25d0
            Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,T2,XH)
            Do jS=1,NS
               P(i,j,iS,jS)=MDot(Xi(:,:,jS),XH)
               P(j,i,iS,jS)=P(i,j,iS,jS)
            End Do
         End Do; End Do
      End Do

      X = Xi(:,:,1)


      T1=0d0
!The following code needs T1 to remain 0
! Generate complete 2-particle density matrix
      P2 = 0d0 
      X=Xi(:,:,1)
      Do i=1,N; Do j=1,N; Do k=1,N; Do l=1,N
        T2=Zero
        T2(i,j,k,l)=1d0
        Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,T2,XH)
        P2(i,j,k,l) = MDot(Xi(:,:,1),XH)
      End Do; End Do; End Do; End Do
 
!     Project onto Site #1  By Brute Force
!     @@TAV Could do the same for two sites relatively easily from Ex2
      X1a=0d0; X2=0d0; X3=0d0; X4=0d0
      Do ii=1,Max1
         i=Ex1(ii,1,1)
         Do jj=1,Max1
            j=Ex1(jj,1,1)
            X1a(i,j)=X(i,j)
         End Do
         Do kk=1,Nstr
            X2(i,kk)=X(i,kk)
            X3(kk,i)=X(kk,i)
         End Do
      End Do
      X2=X2-X1a; X3=X3-X1a
      X4=X-X1a-X2-X3
      
!FCI Overlap    
      S1s(1,1) = MDot(X1a,X1a); S1s(1,2) = MDot(X1a,X2)
      S1s(1,3) = MDot(X1a,X3); S1s(1,4) = MDot(X1a,X4)
      S1s(2,1) = MDot(X2,X1a); S1s(2,2) = MDot(X2,X2)
      S1s(2,3) = MDot(X2,X3); S1s(2,4) = MDot(X2,X4)
      S1s(3,1) = MDot(X3,X1a); S1s(3,2) = MDot(X3,X2)
      S1s(3,3) = MDot(X3,X3); S1s(3,4) = MDot(X3,X4)
      S1s(4,1) = MDot(X4,X1a); S1s(4,2) = MDot(X4,X2)
      S1s(4,3) = MDot(X4,X3); S1s(4,4) = MDot(X4,X4)

!     Normalize States
      X1a=X1a/Sqrt(MDot(X1a,X1a)); X2=X2/Sqrt(MDot(X2,X2));
      X3=X3/Sqrt(MDot(X3,X3)); X4=X4/Sqrt(MDot(X4,X4));

!     Compute Matrix Elements      
      Call UHX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X1a,hnat,hnat, 
     $  Vnat,Vnat,Vnat,X1Ha)
      Call UHX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X2,hnat,hnat, 
     $  Vnat,Vnat,Vnat,X2H)
      Call UHX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X3,hnat,hnat, 
     $  Vnat,Vnat,Vnat,X3H)
      Call UHX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X4,hnat,hnat, 
     $  Vnat,Vnat,Vnat,X4H)

      H1s(1,1) = MDot(X1a,X1Ha); H1s(1,2) = MDot(X1a,X2H)
      H1s(1,3) = MDot(X1a,X3H); H1s(1,4) = MDot(X1a,X4H)
      H1s(2,1) = MDot(X2,X1Ha); H1s(2,2) = MDot(X2,X2H)
      H1s(2,3) = MDot(X2,X3H); H1s(2,4) = MDot(X2,X4H)
      H1s(3,1) = MDot(X3,X1Ha); H1s(3,2) = MDot(X3,X2H)
      H1s(3,3) = MDot(X3,X3H); H1s(3,4) = MDot(X3,X4H)
      H1s(4,1) = MDot(X4,X1Ha); H1s(4,2) = MDot(X4,X2H)
      H1s(4,3) = MDot(X4,X3H); H1s(4,4) = MDot(X4,X4H)


      
!      V = 0d0 !I guess we are calculating the HF energy here? 
!      Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)
      !Efci = MDot(Xi(:,:,1),XH)

      Return
      End Subroutine

      Subroutine Pdm2(N,No,N2,Max1,Max2,Nstr,X,T1,T2,XH,P,P2)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Max1,Max2,Nstr
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,l
      Real*8 X(Nstr,Nstr),XH(Nstr,Nstr),T1(N,N),T2(N,N,N,N)
      Real*8 P2(N,N,N,N),P(N,N)

      Call IString(N,No,N2,Max1,Max2,Ex1,Ex2)

      T2=0d0; P=0d0
      Do i=1,N; Do j=i,N
         T1=0d0
         T1(i,j)=T1(i,j)+0.25d0
         T1(j,i)=T1(j,i)+0.25d0
         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,T2,XH)
         P(i,j)=MDot(X,XH)
         P(j,i)=P(i,j)
      End Do; End Do

      T1=0d0
!The following code needs T1 to remain 0
! Generate complete 2-particle density matrix
      P2 = 0d0 
      Do i=1,N; Do j=1,N; Do k=1,N; Do l=1,N
        T2=0d0
        T2(i,j,k,l)=1d0
        Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,T2,XH)
        P2(i,j,k,l) = MDot(X,XH)
      End Do; End Do; End Do; End Do
      Return
      End Subroutine

      Subroutine Hx(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,Y)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Max1,Max2,Nstr
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl
      Integer*4 AEx1(Max1,N,N),IfZero(N,N),IfSym
      Real*8 X(Nstr,Nstr),h(N,N),V(N,N,N,N),Y(Nstr,Nstr)
      Real*8 Vtmp,VS,VSS,htmp,hS,Tmp,Spin
      Real*8 Xtmp(Max1,Nstr),Ytmp(Max1,NStr)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      Y=Zero; AEx1=Abs(Ex1)

!     Check Spin Symmetry of X
      Spin=One; Y=X-Transpose(X); Tmp=MDot(Y,Y)
      If(Tmp>0.1d0) Spin=-One
      Y=Zero;

! Check for Zero Blocks in V
      IfZero=0
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If(Abs(V(i,j,k,l))>1.0d-10)Go To 20
         End Do; End Do
         IfZero(i,k)=1
 20      Continue
      End Do; End Do

! Check Symmetry of V
      IfSym=1
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If (Abs(V(i,j,k,l)-V(j,i,l,k))>1.0d-10) Then
               IfSym=0
               Go To 40
            End If
         End Do; End Do
      End Do; End Do
 40   Continue

! One Electron Part
      Do i=1,N
         Do j=1,N
            htmp=h(i,j)
            If(Abs(htmp)<1.0d-10)Cycle
            If(i==j) Then
               iiMax=NcR(N-1,No-1)
            Else
               iiMax=NcR(N-2,No-1)
            End If
            Do ii=1,iiMax
               I1=AEx1(ii,i,j)
               I2=AEx1(ii,j,i)
               hS=htmp
               If(I1/=Ex1(ii,i,j))hS=-hS
               Call DaxPy(Nstr,hS,X(1,I1),1,Y(1,I2),1)
            End Do
         End Do
      End Do
 143  Format(4I4,3f14.8)

! Same Spin Two-Electron Part
      ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; kl=0
         Do k=1,N; Do l=k+1,N; kl=kl+1
            Vtmp=(V(i,j,k,l)-V(i,j,l,k)-V(j,i,k,l)+V(j,i,l,k))/Two
            If(abs(Vtmp)<1.0d-10)Go To 100 
            If(i==k.and.j==l) Then
               iiMax=NcR(N-2,No-2)
            Else If(i==k) Then  !i=2,3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(i==l) Then  !i=1-3
               iiMax=NcR(N-3,No-2)
            Else If(j==k) Then  !i=3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(j==l) Then  !i=1,2,3,4,7
               iiMax=NcR(N-3,No-2)
            Else
               iiMax=NcR(N-4,No-2)
            End If
            Do ii=1,iiMax
               I1=Abs(Ex2(ii,ij,kl))
               I2=Abs(Ex2(ii,kl,ij))
               VS=Vtmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
               Call DaxPy(Nstr,VS,X(1,I1),1,Y(1,I2),1)
            End Do
 100        Continue
         End Do; End Do
      End Do; End Do

! Opposite Spin Two-Electron Part
      ik=0
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together elements of Xtmp
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            VS=One; 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            Do jj=1,Nstr
               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
            
! Collect Elements of Ytmp
         Ytmp=Zero
         Do j=1,N; Do l=1,N; jl=jl+1
            If(IfSym==1) Then
               If(ik<jl) Cycle
               If(ik==jl)Vtmp=V(i,j,k,l)/Two
               If(ik>jl)Vtmp=V(i,j,k,l)
            Else
               Vtmp=V(i,j,k,l)/Two
            End If

            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If

            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
               Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
            End Do
         End Do; End Do

! Scatter Elements of Y
         Do ii=1,iiMax
            I1=AEx1(ii,k,i)
            Do jj=1,Nstr
               Y(I1,jj)=Y(I1,jj)+Ytmp(ii,jj)
            End Do
         End Do
      End Do; End Do
! Enforce MS=0
      Y=Y+Spin*Transpose(Y)
      Return
      End Subroutine

      Subroutine UHx(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,
     $     ha,hb,Vaa,Vab,Vbb,Y)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Max1,Max2,Nstr
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl
      Integer*4 AEx1(Max1,N,N),IfZero(N,N),IfSym
      Real*8 X(Nstr,Nstr),Y(Nstr,Nstr)
      Real*8 ha(N,N),hb(N,N),VAA(N,N,N,N),VAB(N,N,N,N),VBB(N,N,N,N)
      Real*8 Vtmp,Vatmp,Vbtmp,VS,VSS,hatmp,hbtmp,hS,Tmp,Spin
      Real*8 Xtmp(Max1,Nstr),Ytmp(Max1,NStr)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      Y=Zero; AEx1=Abs(Ex1);
     
! Check for Zero Blocks in V
      IfZero=0
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If(Abs(VAA(i,j,k,l))>1.0d-10)Go To 20
            If(Abs(VBB(i,j,k,l))>1.0d-10)Go To 20
            If(Abs(VAB(i,j,k,l))>1.0d-10)Go To 20
         End Do; End Do
         IfZero(i,k)=1
 20      Continue
      End Do; End Do

! Check Symmetry of V
      IfSym=1
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If (Abs(Vab(i,j,k,l)-Vab(j,i,l,k))>1.0d-10) Then
               IfSym=0
               Go To 40
            End If
         End Do; End Do
      End Do; End Do
 40   Continue

! One Electron Part
      Do i=1,N
         Do j=1,N
            hatmp=ha(i,j); hbtmp=hb(i,j)
            If(Abs(hatmp)+Abs(hbtmp)<1.0d-10)Cycle
            If(i==j) Then
               iiMax=NcR(N-1,No-1)
            Else
               iiMax=NcR(N-2,No-1)
            End If
            Do ii=1,iiMax
               I1=AEx1(ii,i,j)
               I2=AEx1(ii,j,i)
               hS=hbtmp
               If(I1/=Ex1(ii,i,j))hS=-hS
               Call DaxPy(Nstr,hS,X(1,I1),1,Y(1,I2),1)
               hS=hatmp
               If(I1/=Ex1(ii,i,j))hS=-hS
               Call DaxPy(Nstr,hS,X(I1,1),Nstr,Y(I2,1),Nstr)
            End Do
         End Do
      End Do
 143  Format(4I4,3f14.8)

! Same Spin Two-Electron Part
      ij=0; 
      Do i=1,N; Do j=i+1,N; ij=ij+1; kl=0
         Do k=1,N; Do l=k+1,N; kl=kl+1
            Vatmp=(Vaa(i,j,k,l)-Vaa(i,j,l,k)
     $           -Vaa(j,i,k,l)+Vaa(j,i,l,k))/Two
            Vbtmp=(Vbb(i,j,k,l)-Vbb(i,j,l,k)
     $           -Vbb(j,i,k,l)+Vbb(j,i,l,k))/Two
            If(Abs(Vatmp)+Abs(Vbtmp)<1.0d-10)Go To 100 
            If(i==k.and.j==l) Then
               iiMax=NcR(N-2,No-2)
            Else If(i==k) Then  !i=2,3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(i==l) Then  !i=1-3
               iiMax=NcR(N-3,No-2)
            Else If(j==k) Then  !i=3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(j==l) Then  !i=1,2,3,4,7
               iiMax=NcR(N-3,No-2)
            Else
               iiMax=NcR(N-4,No-2)
            End If
            Do ii=1,iiMax
               I1=Abs(Ex2(ii,ij,kl))
               I2=Abs(Ex2(ii,kl,ij))
               VS=Vatmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
               Call DaxPy(Nstr,VS,X(I1,1),Nstr,Y(I2,1),Nstr)
               VS=Vbtmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
               Call DaxPy(Nstr,VS,X(1,I1),1,Y(1,I2),1)
            End Do
 100        Continue
         End Do; End Do
      End Do; End Do

! Opposite Spin Two-Electron Part
      ik=0; X=Transpose(X)
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together elements of Xtmp
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            VS=One; 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            Do jj=1,Nstr
               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
            
! Collect Elements of Ytmp
         Ytmp=Zero
         Do j=1,N; Do l=1,N; jl=jl+1
            If(IfSym==1) Then
               If(ik<jl) Cycle
               If(ik==jl)Vtmp=Vab(i,j,k,l)/Two
               If(ik>jl)Vtmp=Vab(i,j,k,l)
            Else
               Vtmp=Vab(i,j,k,l)/Two
            End If

            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If

            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
               Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
            End Do
         End Do; End Do

! Scatter Elements of Y
         Do ii=1,iiMax
            I1=AEx1(ii,k,i)
            Do jj=1,Nstr
               Y(jj,I1)=Y(jj,I1)+Ytmp(ii,jj)
            End Do
         End Do
      End Do; End Do

! Opposite Spin Two-Electron Part
      ik=0; X=Transpose(X)
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together elements of Xtmp
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            VS=One; 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            Do jj=1,Nstr
               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
            
! Collect Elements of Ytmp
         Ytmp=Zero
         Do j=1,N; Do l=1,N; jl=jl+1
            If(IfSym==1) Then
               If(ik<jl) Cycle
               If(ik==jl)Vtmp=Vab(i,j,k,l)/Two
               If(ik>jl)Vtmp=Vab(i,j,k,l)
            Else
               Vtmp=Vab(i,j,k,l)/Two
            End If

            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If

            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
               Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
            End Do
         End Do; End Do

! Scatter Elements of Y
         Do ii=1,iiMax
            I1=AEx1(ii,k,i)
            Do jj=1,Nstr
               Y(I1,jj)=Y(I1,jj)+Ytmp(ii,jj)
            End Do
         End Do
      End Do; End Do

      Return
      End Subroutine

      Subroutine GetH0(N,No,N0,N2,Max1,Max2,Nstr,Ex1,Ex2,Istr,h,V,H0)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N0,N2,Max1,Max2,Nstr,m,If0,IY(N0),mmax
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2),Istr(N0),RIstr(Nstr)
      Integer*4 i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl
      Integer*4 AEx1(Max1,N,N),AEx2(Max2,N2,N2),IfZero(N,N)
      Real*8 h(N,N),V(N,N,N,N),H0(N0,N0,N0,N0),H0tmp(N0,N0,N0,N0)
      Real*8 Vtmp,VS,VSS,htmp,hS
      Real*8 Stmp(N0)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      AEx1=Abs(Ex1); AEx2=Abs(Ex2); H0=Zero
! Reverse String Ordering
      RIstr=0
      Do m=1,N0
         RIstr(Istr(m))=m
      End Do

! Remove terms that are not in H0
      Do i=1,N; Do j=1,N; Do ii=1,Max1
         If0=0
         Do m=1,N0
            If(Istr(m)==AEx1(ii,i,j))If0=1
         End Do
         If(If0==0) AEx1(ii,i,j)=0
      End Do; End Do; End Do

      Do i=1,N2; Do j=1,N2; Do ii=1,Max2
         If0=0
         Do m=1,N0
            If(Istr(m)==AEx2(ii,i,j))If0=1
         End Do
         If(If0==0) AEx2(ii,i,j)=0
      End Do; End Do; End Do

! Check for Zero Blocks in V
      IfZero=0
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If(Abs(V(i,j,k,l))>1.0d-10)Go To 20
         End Do; End Do
         IfZero(i,k)=1
 20      Continue
      End Do; End Do

! One Electron Part
      Do i=1,N
         Do j=1,N
            htmp=h(i,j)
            If(Abs(htmp)<1.0d-10)Cycle
            If(i==j) Then
               iiMax=NcR(N-1,No-1)
            Else
               iiMax=NcR(N-2,No-1)
            End If
            Do ii=1,iiMax
               I1=AEx1(ii,i,j)
               I2=AEx1(ii,j,i)
               If(I1==0.or.I2==0) Cycle
               hS=htmp
               If(I1/=Ex1(ii,i,j))hS=-hS
!               Call DaxPy(Nstr,hS,X(1,I1),1,Y(1,I2),1)
               Do m=1,N0
                  H0(m,RIstr(I2),m,RIstr(I1))=
     $                 H0(m,RIstr(I2),m,RIstr(I1))+hS
               End Do
            End Do
         End Do
      End Do
      
! Same Spin Two-Electron Part
      ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; kl=0
         Do k=1,N; Do l=k+1,N; kl=kl+1
            Vtmp=(V(i,j,k,l)-V(i,j,l,k)-V(j,i,k,l)+V(j,i,l,k))/Two
            If(abs(Vtmp)<1.0d-10)Go To 100 
            If(i==k.and.j==l) Then
               iiMax=NcR(N-2,No-2)
            Else If(i==k) Then  !i=2,3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(i==l) Then  !i=1-3
               iiMax=NcR(N-3,No-2)
            Else If(j==k) Then  !i=3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(j==l) Then  !i=1,2,3,4,7
               iiMax=NcR(N-3,No-2)
            Else
               iiMax=NcR(N-4,No-2)
            End If
            Do ii=1,iiMax
               I1=AEx2(ii,ij,kl)
               I2=AEx2(ii,kl,ij)
               If(I1==0.or.I2==0) Cycle
               VS=Vtmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
!               Call DaxPy(Nstr,VS,X(1,I1),1,Y(1,I2),1)
               Do m=1,N0
                  H0(m,RIstr(I2),m,RIstr(I1))=
     $                 H0(m,RIstr(I2),m,RIstr(I1))+VS
               End Do
            End Do
 100        Continue
         End Do; End Do
      End Do; End Do
      
! Opposite Spin Two-Electron Part
      ik=0
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together phases in Stmp
         Stmp=Zero; IY=0
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            I2=AEx1(ii,k,i)
            If(I1==0.or.I2==0)Cycle
            VS=One 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            IY(RIstr(I1))=RIstr(I2)
            Stmp(RIstr(I1))=VS
            Do jj=1,Nstr
!               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
         mmax=0
         Do m=1,N0
            If(IY(m)/=0)mmax=m
            If(IY(m)==0)IY(m)=1
         End Do
         
! Collect Elements of H0
         Do j=1,N; Do l=1,N; jl=jl+1
            If(ik<jl) Cycle
            If(ik==jl)Vtmp=V(i,j,k,l)/Two
            If(ik>jl)Vtmp=V(i,j,k,l)
            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If
            
            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               If(J1==0.or.J2==0)Cycle
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
!                  Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
               Do m=1,mmax
                  H0(IY(m),RIstr(J2),m,RIstr(J1))=
     $                 H0(IY(m),RIstr(J2),m,RIstr(J1))+VS*Stmp(m)
               End Do
            End Do
         End Do; End Do
      End Do; End Do
      Do i=1,N0; Do j=1,N0; Do k=1,N0; Do l=1,N0
         H0tmp(i,j,k,l)=H0(i,j,k,l)+H0(j,i,l,k)
      End Do; End Do; End Do; End Do
      H0=H0tmp
      Return
      End Subroutine

      Subroutine GetIstr(N,No,Nstr,N0,Hd,Istr)
      Integer*4 N,No,Nstr,N0,Istr(N0),itmp,a,b
      Real*8 Hd(Nstr,Nstr),Min

      Do a=1,N0
         Min=500.0d0; Itmp=1
         Do i=1,Nstr
            If(Hd(i,i)< Min) Then
               Min=Hd(i,i)
               Itmp=i
            End If            
         End Do
      
         Istr(a)=Itmp

         Hd(Istr(a),Istr(a))=Hd(Istr(a),Istr(a))+1000.0d0

      End Do

      Do a=1,N0
         Hd(Istr(a),Istr(a))=Hd(Istr(a),Istr(a))-1000.0d0
      End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     GetHd(N,No,Nstr,h,V,Iocca,Ioccb,Hd,IRecur)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Nstr,Iocca(No),Ioccb(No),IRecur
      Integer*4 i,j,imin,jmin,k,ka,kb,l,la,lb,Isigna,Isignb
      Real*8 h(N,N),V(N,N,N,N),Hd(Nstr,Nstr),tmp

      If(IRecur==1) Then
         imin=1
         jmin=1
      Else
         imin=Iocca(IRecur-1)+1
         jmin=Ioccb(IRecur-1)+1
      End If

      Do i=imin,N
         Iocca(IRecur)=i
         Do j=jmin,N
            Ioccb(IRecur)=j

            If(IRecur==No) Then
               tmp=0.0d0
               Do k=1,No
                  ka=Iocca(k); kb=Ioccb(k)
! Spin Contaminated Elements
                  tmp=tmp+h(ka,ka)+h(kb,kb)
                  Do l=1,No
                     la=Iocca(l); lb=Ioccb(l)
                     tmp=tmp+.50d0*(
     $                    +V(ka,la,ka,la)-V(ka,la,la,ka)
     $                    +V(ka,lb,ka,lb)+V(kb,la,kb,la)
     $                    +V(kb,lb,kb,lb)-V(kb,lb,lb,kb))
                  End Do
               End Do
!               Write(6,*)'@@tmp',tmp,Iocca,Ioccb
               Hd(Index(No,Iocca),Index(No,Ioccb))=tmp
               
            Else
               Call GetHd(N,No,Nstr,h,V,Iocca,Ioccb,Hd,IRecur+1)
            End If
         End Do
      End Do

      Return
      End Subroutine


      Subroutine IString(N,No,N2,Max1,Max2,Ex1,Ex2)
      Implicit None
      Integer*4 N,No,N2,Max1,Max2
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,a,b,ij,ab,IRecur,m,Ind,Iocc(No)

! Find Strings that differ by a Single Excitation
      Ex1=0
      Do i=1,N; Do a=1,N
         m=1; Ind=1
         Call RecurEx1(N,No,i,a,Max1,Iocc,Ex1,0,m,Ind)
      End Do; End Do

! Find Strings that differ by a Double Excitation
      Ex2=0; ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; ab=0; Do a=1,N; Do b=a+1,N
         ab=ab+1; m=1; Ind=1
         Call RecurEx2(N,No,N2,i,j,a,b,ij,ab,
     $        Max2,Iocc,Ex2,0,m,Ind)
      End Do; End Do; End Do; End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     RecurEx2(N,No,N2,i,j,a,b,ij,ab,Max2,Iocc,Ex2,IRecur,m,Ind)
      Implicit None
      Integer*4 N,No,N2,i,j,a,b,ij,ab,Max2,Iocc(No),Ex2(Max2,N2,N2)
      Integer*4 IRecur,m,Ind,ii,jj,IsiO,IsjO,IsaO,IsbO
      Integer*4 Isubst(No),iind,jind,aind,bind,Isign,iimin

      If(IRecur==0)iimin=1
      If(IRecur/=0)iimin=Iocc(IRecur)+1

      Do ii=iimin,N
         Iocc(IRecur+1)=ii
         
         If(IRecur==No-1) Then
! End Recursion Condition
            IsiO=0; IsjO=0; IsaO=0; IsbO=0

! Is abj*i* |Iocc> Non-Zero?
            Do jj=1,No
               If(Iocc(jj)==i) Then
                  IsiO=1
               Else If(Iocc(jj)==j) Then
                  IsjO=1
               Else If(Iocc(jj)==a) Then
                  IsaO=1
               Else If(Iocc(jj)==b) Then
                  IsbO=1
               End If
            End Do
            If(IsiO==1 .and. IsjO==1 .and. IsaO==0 .and. IsbO==0) Then
               Isubst=Iocc
               Do jj=1,No
                  If(Isubst(jj)==i) Then
                     Isubst(jj)=a
                  Else If(Isubst(jj)==j) Then
                     Isubst(jj)=b
                  End If
               End Do

               Call ISort(No,Isubst,Isign)

               Ex2(m,ab,ij)=Isign*Ind
               m=m+1
            End If
            Ind=Ind+1
         Else
            Call RecurEx2(N,No,N2,i,j,a,b,ij,ab,
     $           Max2,Iocc,Ex2,IRecur+1,m,Ind)
         End If

      End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     RecurEx1(N,No,i,a,Max1,Iocc,Ex1,IRecur,m,Ind)
      Implicit None
      Integer*4 N,No,i,a,Max1,Iocc(No),Ex1(Max1,N,N),IRecur,m,Ind,ii
      Integer*4 jj,IsiO,IsaO,Isubst(No),aind,iind,Isign,iimin

      If(IRecur==0)iimin=1
      If(IRecur/=0)iimin=Iocc(IRecur)+1

      Do ii=iimin,N
         Iocc(IRecur+1)=ii
         
         If(IRecur==No-1) Then
! End Recursion Condition
            IsiO=0; IsaO=0
            Do jj=1,No
               If(Iocc(jj)==i) Then 
                  IsiO=1
               Else If(Iocc(jj)==a) Then
                  IsaO=1
               End If
            End Do

            If(IsiO==1 .and. IsaO==0) Then
               Isubst=Iocc
               Do jj=1,No
                  If(i==Isubst(jj)) Isubst(jj)=a
               End Do

               Call Isort(No,Isubst,Isign)

               Ex1(m,a,i)=Isign*Ind
               m=m+1
            End If
            Ind=Ind+1
         Else
            Call RecurEx1(N,No,i,a,Max1,Iocc,Ex1,IRecur+1,m,Ind)
         End If

      End Do

      Return
      End Subroutine

      Subroutine GSOrtho(N,C)
! Takes a set of Non-Orthogonal, Non-Normalized vectors and makes
!   them orthonormal by the Gramm-Schmidt procedure
!   Vectors are C(:,p)
      Implicit None
      Integer*4 N,p,q
      Real*8 C(N,N),pNorm,DDot

      Do p=1,N
         Do q=1,p-1
            C(:,p)=C(:,p)-DDot(N,C(:,p),1,C(:,q),1)*C(:,q)
         End Do
         pNorm=0.0d0
         Do q=1,N
            pNorm=pNorm+C(q,p)*C(q,p)
         End Do

         C(:,p)=C(:,p)/Sqrt(pNorm)
      End Do

      Return
      End Subroutine

      Subroutine ISort(N,X,Isign)
!  Sorts an Integer arrary X by straight insertion
!     Isign is the sign of the Permutation needed to bring it in order
      Implicit None
      Integer*4 N,X(N),i,j,tmp,Isign

      Isign=1
      Do i=2,N
         tmp=X(i)
         Do j=i-1,1,-1
            If(X(j).le.tmp)Go To 100
! Transposing indicies ...
            X(j+1)=X(j)
! Gives us a minus sign.
            Isign=-Isign
         End Do
         j=0
 100     X(j+1)=tmp
      End Do

      Return
      End Subroutine
      Subroutine Sort(N,X,Isign)
!  Sorts a Real arrary X by straight insertion
!     Isign is the sign of the Permutation needed to bring it in order
      Implicit None
      Integer*4 N,i,j,Isign
      Real*8 X(N),tmp

      Isign=1
      Do i=2,N
         tmp=X(i)
         Do j=i-1,1,-1
            If(X(j).le.tmp)Go To 100
! Transposing indicies ...
            X(j+1)=X(j)
! Gives us a minus sign.
            Isign=-Isign
         End Do
         j=0
 100     X(j+1)=tmp
      End Do

      Return
      End Subroutine

      Subroutine GetEig2(H,E,U)
      Implicit None
! Diagonalize a 2x2 in High Precision
      Real*8 H(2,2),E(2),U(2,2)
      Real*16 HQ(2,2),EQ(2),UQ(2,2)
      Real*16 Two,fac

      Two=2.0q0
      HQ=H
      
      fac=(HQ(1,1)-HQ(2,2))*(HQ(1,1)-HQ(2,2))/4.0q0 + HQ(1,2)*HQ(1,2)
      fac=sqrt(fac)
      EQ(1)=(HQ(1,1)+HQ(2,2))/Two-fac
      EQ(2)=(HQ(1,1)+HQ(2,2))/Two+fac
      
      UQ=0q0; UQ(1,1)=1q0; UQ(2,2)=1q0
      If(Abs(EQ(1)-HQ(1,1)).gt.1q-18) Then
         fac= 1.0q0 + HQ(1,2)/(HQ(1,1)-EQ(1))*HQ(1,2)/(HQ(1,1)-EQ(1))
         fac=sqrt(fac)
         
         UQ(1,1)= -HQ(1,2)/fac/(HQ(1,1)-EQ(1))
         UQ(2,1)= 1.0q0/fac
         
         UQ(1,2)=-UQ(2,1)
         UQ(2,2)=UQ(1,1)
      End If

      E=EQ; U=UQ

      Return 
      End

      Integer*4 Function MyMod(N,M)
      Integer*4 M,N,I

      I=N
      Do While (.not.(I.le.M.and.I.gt.0))
         If(N.gt.M)I=I-M
         If(N.le.0)I=I+M
      End Do

      MyMod=I

      End Function
