      Subroutine Diagonalize(N,A,Eigs,U,ierr)
      Implicit None
      Integer N,ierr,LWork
      Real*8 A(N,N),Eigs(N),U(N,N),Work(15*N+10)
! Get Eigenvalues and eigenvectors of a Real Symmetric Matrix

      LWork=15*N+10; U=A
      Call DSyEV('V','U',N,U,N,Eigs,Work,LWork,ierr)

      end subroutine

      Subroutine DgeMatMul(L,M,N,Y,Z,X,alpha,beta)
      Implicit None
!  X(L,N)= alpha*Y(L,M)*Z(M,N)+beta*(X(L,N)

      Integer*4 L,M,N
      Real*8 Y(L,M),Z(M,N),X(L,N),alpha,beta

!  ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

      Call Dgemm('N','N',L,N,M,alpha,Y(1,1),L,Z(1,1),M,beta,X(1,1),L)

      Return
      End Subroutine


      Subroutine Axb(N,LDA,A,b,M,ierr)
      Implicit None
! Solves Ax=b in the most stable way possible
      Integer*4 N,LDA,M,Ipiv(N),IWork(N),ierr
      Real*8 A(LDA,N),b(N,M)
      Real*8 AF(N,N),R(N),C(N),X(N,M),Rcond,Ferr,Berr,Work(4*N)
      Character*1 EQ
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

! Clear out all the temporary arrays
      Ipiv=0; IWork=0
      AF=zero; R=zero; C=zero; X=zero; Rcond=zero; Ferr=zero 
      Berr=zero; Work=zero


!  Solve  Ax=b
!    In the first spot:
!     'N' - no equilibration
!     'E' - equilibrate A

      Call DGeSVX('N','N',N,M,A,LDA,AF,N,Ipiv,EQ,R,
     $               C,B,N,X,N,Rcond,Ferr,Berr,Work,IWork,ierr)

! Put the result back in b
      b=x
      
      Return
      End

      Real*8 Function factorial(n)
      Implicit None
! factorial - computes n!

      Integer*4 i,n

      factorial=1.0d0
      Do i=1,n
        factorial=factorial*i
      End Do

      Return
      End Function






