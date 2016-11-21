module EFLib

  USE mod_lec_fic
  USE SysLinCSR
  USE Precond
  USE BLASS

  IMPLICIT NONE

  CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Compute_LocalStiffness_Matrix_And_LocalRHS(mesh,K,ALOC,RHSLOC)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  IMPLICIT NONE

  type(msh),intent(in) :: mesh
  INTEGER, intent(in)  :: K

  REAL (kind=8), dimension(3,3), intent(out) :: ALOC
  REAL (kind=8), dimension(3), intent(out) :: RHSLOC

  REAL (kind=8), dimension(2,3) :: GRADPHIREF,GRADPHI
  INTEGER, dimension(3) :: TRI
  REAL (kind=8), dimension(2) :: P1,P2,P3
  REAL (kind=8) :: x1,x2,x3,y1,y2,y3
  INTEGER :: k1,k2,k3
  REAL (kind=8), dimension(2,2) :: MK,MKT,INVMKT
  REAL (kind=8) :: DETMK 

! Reference Stiffness Matrix
!ccccccccccccccccccccccccccc
  GRADPHIREF=reshape((/-1,-1,1,0,0,1/),(/2,3/))

! Nodes coordinates
!cccccccccccccccccc 

   TRI=mesh%triangles(k,1:3)
   k1=TRI(1); k2=TRI(2); k3=TRI(3);
   P1=mesh%POS(k1,1:2);     P2=mesh%POS(k2,1:2);     P3=mesh%POS(k3,1:2)

   x1=P1(1); y1=P1(2)
   x2=P2(1); y2=P2(2)
   x3=P3(1); y3=P3(2)

! Computation of MK, MKT, INVMKT
!cccccccccccccccccccccccccccccccc

   MK(1,1)=x2-x1 ; MK(1,2)=x3-x1 
   MK(2,1)=y2-y1 ; MK(2,2)=y3-y1 

   DETMK=MK(1,1)*MK(2,2)-MK(2,1)*MK(1,2) 

   MKT=TRANSPOSE(MK) 

   INVMKT(1,1)=MKT(2,2)/DETMK;
   INVMKT(1,2)=-MKT(1,2)/DETMK;
   INVMKT(2,1)=-MKT(2,1)/DETMK;
   INVMKT(2,2)=MKT(1,1)/DETMK;
  
! ALOC Computation
!ccccccccccccccccc

   GRADPHI=MATMUL(INVMKT,GRADPHIREF)
   ALOC=ABS(DETMK)/2*MATMUL(TRANSPOSE(GRADPHI),GRADPHI)
  
! RHSLOC Computation
!cccccccccccccccccccc

   RHSLOC(1)=source(x1,y1)
   RHSLOC(2)=source(x2,y2)
   RHSLOC(3)=source(x3,y3)
   RHSLOC=RHSLOC*ABS(DETMK)/6
   
  end subroutine Compute_LocalStiffness_Matrix_And_LocalRHS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Assembling_Global_Stiffness_Matrix_And_RHS(mesh,IND,VAL,NBE,RHS)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  IMPLICIT NONE

  type(msh),intent(in) :: mesh

  INTEGER, dimension(:,:), allocatable, intent(out) :: IND
  REAL (kind=8), dimension(:,:), allocatable, intent(out) :: VAL
  INTEGER, dimension(:), allocatable, intent(out) :: NBE
  REAL (kind=8), dimension(:), allocatable, intent(out) :: RHS

  INTEGER :: NT,NN,K,L,LIG,J,I
  INTEGER, dimension(3) :: TRI
  REAL (kind=8), dimension(3,3) :: ALOC
  REAL (kind=8), dimension(3)   :: RHSLOC
  
  NT=mesh%nbTriangles
  NN=mesh%nbNod
! The value "10" is the maximal number of neighbours that a node can have.
  allocate(IND(NN,10),VAL(NN,10),NBE(NN))
  allocate(RHS(NN))
  IND=0; VAL=0; NBE=0; RHS=0

  DO K=1,NT
    CALL Compute_LocalStiffness_Matrix_And_LocalRHS(mesh,K,ALOC,RHSLOC)
    TRI=mesh%triangles(K,1:3)
    DO L=1,3
     LIG=TRI(L) ! number of the current node
     RHS(LIG)=RHS(LIG)+RHSLOC(L)
     IF (NBE(LIG).eq.0) then
! If the node was never been encoutered
! One put every nodes of the triangle in the connectivity table
       IND(LIG,1:3)=TRI 
! The number of elements in the row corresponding to the current node
! becomes 3
       NBE(LIG)=3
! Assembling
       VAL(LIG,1:3)= ALOC(L,1:3)
     ELSE
! The node was already been visited in a previous triangle
! It is therefore needed to examine every neighbouring nodes 
! (including itself)
       OUT: DO I=1,3
        ! This loop allows to avoid to add already existing nodes
        ! in the connectivity table. If a connectivity already
        ! exists, one adds to the matrix component the local matrix.
         DO J=1,NBE(LIG)
          IF (IND(LIG,J).EQ.TRI(I)) THEN
            VAL(LIG,J)=VAL(LIG,J)+ALOC(L,I) 
            CYCLE OUT
          END IF
         END DO
        ! If no connection was present, we create it.
        ! One assembles the element thanks to the local matrix.
              NBE(LIG)=NBE(LIG)+1
              IND(LIG,NBE(LIG))=TRI(I)
              VAL(LIG,NBE(LIG))=ALOC(L,I)
       END DO OUT
    END IF
   END DO
  END DO

end subroutine Assembling_Global_Stiffness_Matrix_And_RHS

 
!ccccccccccccccccccccccccccccccccccccccccccccxxxccccccc
  subroutine Put_In_CSR_Format(IND,VAL,NBE,NNE,A,JA,IA)
!cccccccccccccccccccccccccccccccccccccccccxxxcccccccccc

  IMPLICIT NONE  

  INTEGER, dimension(:,:), intent(in), allocatable :: IND
  REAL (kind=8), dimension(:,:),intent(in), allocatable :: VAL
  INTEGER, dimension(:),intent(in), allocatable :: NBE

  REAL (kind=8), dimension(:), allocatable, intent(out) :: A
  INTEGER, dimension(:), allocatable, intent(out) :: IA, JA
  INTEGER, intent(out) :: NNE

  INTEGER :: J,LIG
  INTEGER, dimension(:), allocatable :: iwork
  INTEGER :: NN  

  NN=size(IND,1)
  NNE=SUM(NBE)
!  write(6,*) 'Number of non nul elements in the matrix: ',NNE
  ALLOCATE(A(NNE),JA(NNE),IA(NN+1))
  A=0.; JA=0; IA=0
  J=1
  DO LIG=1,NN
     A(J:J+NBE(LIG)-1)=VAL(LIG,1:NBE(LIG))
     JA(J:J+NBE(LIG)-1)=IND(LIG,1:NBE(LIG))
     IF (LIG==1) THEN
        IA(LIG)=1
     ELSE
        IA(LIG)=IA(LIG-1)+NBE(LIG-1)
     END IF
     J=J+NBE(LIG)
  END DO
  IA(NN+1)=IA(NN)+NBE(NN)

  ALLOCATE(iwork(max(nn+1,2*nne)))
  CALL CSORT(nn,a,ja,ia,iwork,.true.) 
  DEALLOCATE(iwork)

end subroutine Put_In_CSR_Format

!ccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Boundary_Conditions(mesh,IA,JA,A,RHS)
!cccccccccccccccccccccccccccccccccccccccccccccccccc
!! Take Dirichlet BC into account
!! Identification of nodes
!! For a node I on the boundary
!! It is needed to:
!!      - put 0 on every component on the row I excepted for column I
!!         which takes value 1
!!      - put 0 on every component on the column I excepted for row I
!!         which takes value 1
!! WHERE(JA.EQ.I) A=0.0
!! A(IA(I):IA(I+1)-1)=0.0
!! GETELM (I,I,a,ja,ia,iadd,.TRUE.)
!! A(IADD)=1.0
!! - One compute the complete stiffness matrix and the nodal forces Bi
!! - One multiplies the k-th column of the stiffness matrix by the
!!  value Tk, and one subtracts it to the nodal forces vector.
!!   B = B - A*(0,...,0,Tk,0,...0)
!!   call amux(n,x,y,a,ja,ia)
!! - The k-th row and k-th column are replaced by 0
!!   WHERE(JA.EQ.I) A=0.0
!!   A(IA(I):IA(I+1)-1)=0.0
!! - The component Akk is replaced by 1.
!!   GETELM (I,I,a,ja,ia,iadd,.TRUE.)
!!   A(IADD)=1.0
!! - The component Bk is replaced by Tk.
!!   B(k)=Tk

     IMPLICIT NONE

     type(msh),intent(in) :: mesh

     INTEGER, dimension(:), allocatable, intent(inout) :: IA,JA
     REAL (kind=8), dimension(:), allocatable, intent(inout) :: A
     REAL (kind=8), dimension(:), allocatable, intent(inout) :: RHS

     INTEGER, dimension(:,:), allocatable :: boundary_nodes
     INTEGER, dimension(:), allocatable :: vu
     REAL (kind=8), dimension(:), allocatable ::V1,V2
     REAL(kind=8) :: cl,val_cl,g1
     INTEGER :: k,j,pp1,pp2,type,Zone_Phys_Bord,lig,NN,iadd

   NN=mesh%nbNod
   allocate(boundary_nodes(mesh%nbLines,2))
   allocate(vu(NN))
   boundary_nodes=0
   vu=0
   k=0
   do j=1,mesh%nbLines
    pp1=mesh%LINES(j,1)
    pp2=mesh%LINES(j,2)
    type=mesh%LINES(j,3)
    if (vu(pp1).eq.0) then
        k=k+1
        boundary_nodes(k,1)=pp1
        boundary_nodes(k,2)=type
        vu(pp1)=k
     end if
     if (vu(pp2).eq.0) then
        k=k+1
        boundary_nodes(k,1)=pp2
        boundary_nodes(k,2)=type
        vu(pp2)=k
     end if
  end do

 allocate(V1(NN),V2(NN))
  V1=0.0 ; V2=0.0
  cl=0.0 ! MODIFIER ICI POUR LES CL DIRICHLET...
  Zone_Phys_Bord=1001
  do j=1,mesh%nbLines
     !! Identification of the boundary node to treat
     lig=boundary_nodes(j,1)
     type=boundary_nodes(j,2)
     if (type.eq.Zone_Phys_Bord) then
       val_CL=CL
     end if
     !! - One multiplies the k-th column of the stiffness matrix by the
     !!  value Tk, and one subtracts it to the nodal forces vector.
     !!   B = B - A*(0,...,0,Tk,0,...0)
       V1(lig)=Val_CL
       call amux(nn,v1,v2,a,ja,ia)
       !print *,maxval(v2)
        rhs=rhs-v2
        v1(lig)=0.0
     !! - The k-th row and k-th column are replaced by 0
     !!   WHERE(JA.EQ.I) A=0.0
        WHERE(JA.EQ.LIG) A=0.0
        A(IA(LIG):IA(LIG+1)-1)=0.0
     !! - The component Akk is replaced by 1.
       g1=GETELM (lig,lig,a,ja,ia,iadd,.true.)
       A(IADD)=1.0
     !! - The component Bk is replaced by Tk.
       RHS(LIG)=Val_CL
     end do

  deallocate(boundary_nodes,vu,v1,v2)

end subroutine Boundary_Conditions

!ccccccccccccccccccccccccccccc
function source(x,y) result(z)
!ccccccccccccccccccccccccccccc
      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: z
         z=-2*(y-3)*(y+3)-2*(x-3)*(x+3)
      return
end function source
 
!ccccccccccccccccccccccccccccc
function solexacte(x,y) result(z)
!ccccccccccccccccccccccccccccc
      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: z
         z=(x-3)*(x+3)*(y-3)*(y+3)
      return
end function solexacte

!ccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Compute_Error1(mesh,U,err)
!cccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      type(msh),intent(in) :: mesh
      real (kind=8), dimension(:), allocatable, intent(in) :: U
      real (kind=8), intent(out) :: err     

      real (kind=8) :: error, valnorm, solexact, solnum, area, XG, YG
      integer :: K
      INTEGER, dimension(3) :: TRI
      REAL (kind=8), dimension(2,3) :: GRADPHIREF,GRADPHI
      REAL (kind=8), dimension(2) :: P1,P2,P3
      REAL (kind=8) :: x1,x2,x3,y1,y2,y3
      INTEGER :: k1,k2,k3
      REAL (kind=8), dimension(2,2) :: MK

      error=0.0
      valnorm=0.0
      DO K=1,mesh%nbTriangles
! Computation of the area of the triangle
        TRI=mesh%triangles(K,1:3)
        k1=TRI(1); k2=TRI(2); k3=TRI(3);
        P1=mesh%POS(k1,1:2);     P2=mesh%POS(k2,1:2);     P3=mesh%POS(k3,1:2)
        x1=P1(1); y1=P1(2)
        x2=P2(1); y2=P2(2)
        x3=P3(1); y3=P3(2)
        MK(1,1)=x2-x1 ; MK(1,2)=x3-x1  
        MK(2,1)=y2-y1 ; MK(2,2)=y3-y1 
        area=ABS(MK(1,1)*MK(2,2)-MK(2,1)*MK(1,2))/2 
! Computation of the local error
          XG=SUM(mesh%POS(TRI(1:3),1))/3.0
          YG=SUM(mesh%POS(TRI(1:3),2))/3.0
          solexact=solexacte(XG,YG)
          solnum=SUM(U(TRI(1:3)))/3.0
          error=error+area*(solexact-solnum)**2
          valnorm=valnorm+area*solexact**2
      END DO
      err=sqrt(error)/sqrt(valnorm)
end subroutine Compute_Error1

!ccccccccccccccccccccccccccccc
function phi1(x,y) result(z)
!ccccccccccccccccccccccccccccc
      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: z
         z=1-x-y
      return
end function phi1

!ccccccccccccccccccccccccccccc
function phi2(x,y) result(z)
!ccccccccccccccccccccccccccccc
      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: z
         z=x
      return
end function phi2


!ccccccccccccccccccccccccccccc
function phi3(x,y) result(z)
!ccccccccccccccccccccccccccccc
      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: z
         z=y
      return
end function phi3

!ccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Compute_Error2(mesh,U,err)
!cccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      type(msh),intent(in) :: mesh
      real (kind=8), dimension(:), allocatable, intent(in) :: U
      real (kind=8), intent(out) :: err     

      real (kind=8) :: error, valnorm, area, errorK, valnormK
      integer :: K
      INTEGER, dimension(3) :: TRI
      REAL (kind=8), dimension(2,3) :: GRADPHIREF,GRADPHI
      REAL (kind=8), dimension(2) :: P1,P2,P3
      REAL (kind=8) :: x1,x2,x3,y1,y2,y3
      INTEGER :: k1,k2,k3
      REAL (kind=8), dimension(2,2) :: MK

      INTEGER :: NPOINTSINTG,I,J
      REAL (kind=8), dimension(3) :: POIDS
      REAL (kind=8), dimension(3) :: POINTS
      REAL(kind=8) :: XJ,ZI,solnumij,solexij,XX,YY

      NPOINTSINTG=3
      POIDS=(/ 5.0/18.0, 8.0/18.0, 5.0/18.0 /)
      POINTS=(/ 0.5*(1.0-sqrt(3.0/5.0)), 0.5, 0.5*(1.0+sqrt(3.0/5.0)) /)

      error=0.0
      valnorm=0.0

      DO K=1,mesh%nbTriangles
! Computation of the area of the triangle
        TRI=mesh%triangles(K,1:3)
        k1=TRI(1); k2=TRI(2); k3=TRI(3);
        P1=mesh%POS(k1,1:2);     P2=mesh%POS(k2,1:2);     P3=mesh%POS(k3,1:2)
        x1=P1(1); y1=P1(2)
        x2=P2(1); y2=P2(2)
        x3=P3(1); y3=P3(2)
        MK(1,1)=x2-x1 ; MK(1,2)=x3-x1  
        MK(2,1)=y2-y1 ; MK(2,2)=y3-y1 
        area=ABS(MK(1,1)*MK(2,2)-MK(2,1)*MK(1,2))/2 

        errorK=0.0
        valnormK=0.0
! Computation of the local error
          DO I=1,NPOINTSINTG
           DO J=1,NPOINTSINTG
            XJ=POINTS(J)
            ZI=POINTS(I)
            solnumij=  U(k1)*phi1(XJ,(1-XJ)*ZI) &
                     + U(k2)*phi2(XJ,(1-XJ)*ZI) &
                     + U(k3)*phi3(XJ,(1-XJ)*ZI) 
            XX=x1+MK(1,1)*XJ+MK(1,2)*((1-XJ)*ZI)
            YY=y1+MK(2,1)*XJ+MK(2,2)*((1-XJ)*ZI)
            solexij = solexacte(XX,YY)
            errorK= errorK+(1-XJ)*POIDS(I)*POIDS(J)*(solexij-solnumij)**2
            valnormK=valnormK+(1-XJ)*POIDS(I)*POIDS(J)*solexij**2
           END DO
         END DO
         error=error+2*area*errorK
         valnorm=valnorm+2*area*valnormK
      END DO
      err=sqrt(error)/sqrt(valnorm)
end subroutine Compute_Error2

end module EFLib
