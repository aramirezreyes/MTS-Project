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
  REAL (kind=8), dimension(3,3) :: ALOC_stiff, ALOC_mass
  REAL (kind=8), dimension(2,3) :: GRADPHIREF,GRADPHI
  REAL (kind=8), dimension(3,3) :: PHIREF
  INTEGER, dimension(3) :: TRI
  REAL (kind=8), dimension(2) :: P1,P2,P3
  REAL (kind=8) :: x1,x2,x3,y1,y2,y3,THETA,MU
  INTEGER :: k1,k2,k3
  REAL (kind=8), dimension(2,2) :: MK,MKT,INVMKT
  REAL (kind=8) :: DETMK 

! Reference Stiffness Matrix
!ccccccccccccccccccccccccccc
  GRADPHIREF=reshape((/-1,-1,1,0,0,1/),(/2,3/))
! Reference mass matrix
!cccccccccccccccccccccccccccc
  PHIREF = reshape((/1./12.,1./24.,1./24.,1./24.,1./12.,1./24.,1./24.,1./24.,1./12./),(/3,3/))  
! Nodes coordinates
!cccccccccccccccccc 

   TRI=mesh%triangles(k,1:3)
   k1=TRI(1); k2=TRI(2); k3=TRI(3);
   P1=mesh%POS(k1,1:2);     P2=mesh%POS(k2,1:2);     P3=mesh%POS(k3,1:2)

   x1=P1(1); y1=P1(2)
   x2=P2(1); y2=P2(2)
   x3=P3(1); y3=P3(2)

!Adjustements of the parameters THETA and MU
!ccccccccccccccccccccccccccccccccccccccccccc
   IF ((mesh%triangles(k,4).EQ.030)) then
      print*, "Area 3"
      MU = 0.5
      THETA = 1.4
   ELSE
      IF(mesh%triangles(k,4).EQ.010) then
      print*, "Area 1"
         THETA = 1
      ELSE IF(mesh%triangles(k,4).EQ.020) then
      print*, "Area 2"
         THETA = 1.2
      END IF
      MU = 0
   END IF


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
   ALOC_stiff=THETA * ABS(DETMK)/2*MATMUL(TRANSPOSE(GRADPHI),GRADPHI)
   ALOC_mass = MU * ABS(DETMK)*PHIREF
   
   ALOC = ALOC_stiff + ALOC_mass
  
! RHSLOC Computation
!cccccccccccccccccccc
!CC is this approximation correct????????????????????????????????????????????????????????????
   RHSLOC=-3*MU*ABS(DETMK)/2
   
  end subroutine Compute_LocalStiffness_Matrix_And_LocalRHS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  Compute_Local_Border_Terms(mesh,K,BLOC)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE

  type(msh),intent(in) :: mesh
  INTEGER, intent(in)  :: K !number of border segment

  REAL (kind=8), dimension(2,2), intent(out) :: BLOC
  REAL (kind=8), dimension(2) :: P1,P2
  REAL (kind=8) :: x1,x2,y1,y2,MU,e
  INTEGER, dimension(2) :: SEG
  INTEGER :: k1,k2,label
  BLOC=reshape((/2,1,1,2/),(/2,2/))
  SEG = mesh%lines(K,1:2)
  k1 = SEG(1); k2 = SEG(2)
  P1 = mesh%pos(k1,1:2)
  P2 = mesh%pos(k2,1:2)
  x1 = P1(1) ; y1 = P1(2) ; x2 = P2(1) ; y2 = P2(2)
  LABEL = mesh%lines(K,3)
  IF (LABEL .EQ. 200) THEN
     MU = 0.5
  ELSE
     MU = 0
  END IF
  e = SQRT((x2-x1)**2 + (y2-y1)**2)
  BLOC = (MU*e/12)*BLOC
  end subroutine Compute_Local_Border_Terms



  
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
!cccccccccccccccccccccccccc Changes made
 allocate(V1(NN),V2(NN))
  V1=0.0 ; V2=0.0
  cl=0.
  Zone_Phys_Bord=000
  do j=1,mesh%nbLines
     !! Identification of the boundary node to treat
     lig=boundary_nodes(j,1)
     type=boundary_nodes(j,2)
     if (type.eq.Zone_Phys_Bord) then
       val_CL=CL
     
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
    end if !ccccccccccccccccccccccccccccccccccccccccc end of changes made
     end do

  deallocate(boundary_nodes,vu,v1,v2)

end subroutine Boundary_Conditions

!ccccccccccccccccccccccccccccc
function source(x,y) result(z)
!ccccccccccccccccccccccccccccc

      IMPLICIT NONE

      REAL(kind=8), intent(in) :: x,y
      REAL(kind=8) :: r,z

      r=sqrt(x**2.+y**2.)
      z=0.
      if ((r>1.5).and.(r<2)) then
         z=z+10000*((r-1.5)**2.)*(r-2)**2.
      end if
      if (r<1.) then
         z=z-250*(r**2.)*((r-1)**2.)
      end if
      return

end function source

end module EFLib
