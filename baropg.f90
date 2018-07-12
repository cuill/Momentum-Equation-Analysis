!==============================================================================|
!     CALCULATE THE BAROCLINIC PRESSURE GRADIENT IN SIGMA COORDINATES          |
!==============================================================================|

   SUBROUTINE BAROPG(RHO1,RHO,VX,VY,ZZ,NV,N,M,KBM1,ART,DT,DT1,DPBCX,DPBCY)

   IMPLICIT NONE
   INTEGER :: N,M,KBM1
   INTEGER :: NV(N,3)
   REAL*8 :: RHO1(M,KBM1),RHO(N,KBM1),VX(M),VY(M),ZZ(KBM1)
   REAL*8 :: ART(N),DT(M),DT1(N),DPBCX(N,KBM1),DPBCY(N,KBM1)
   REAL*8 :: RIJK(N,3,KBM1), DRIJK1(N,3,KBM1), DRIJK2(N,KBM1)
   REAL*8 :: TEMP,RAMP1,DIJ,DRHO1,DRHO2,GRAV
   INTEGER  :: I,K,J,J1,J2,IJK

!----------SUBTRACT MEAN DENSITY TO MINIMIZE ROUNDOFF ERROR--------------------!

!   RHO1(:,1:KBM1) = RHO1(:,1:KBM1) - RMEAN1(:,1:KBM1)
!   RHO = RHO - RMEAN 

!----------INITIALIZE ARRAYS---------------------------------------------------!

 !  RMEAN(0,:) = 0.0d0
 !  RHO(0,:)   = 0.0d0
   RIJK       = 0.0d0
   DRIJK1     = 0.0d0
   DRIJK2     = 0.0d0
   GRAV       = 9.8d0

!----------CALCULATE AVERAGE DENSITY ON EACH EDGE------------------------------!

   DO K=1,KBM1
     DO I=1,N
       DO J=1,3
         J1=J+1-INT((J+1)/4)*3
         J2=J+2-INT((J+2)/4)*3
         RIJK(I,J,K)=0.5d0*(RHO1(NV(I,J1),K)+RHO1(NV(I,J2),K))
       END DO
     END DO
   END DO

   DO I=1,N
     DO J=1,3
       DRIJK1(I,J,1)=RIJK(I,J,1)*(-ZZ(1))
       DO K=2,KBM1
         DRIJK1(I,J,K)=0.5d0*(RIJK(I,J,K-1)+RIJK(I,J,K))*(ZZ(K-1)-ZZ(K))
         DRIJK1(I,J,K)=DRIJK1(I,J,K)+DRIJK1(I,J,K-1)
       END DO
     END DO
   END DO

   DO I=1,N
     DRIJK2(I,1)=0.0d0
     DO K=2,KBM1
       DRIJK2(I,K)=0.5d0*(ZZ(K-1)+ZZ(K))*(RHO(I,K)-RHO(I,K-1))
       DRIJK2(I,K)=DRIJK2(I,K-1)+DRIJK2(I,K)
     END DO
   END DO

   DO I = 1, N
     DO K=1,KBM1
        DO J = 1, 3
          J1=J+1-INT((J+1)/4)*3
          J2=J+2-INT((J+2)/4)*3
     !     IJK=NBE(I,J)
          DIJ=0.5d0*(DT(NV(I,J1))+DT(NV(I,J2)))

! right-hand side 
          DRHO1=(VY(NV(I,J1))-VY(NV(I,J2)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VY(NV(I,J1))-VY(NV(I,J2)))*DIJ*DRIJK2(I,K)
          DPBCX(I,K)=DPBCX(I,K)+DRHO1+DRHO2

	  DRHO1=(VX(NV(I,J2))-VX(NV(I,J1)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VX(NV(I,J2))-VX(NV(I,J1)))*DIJ*DRIJK2(I,K)
          DPBCY(I,K)=DPBCY(I,K)+DRHO1+DRHO2

       END DO
     END DO
   END DO

!----------MULTIPLY BY GRAVITY AND ELEMENT DEPTH AND DIVIDED BY ART-------------------------------!

   DO K=1,KBM1
     DPBCX(:,K)=DPBCX(:,K)*DT1(:)*GRAV/ART(:)
     DPBCY(:,K)=DPBCY(:,K)*DT1(:)*GRAV/ART(:)
   END DO

   END SUBROUTINE BAROPG
