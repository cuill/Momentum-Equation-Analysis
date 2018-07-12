!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_GCN(U,V,XC,YC,N,M,ART,NE,IEC,IENODE,NBE,DT,DT1,EGF,XIJC,YIJC,A1U,A2U, &
              DLTYC,DLTXC,ADVXU,ADVXV,ADVYU,ADVYV,VISCX,VISCY,DPBPX,DPBPY)

!==============================================================================!
! this subroutine calculate horizontal advective, barotropic pressure gradient,!
!  horizontal diffusion in x and y momentum equations except vertical diffusion!
!  terms for internal mode.
! ADVXU,ADVXV,ADVYU,ADVYV,VISCX,VISCY,DPBPX,DPBPY
! DT - depth at previous time step at nodes
! DT1 - depth at previous time step at elements
!==============================================================================!

   IMPLICIT NONE
   
   INTEGER :: N,M,NE
   INTEGER :: IEC(NE,2),IENODE(NE,2),NBE(N,3)
   REAL*8 :: XC(N),YC(N),U(N),V(N),ART(N),DT(M),EGF(M),DT1(N)
   REAL*8 :: XIJC(NE),YIJC(NE),DLTYC(NE),DLTXC(NE),A1U(N,4),A2U(N,4)
   REAL*8 :: ADVXU(N),ADVXV(N),ADVYU(N),ADVYV(N),VISCX(N),VISCY(N),DPBPX(N),DPBPY(N)
   REAL*8 :: CC_EHV(N)
   REAL*8 :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL*8 :: XADVU,XADVV,YADVU,YADVV,TXXIJ,TYYIJ,TXYIJ
   REAL*8 :: VISCOF,VISCOF1,VISCOF2
   REAL*8 :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL*8 :: DIJ,ELIJ,TMPA,TMPB,TMP
   REAL*8 :: HORCON,FACT,FM1,HPRNU,EXFLUX,EXFLUXU,EXFLUXV,GRAV
   REAL*8 :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,J
   
   HORCON = 4.000E-1
   FACT = 1.0d0
   FM1  = 0.0d0
   HPRNU   = 1.000d0
   GRAV = 9.8d0
   CC_EHV = 1.0d0

!-----Loop Over Edges and Accumulate Flux--------------------------------------!

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ=0.5d0*(DT(J1)+DT(J2))
     ELIJ=0.5d0*(EGF(J1)+EGF(J2))

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
    XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)

!     DO K=1,KBM1
       COFA1=A1U(IA,1)*U(IA)+A1U(IA,2)*U(K1)+A1U(IA,3)*U(K2)+A1U(IA,4)*U(K3)
       COFA2=A2U(IA,1)*U(IA)+A2U(IA,2)*U(K1)+A2U(IA,3)*U(K2)+A2U(IA,4)*U(K3)
       COFA5=A1U(IA,1)*V(IA)+A1U(IA,2)*V(K1)+A1U(IA,3)*V(K2)+A1U(IA,4)*V(K3)
       COFA6=A2U(IA,1)*V(IA)+A2U(IA,2)*V(K1)+A2U(IA,3)*V(K2)+A2U(IA,4)*V(K3)

       UIJ1=U(IA)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA)+COFA5*XIJA+COFA6*YIJA

       COFA3=A1U(IB,1)*U(IB)+A1U(IB,2)*U(K4)+A1U(IB,3)*U(K5)+A1U(IB,4)*U(K6)
       COFA4=A2U(IB,1)*U(IB)+A2U(IB,2)*U(K4)+A2U(IB,3)*U(K5)+A2U(IB,4)*U(K6)
       COFA7=A1U(IB,1)*V(IB)+A1U(IB,2)*V(K4)+A1U(IB,3)*V(K5)+A1U(IB,4)*V(K6)
       COFA8=A2U(IB,1)*V(IB)+A2U(IB,2)*V(K4)+A2U(IB,3)*V(K5)+A2U(IB,4)*V(K6)

       UIJ2=U(IB)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB)+COFA7*XIJB+COFA8*YIJB

       UIJ=0.5d0*(UIJ1+UIJ2)
       VIJ=0.5d0*(VIJ1+VIJ2)
! right-hand side
       EXFLUX = -DIJ*UIJ*DLTYC(I) + DIJ*VIJ*DLTXC(I)
       EXFLUXU = -DIJ*UIJ*DLTYC(I)
       EXFLUXV = DIJ*VIJ*DLTXC(I)

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5d0*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5d0*(COFA4+COFA7)**2)
!---huang add 9---
!       VISCOF=FACT*0.5d0*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=HORCON*(FACT*0.5d0*(VISCOF1*CC_EHV(IA)+VISCOF2*CC_EHV(IB)) + FM1*0.5d0*(CC_EHV(IA)+CC_EHV(IB)))/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5d0*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))
       VISCX(IA) = VISCX(IA) + FXX
       VISCY(IA) = VISCY(IA) + FYY
       VISCX(IB) = VISCX(IB) - FXX
       VISCY(IB) = VISCY(IB) - FYY

! At the beginning, I made a mistake by using EXFLUXU/EXFLUXV in function SIGN
       XADVU=EXFLUXU*((1.0d0-SIGN(1.0d0,EXFLUX))*UIJ2+(1.0d0+SIGN(1.0d0,EXFLUX))*UIJ1)*0.5d0
       XADVV=EXFLUXV*((1.0d0-SIGN(1.0d0,EXFLUX))*UIJ2+(1.0d0+SIGN(1.0d0,EXFLUX))*UIJ1)*0.5d0
       YADVU=EXFLUXU*((1.0d0-SIGN(1.0d0,EXFLUX))*VIJ2+(1.0d0+SIGN(1.0d0,EXFLUX))*VIJ1)*0.5d0
       YADVV=EXFLUXV*((1.0d0-SIGN(1.0d0,EXFLUX))*VIJ2+(1.0d0+SIGN(1.0d0,EXFLUX))*VIJ1)*0.5d0
       ADVXU(IA) = ADVXU(IA) + XADVU
       ADVXV(IA) = ADVXV(IA) + XADVV
       ADVYU(IA) = ADVYU(IA) + YADVU
       ADVYV(IA) = ADVYV(IA) + YADVV
       ADVXU(IB) = ADVXU(IB) - XADVU
       ADVXV(IB) = ADVXV(IB) - XADVV
       ADVYU(IB) = ADVYU(IB) - YADVU
       ADVYV(IB) = ADVYV(IB) - YADVV

! right-hand side
       DPBPX(IA)=DPBPX(IA)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
       DPBPY(IA)=DPBPY(IA)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
       DPBPX(IB)=DPBPX(IB)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
       DPBPY(IB)=DPBPY(IB)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
!     END DO ! loop k
   END DO  ! loop NE

   DO I=1,N
      ADVXU(I) = ADVXU(I)/ART(I)
      ADVXV(I) = ADVXV(I)/ART(I)
      ADVYU(I) = ADVYU(I)/ART(I)
      ADVYV(I) = ADVYV(I)/ART(I)
      VISCX(I) = VISCX(I)/ART(I)
      VISCY(I) = VISCY(I)/ART(I)
      DPBPX(I) = DPBPX(I)/ART(I)
      DPBPY(I) = DPBPY(I)/ART(I)
   END DO

   RETURN
   END SUBROUTINE ADV_UV_EDGE_GCN
!==============================================================================!
