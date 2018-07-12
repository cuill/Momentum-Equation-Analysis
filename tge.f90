   SUBROUTINE TRIANGLE_GRID_EDGE(vx,vy,xc,yc,nv,N,NCT,M,NE,NTSN,NBSN,NTVE,NBVE, &
              NBVT,IEC,IENODE,NBE,XIJC,YIJC,DLTXC,DLTYC)

!==============================================================================|
!   FIND NEIGHBORING ELEMENTS, MARK BOUNDARY NODES AND ELEMENTS                |
!									       |
!   NBE(N,3) :  NBE(I,J) = IDENTITY OF NBOR ELMNT TO TRI I ON EDGE J           |
!   IBCE(N)  :  DESCRIBED IN SUBROUTINE HEADING			               |	
!   ISONB(M):  DESCRIBED IN SUBROUTINE HEADING			               |	
!==============================================================================|
   IMPLICIT NONE
   INTEGER N,M,NCT,N1,N2,N3,J1,J2,J3,I,J,NC,MX_NBR_ELEM
   INTEGER NBE(N,3),NV(N,3),NTSN(M),NBSN(M,11),NTVE(M),NBVE(M,9),NBVT(M,9),ISONB(M)
   INTEGER, ALLOCATABLE :: ISBCE(:)
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,TEMP2,NB_TMP,ISET
   INTEGER II,JJ,NTMP,NCNT,INEY,JN,NE,NCV,NCV_I
   INTEGER ITMP1,ITMP2,ITMP3,JJB,IBCETMP,NCTMP,NCETMP,NPT
   INTEGER, ALLOCATABLE :: CELLS(:,:),CELLCNT(:)
   REAL*8  vx(M),vy(M),xc(N),yc(N)
   INTEGER, ALLOCATABLE :: ISBC(:)
   INTEGER :: NIEC(NCT,2)
   INTEGER :: NTRG(NCT)
   INTEGER :: IEC(NE,2)
   INTEGER :: IENODE(NE,2)
   REAL*8 :: DLTXC(NE)
   REAL*8 :: DLTYC(NE)
   REAL*8 :: DLTXYC(NE)
   REAL*8 :: DLTXE(NCT)
   REAL*8 :: DLTYE(NCT)
   REAL*8 :: DLTXYE(NCT)
   REAL*8,ALLOCATABLE :: SITAC(:)
   REAL*8 :: SITAE(NCT)
   REAL*8 :: XIJC(NE)
   REAL*8 :: YIJC(NE)
   REAL*8 :: XIJE(NCT,2)
   REAL*8 :: YIJE(NCT,2)
   REAL*8  DTMP

!----------------------------INITIALIZE----------------------------------------!
   
   ALLOCATE(ISBCE(N))
   ISBCE = 0
!   ISONB = 0
   NBE   = 0

!
!----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!
!
   ALLOCATE(CELLS(M,50)) ; CELLS = 0
   ALLOCATE(CELLCNT(M))  ; CELLCNT = 0
   print*,'I am here 1!'
   DO I=1,N
     N1 = NV(I,1) ; CELLCNT(N1) = CELLCNT(N1)+1
     N2 = NV(I,2) ; CELLCNT(N2) = CELLCNT(N2)+1
     N3 = NV(I,3) ; CELLCNT(N3) = CELLCNT(N3)+1
     CELLS(NV(I,1),CELLCNT(N1)) = I
     CELLS(NV(I,2),CELLCNT(N2)) = I
     CELLS(NV(I,3),CELLCNT(N3)) = I
   END DO
   if(maxval(cellcnt) > 50)write(*,*)'bad',maxval(cellcnt)
   DO I=1,N
     N1 = NV(I,1)
     N2 = NV(I,2)
     N3 = NV(I,3)
     DO J1 = 1,CELLCNT(N1) 
     DO J2 = 1,CELLCNT(N2) 
       IF((CELLS(N1,J1) == CELLS(N2,J2)).AND. CELLS(N1,J1) /= I)NBE(I,3) = CELLS(N1,J1)
     END DO
     END DO
     DO J2 = 1,CELLCNT(N2) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N2,J2) == CELLS(N3,J3)).AND. CELLS(N2,J2) /= I)NBE(I,1) = CELLS(N2,J2)
     END DO
     END DO
     DO J1 = 1,CELLCNT(N1) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N1,J1) == CELLS(N3,J3)).AND. CELLS(N1,J1) /= I)NBE(I,2) = CELLS(N3,J3)
     END DO
     END DO
   END DO
   DEALLOCATE(CELLS,CELLCNT)
   print*,'I am here 2!'
!
!--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
!
!   NFLAG = 0
!   DO I=1,N
!     IF(SUM(NBE(I,1:3))==0)THEN 
!       NFLAG = 1
!       WRITE(*,*)'ELEMENT ',I,' AT ',XC(I),YC(I),' HAS NO NEIGHBORS'
!       CALL STOP
!     END IF
!   END DO
!   IF(NFLAG == 1) CALL STOP
     
!
!----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
!

   DO I=1,N 
     IF(MIN(NBE(I,1),NBE(I,2),NBE(I,3))==0)THEN    !!ELEMENT ON BOUNDARY
       ISBCE(I) = 1
       IF(NBE(I,1) == 0)THEN 
         ISONB(NV(I,2)) = 1 ; ISONB(NV(I,3)) = 1
       END IF
       IF(NBE(I,2) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,3)) = 1
       END IF
       IF(NBE(I,3) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,2)) = 1
       END IF
     END IF
   END DO
        
   print*,'I am here 3!'
!==============================================================================|
!             DEFINE NTVE, NBVE, NBVT                                          !
!                                                                              !
! ntve(1:m):           total number of the surrounding triangles               !
!                      connected to the given node                             !                                        
! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
!                      triangles with a common node (counted clockwise)        !
! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
!                      each individual surrounding triangle(counted            !
!                      clockwise)                                              !
! ntsn(1:m):           total number of surrounding nodes                       !
! nbsn(1:m, ntsn):     the identification number of surrounding nodes          !
!                      (counted clockwise)                                     !
! nbse(1:m,2*ntsn):    the identification number of control volume s           !
!                      edges between two neighbor nodes                        !
!==============================================================================|

!
!----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
!
! To speed up code, define MX_NBR_ELEM equal to 8 instead of searching.
   MX_NBR_ELEM =8
!   DO I=1,M
!     NCNT = 0
!     DO J=1,N
!       IF( FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0) &
!         NCNT = NCNT + 1
!     END DO
!     MX_NBR_ELEM = MAX(MX_NBR_ELEM,NCNT)
!   END DO
!   WRITE(IPT,*) 'MAXIMUM NUMBER OF NEIGHBOR ELEMENTS',MX_NBR_ELEM

!
!----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
! 
!   ALLOCATE(NBVE(M,MX_NBR_ELEM+1))
!   ALLOCATE(NBVT(M,MX_NBR_ELEM+1))
!   ALLOCATE(NBSN(M,MX_NBR_ELEM+3)) !!MHB
!
!--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
!--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
!--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
!
       
   DO I=1,M
     NCNT=0
     DO J=1,N
!       IF( (NV(J,1)-I) == 0 .OR.  (NV(J,2)-I) == 0 .OR. (NV(J,3)-I) == 0)THEN 
        IF (FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0)THEN
         NCNT = NCNT+1
         NBVE(I,NCNT)=J
         IF((NV(J,1)-I) == 0) NBVT(I,NCNT)=1
         IF((NV(J,2)-I) == 0) NBVT(I,NCNT)=2
         IF((NV(J,3)-I) == 0) NBVT(I,NCNT)=3
       END IF
     ENDDO
     NTVE(I)=NCNT
   ENDDO

!
!--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
!--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
!--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!
!
   ALLOCATE(NB_TMP(M,MX_NBR_ELEM+1))
   DO I=1,M
     IF(ISONB(I) == 0) THEN
       NB_TMP(1,1)=NBVE(I,1)
       NB_TMP(1,2)=NBVT(I,1)
       DO J=2,NTVE(I)+1
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=2,NTVE(I)+1
         NBVE(I,J)=NB_TMP(J,1)
       ENDDO

       DO J=2,NTVE(I)+1
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NTMP=NTVE(I)+1
!       IF(NBVE(I,1) /= NBVE(I,NTMP)) THEN
!          PRINT*, I,'NBVE(I) NOT CORRECT!!'
!          CALL STOP
!       ENDIF
!       IF(NBVT(I,1) /= NBVT(I,NTMP)) THEN
!          PRINT*, I,'NBVT(I) NOT CORRECT!!'
!          CALL STOP
!       END IF

       NTSN(I)=NTVE(I)

       DO J=1,NTSN(I)
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       NTSN(I)=NTSN(I)+1
       NBSN(I,NTSN(I))=NBSN(I,1)

     ELSE 
           JJB=0

       DO J=1,NTVE(I)
         JJ=NBVT(I,J)
         IF(NBE(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3) == 0) THEN
           JJB=JJB+1
           NB_TMP(JJB,1)=NBVE(I,J)
           NB_TMP(JJB,2)=NBVT(I,J)
         END IF
       ENDDO

!       IF(JJB /= 1) THEN
!         PRINT*, 'ERROR IN ISONB !,I,J', I,J
!         PAUSE
!       END IF

       DO J=2,NTVE(I)
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=1,NTVE(I)
         NBVE(I,J)=NB_TMP(J,1)
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NBVE(I,NTVE(I)+1)=0
       NTSN(I)=NTVE(I)+1
       NBSN(I,1)=I

       DO J=1,NTSN(I)-1
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J+1)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       J=NTSN(I)
       II=NBVE(I,J-1)
       JJ=NBVT(I,J-1)
       NBSN(I,J+1)=NV(II,JJ+2-INT((JJ+2)/4)*3)
       NTSN(I)=NTSN(I)+2
       NBSN(I,NTSN(I))=I
     END IF
   END DO
   DEALLOCATE(NB_TMP)
!   IF(MX_NBR_ELEM+3 -MAXVAL(NTSN) < 0)THEN
!      WRITE(*,*)'CHECK NTSN/NBSN',MAXVAL(NTSN),MX_NBR_ELEM+3
!      CALL STOP
!   END IF

!==============================================================================!
!  Define the parameters of each triangular edge                               !
!                                                                              !
!  ne           :    number of unique element edges                            !
!  iec(1:ne,1:2):    counting number identifying two connected cells           !
!  isbc(1:ne):       0: triangle s edge in the interior                        !
!                    1: triangle s edge on the boundary                        !
!  ienode(1:ne,1:2): the identification number of two end points of a          !
!                    edge                                                      !
!  xijc(1:ne):       the x-coordinate location of the middle points            !
!                    of a edge                                                 !
!  yijc(1:ne):       the y-coordinate location of the middle points            !
!                    of a edge                                                 !
!  dltxyc(1:ne):     length of the edge                                        !
!  dltxc(1:ne):      vx(ienode(i,2))-vx(idnode(i,1))                           !
!  dltyc(1:ne):      vy(ienode(i,2))-vy(idnode(i,1))                           !
!  sitac(1:ne):      arctg(dltyc,dltxc)                                        !
!==============================================================================!

   print*,'NE1=',NE
   ALLOCATE(ISET(N,3),TEMP((N)*3,2),TEMP2((N)*3,2))
   ISET = 0
   NE = 0
   TEMP = 0
   TEMP2 = 0
   DO I=1,N
     DO J=1,3
       IF(ISET(I,J) == 0)THEN
         NE   = NE + 1
         INEY = NBE(I,J)
         ISET(I,J) = 1
         DO JN=1,3
           IF(I == NBE(INEY,JN)) ISET(INEY,JN) = 1
         END DO
         TEMP(NE,1) = I ; TEMP(NE,2) = INEY
         TEMP2(NE,1) = NV(I,J+1-INT((J+1)/4)*3)
         TEMP2(NE,2) = NV(I,J+2-INT((J+2)/4)*3)
       END IF
     END DO
   END DO
   print*,'NE=',NE
   DEALLOCATE(ISET)
     
!--ALLOCATE ARRAYS REQUIRING NUMBER OF EDGES-----------------------------------!
!
!   ALLOCATE(IEC(NE,2))
!   ALLOCATE(IENODE(NE,2))
!   ALLOCATE(XIJC(NE))
!   ALLOCATE(YIJC(NE))
!   ALLOCATE(DLTXYC(NE))
!   ALLOCATE(DLTXC(NE))
!   ALLOCATE(DLTYC(NE))
   ALLOCATE(SITAC(NE))
   ALLOCATE(ISBC(NE))

   IEC(:,1) = TEMP(1:NE,1)
   IEC(:,2) = TEMP(1:NE,2)
   IENODE(:,1) = TEMP2(1:NE,1)
   IENODE(:,2) = TEMP2(1:NE,2)


   DEALLOCATE(TEMP,TEMP2)

!
!------MARK ELEMENT EDGES THAT ARE ON THE BOUNDARY-----------------------------!
!
   ISBC = 0
   DO I=1,NE
     IF((IEC(I,1) == 0) .OR. (IEC(I,2) == 0)) ISBC(I) = 1
   END DO

!
!------CALCULATE ELEMENT EDGE METRICS------------------------------------------!
!

   DO I=1,NE
     DLTXC(I) =  VX(IENODE(I,2))-VX(IENODE(I,1))
     DLTYC(I) =  VY(IENODE(I,2))-VY(IENODE(I,1))
     XIJC(I)  = (VX(IENODE(I,1))+VX(IENODE(I,2)))/2.0
     YIJC(I)  = (VY(IENODE(I,1))+VY(IENODE(I,2)))/2.0
     DLTXYC(I)= SQRT(DLTXC(I)**2+DLTYC(I)**2)
     SITAC(I) = ATAN2(DLTYC(I),DLTXC(I))
   END DO

!
!----TRAVERSE  BOUNDARY NODE NUMBERS AND SET ISONB(NODE)=2---------------------!
!
   DO I=1,214
     ISONB(I)=2
   ENDDO

!
!----DETERMINE IF ELEMENT IS ON OPEN BOUNDARY (CONTAINS EDGE ON OPEN BOUNDARY)-!
!
   IBCETMP=0
   DO I=1,N
     ITMP1=ISONB(NV(I,1))
     ITMP2=ISONB(NV(I,2))
     ITMP3=ISONB(NV(I,3))
!JQI< 07/23/04
!JQI     IF(SUM(ISONB(NV(I,1:3))) >= 4) THEN
     IF(SUM(ISONB(NV(I,1:3))) == 4) THEN
       ISBCE(I)=2
       IBCETMP =IBCETMP+1
     ELSE IF(SUM(ISONB(NV(I,1:3))) > 4) THEN
       PRINT*,'SORRY, THE BOUNDARY CELL',I,'IS NOT GOOD FOR MODEL.'
     END IF
   END DO

   DO I=1,N
     IF((NBE(I,1)+NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,2) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
   ENDDO

!==============================================================================!
!  xije(1:nc,1:2):  the x coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  yije(1:nc,1:2):  the y coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  niec(1:nc,1:2):  the counting number of left and right nodes                !
!                   conected to this control volumes edge from                 !
!                   starting point to ending point                             !
!  dltxe(1:nc):     the x distance of individual edges                         !
!  dltye(1:nc)      the y distance of individual edges                         !
!  dltxye(1:nc):    the length of individual edges                             !
!  ntrg(1:nc)  :    element associated with this control volume edge           !
!==============================================================================!
   NCTMP  = 0
   NCETMP = 0

   DO I=1,NE
     IF(ISBC(I) == 0) THEN
       IF(IEC(I,1) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       XIJE(NPT,1) = XC(IEC(I,1))
       YIJE(NPT,1) = YC(IEC(I,1))
       XIJE(NPT,2) = XIJC(I)
       YIJE(NPT,2) = YIJC(I)
       NIEC(NPT,1) = IENODE(I,1)
       NIEC(NPT,2) = IENODE(I,2)
       NTRG(NPT)   = IEC(I,1)
       DLTXE(NPT)  = XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)  = YIJE(NPT,2)-YIJE(NPT,1)

       DTMP        = DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT) = SQRT(DTMP)
       SITAE(NPT)  = ATAN2(DLTYE(NPT),DLTXE(NPT))

       IF(IEC(I,2) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       XIJE(NPT,1)=XC(IEC(I,2))
       YIJE(NPT,1)=YC(IEC(I,2))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,2)
       NIEC(NPT,2)=IENODE(I,1)
       NTRG(NPT)=IEC(I,2)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))

     ELSE IF(ISBC(I) == 1) THEN
       IF(IEC(I,1) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       IF(IEC(I,1) == 0) THEN
         PRINT*, I,'IEC(I,1)===0'  
!         CALL PSTOP
       END IF
       XIJE(NPT,1)=XC(IEC(I,1))
       YIJE(NPT,1)=YC(IEC(I,1))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,1)
       NIEC(NPT,2)=IENODE(I,2)
       NTRG(NPT)=IEC(I,1)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))
     ELSE
       WRITE(*,*) 'ISBC(I) NOT CORRECT, I==',I
!       CALL PSTOP
     END IF
   ENDDO

   NCV_I  = NCTMP
   NCV    = NCETMP+NCTMP

   IF(NCV /= 3*(N)) THEN
     PRINT*,'NCV IS NOT CORRECT, PLEASE CHECK THE SETUP'
!     CALL PSTOP
   END IF
   IF(NCV_I /= 3*N) THEN
     PRINT*,'NCV_I IS NOT CORRECT, PLEASE CHECK THE SETUP'
!      CALL PSTOP
   END IF

   DO I=1,NCV_I
     IF(NIEC(I,1) > M .OR. NIEC(I,2) > M)THEN
       write(*,*)'problemas',niec(i,1),niec(i,2),m
!       CALL PSTOP
     END IF
   END DO

!   open(22,file='tge_infor.dat')
!   do I=1,NCTMP
!      write(22,'(6(f15.4,2x),3(i8,2x))') dltxe(i),dltye(i), &
!           xije(i,1),xije(i,2),yije(i,1),yije(i,2),niec(i,1),niec(i,2),ntrg(i)
!   end do
!   close(22)

   RETURN
   END SUBROUTINE TRIANGLE_GRID_EDGE
!==============================================================================!
