#include "CTRL_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C--   File ctrl_map_genarr.F: Generic arrays control vector
C--    Contents:
C--    o CTRL_MAP_GENARR2D
C--    o CTRL_MAP_GENARR3D

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C     !ROUTINE: CTRL_MAP_GENARR2D
C     !INTERFACE:
      SUBROUTINE CTRL_MAP_GENARR2D( fld, iarr, myThid )

C     !DESCRIPTION: \bv
C     *=================================================================
C     | SUBROUTINE CTRL_MAP_GENARR2D
C     | Add the generic 2D-arrays of the
C     | control vector to the model state and update the tile halos.
C     | The control vector is defined in the header file "CTRL.h".
C     *=================================================================
C     \ev

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

#include "CTRL_SIZE.h"
#include "CTRL.h"
#include "OPTIMCYCLE.h"
#include "CTRL_GENARR.h"
#include "CTRL_DUMMY.h"
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES 
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h" /* for scaling factors */
#endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
      _RL     fld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER iarr
      INTEGER myThid

#ifdef ALLOW_GENARR2D_CONTROL
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
      INTEGER bi,bj
      INTEGER i,j
      INTEGER numsmo, k2
      LOGICAL dowc01
      LOGICAL dosmooth
      LOGICAL doscaling
      LOGICAL dolog10ctrl
      _RL log10initval
      _RL ln10
      _RL xx_gen  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      CHARACTER*(MAX_LEN_FNAM) fnamegenIn
      CHARACTER*(MAX_LEN_FNAM) fnamegenOut
      CHARACTER*(MAX_LEN_FNAM) fnamebase
#ifdef MULTISCALE_COUPLING_TAPES 
      character*(80) fnamegenInlo
      character*(80) fnamegenOutlo
#endif
      INTEGER ilgen
      LOGICAL doglobalread
      LOGICAL ladinit
      _RS mask2D (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER ilDir
#ifdef MULTISCALE_COUPLING_TAPES 
catn: first attempt, hard-coding 365
      _RL     fld_lo(1:sNx/MSfacX,1:sNy/MSfacY,nSx,nSy)
      _RL     xx_gen_lo(1:sNx/MSfacX,1:sNy/MSfacY,nSx,nSy)
      _RL     wgenarr2d_lo(1:sNx/MSfacX,1:sNy/MSfacY,nSx,nSy)
      integer nArr,tmpcount,tmpind1,tmpind2,k,l,m,helpi,helpj,il
#endif
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('CTRL_MAP_GENARR2D', myThid )
#endif
C--   Now, read the control vector.
      doglobalread = .FALSE.
      ladinit      = .FALSE.

#ifndef ALLOW_TAPENADE
      CALL CTRL_ASSIGN(xx_gen, 1, zeroRL, myThid)
#endif

      dosmooth=.FALSE.
      dowc01  = .FALSE.
      doscaling=.TRUE.
      dolog10ctrl=.FALSE.
      log10initval=0. _d 0
      ln10 = LOG(10. _d 0)

#ifdef MULTISCALE_COUPLING_TAPES
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j = 1,sNy/MSfacY
         DO i = 1,sNx/MSfacY
          fld_lo(i,j,bi,bj)=0. _d 0
          xx_gen_lo(i,j,bi,bj)=0. _d 0
          wgenarr2d_lo(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* MULTISCALE_COUPLING_TAPES */

      numsmo = 1
      DO k2 = 1, maxCtrlProc
        IF ( xx_genarr2d_preproc(k2,iarr).EQ.'WC01' ) THEN
          dowc01 = .TRUE.
          IF (xx_genarr2d_preproc_i(k2,iarr).NE.0)
     &      numsmo = xx_genarr2d_preproc_i(k2,iarr)
        ENDIF
        IF ( (.NOT.dowc01).AND.
     &       (xx_genarr2d_preproc(k2,iarr).EQ.'smooth') ) THEN
          dosmooth = .TRUE.
          IF (xx_genarr2d_preproc_i(k2,iarr).NE.0)
     &      numsmo = xx_genarr2d_preproc_i(k2,iarr)
        ENDIF
        IF ( xx_genarr2d_preproc(k2,iarr).EQ.'noscaling' ) THEN
          doscaling = .FALSE.
        ENDIF
        IF (xx_genarr2d_preproc_c(k2,iarr).EQ.'log10ctrl') THEN
C             fld := log10(xx_gen)
C               with initial guess for xx_gen set to log10initval
C               passed to data.ctrl, (default = 0.0)
          dolog10ctrl=.TRUE.
          log10initval = xx_genarr2d_preproc_r(k2,iarr)
        ENDIF
      ENDDO

      fnamebase = xx_genarr2d_file(iarr)
      ilgen = ILNBLNK( fnamebase )
C     Find ctrlDir (w/o trailing blanks) length
      ilDir = ilnblnk(ctrlDir)
      WRITE(fnamegenIn,'(2A,I10.10)')
     & ctrlDir(1:ilDir)//fnamebase(1:ilgen),'.',optimcycle
      WRITE(fnamegenOut,'(2A,I10.10)')
     & ctrlDir(1:ilDir)//fnamebase(1:ilgen),'.effective.',optimcycle
      write(*,'(A)') 'genarr2d fnamegenIn: ', fnamegenIn(1:il)
      il=ilnblnk( fnamegenOut )
      write(*,'(A)') 'genarr2d fnamegenOut: ', fnamegenOut(1:il)
#ifdef MULTISCALE_COUPLING_TAPES
      write(fnamegenInlo(1:80),'(2a,i10.10)')
     &     fnamebase(1:ilgen),'_lo.',optimcycle
      write(fnamegenOutlo(1:80),'(2a,i10.10)')
     &     fnamebase(1:ilgen),'_lo.effective.',optimcycle
      il=ilnblnk( fnamegenInlo )
      write(*,'(A)') 'genarr2d fnamegenInlo: ', fnamegenInlo(1:il)
      il=ilnblnk( fnamegenOutlo )
      write(*,'(A)') 'genarr2d fnamegenOutlo: ', fnamegenOutlo(1:il)
#endif /* MULTISCALE_COUPLING_TAPES */

      CALL READ_REC_3D_RL( xx_genarr2d_weight(iarr), ctrlprec, 1,
     &          wgenarr2d(1-OLx,1-OLy,1,1,iarr), 1, 1, myThid )

#ifdef ALLOW_AUTODIFF
      CALL ACTIVE_READ_XY( fnamegenIn, xx_gen, 1, doglobalread,
     &     ladinit, optimcycle, myThid, xx_genarr2d_dummy(iarr) )
#else
      CALL READ_REC_XY_RL( fnamegenIn, xx_gen, 1, 1, myThid)
#endif

C --- Get appropriate mask for ctrl variable
      CALL CTRL_GET_MASK2D(xx_genarr2d_file(iarr), mask2D, myThid)

#ifdef MULTISCALE_COUPLING_TAPES
      write(*,*) 'genarr2d: i,j,tmpcount,xx_gen_lo'
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        do j=1,sNy/MSfacY            !1:16
         do i=1,sNx/MSfacX           !1:20
          tmpcount=0
          do l=1,MSfacY              !1:4
           do m=1,MSfacX             !1:4
            tmpind1 = (i-1)*MSfacX+m !1:64
            tmpind2 = (j-1)*MSfacY+l !1:80
            if(mask2D(tmpind1,tmpind2,bi,bj).ne. 0.) then
             tmpcount=tmpcount+1
             xx_gen_lo(i,j,bi,bj)=
     &       xx_gen_lo(i,j,bi,bj)+
     &       xx_gen(tmpind1,tmpind2,bi,bj)

             fld_lo(i,j,bi,bj)=
     &       fld_lo(i,j,bi,bj)+
     &       fld(tmpind1,tmpind2,bi,bj)
            endif
           enddo   !m
          enddo    !l
          if(tmpcount .gt. 0) then
           xx_gen_lo(i,j,bi,bj)=xx_gen_lo(i,j,bi,bj)/tmpcount
           fld_lo(i,j,bi,bj)=fld_lo(i,j,bi,bj)/tmpcount
          endif
          write(*,*) i,j,tmpcount,xx_gen_lo(i,j,bi,bj),fld_lo(i,j,bi,bj)
         enddo     !i
        enddo      !j
       ENDDO        !bi
      ENDDO         !bj
#endif /* MULTISCALE_COUPLING_TAPES */

C --- Do any smoothing
#ifdef ALLOW_SMOOTH
      IF (useSMOOTH) THEN
       IF (dowc01) CALL SMOOTH_CORREL2D(xx_gen,mask2D,numsmo,myThid)
       IF (dosmooth) CALL SMOOTH2D(xx_gen,mask2D,numsmo,myThid)
      ENDIF
#endif

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        IF (doscaling) THEN
         DO j = 1,sNy
          DO i = 1,sNx
C scale param adjustment
           IF ( wgenarr2d(i,j,bi,bj,iarr).GT.0. ) THEN
            xx_gen(i,j,bi,bj) = xx_gen(i,j,bi,bj)
     &               / SQRT( wgenarr2d(i,j,bi,bj,iarr) )
            IF (dolog10ctrl) THEN
             xx_gen(i,j,bi,bj) = xx_gen(i,j,bi,bj) + log10initval
             xx_gen(i,j,bi,bj) = EXP(ln10 * xx_gen(i,j,bi,bj))
            ENDIF ! dolog10ctrls
           ELSE
            xx_gen(i,j,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDIF ! doscaling

#ifdef MULTISCALE_COUPLING_TAPES
        DO j=1,sNy/MSfacY            !1:16
         DO i=1,sNx/MSfacX           !1:20
          IF (doscaling) then
           tmpcount=0
           DO l=1,MSfacY              !1:4
            DO m=1,MSfacX             !1:4
             tmpind1 = (i-1)*MSfacX+m   !1:64
             tmpind2 = (j-1)*MSfacY+l   !1:80
             IF ((mask2D(tmpind1,tmpind2,bi,bj).NE.0.).AND.
     &            (wgenarr2d(tmpind1,tmpind2,bi,bj,iarr).GT.0.))
     &        tmpcount=tmpcount+1
            ENDDO   !m
           ENDDO    !l
           IF(tmpcount.GT.0) THEN   !pick at least 4 pts 
            DO l=1,MSfacY                !1:4
             DO m=1,MSfacX               !1:4
              tmpind1 = (i-1)*MSfacX+m   !1:64
              tmpind2 = (j-1)*MSfacY+l   !1:80
              wgenarr2d_lo(i,j,bi,bj)  =
     &          max( wgenarr2d_lo(i,j,bi,bj),
     &          wgenarr2d(tmpind1,tmpind2,bi,bj,iarr) )
             ENDDO   !m
            ENDDO    !l
            fld_lo(i,j,bi,bj)=fld_lo(i,j,bi,bj)+
     &        xx_gen_lo(i,j,bi,bj)
     &              /sqrt(wgenarr2d_lo(i,j,bi,bj))
            write(*,*) i,j,tmpcount,fld_lo(i,j,bi,bj)
           ELSE
            fld_lo(i,j,bi,bj)= 0. _d 0
           ENDIF    !tmpcount gt 0
          ELSE
           fld_lo(i,j,bi,bj)=fld_lo(i,j,bi,bj)+xx_gen_lo(i,j,bi,bj)
          ENDIF     !doscaling
         ENDDO      !i
        ENDDO       !j
#endif /* MULTISCALE_COUPLING_TAPES */        

C add to model parameter
C or in case of log10ctrl, fld = 10^(xx_gen)
        IF (dolog10ctrl) THEN
         DO j = 1,sNy
          DO i = 1,sNx
           fld(i,j,bi,bj) = xx_gen(i,j,bi,bj)*mask2D(i,j,bi,bj)
          ENDDO
         ENDDO
        ELSE
         DO j = 1,sNy
          DO i = 1,sNx
           fld(i,j,bi,bj) = fld(i,j,bi,bj)
     &           + xx_gen(i,j,bi,bj)*mask2D(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF ! dolog10ctrls
C     end bi/bj-loops
       ENDDO
      ENDDO

C avoid param out of [boundsVec(1) boundsVec(4)]
      CALL CTRL_BOUND_2D(fld,mask2D,xx_genarr2d_bounds(1,iarr),myThid)

      CALL EXCH_XY_RL( fld, myThid )

#ifdef ALLOW_TAPENADE
C$AD DO-NOT-DIFF
#endif
      CALL WRITE_REC_3D_RL( fnamegenOut, ctrlprec, 1,
     &                      fld, 1, optimcycle, myThid )
#ifdef ALLOW_TAPENADE
C$AD END-DO-NOT-DIFF
#endif

!write out 
#ifdef MULTISCALE_COUPLING_TAPES
        write(*,*) 'nArr for fnamegenInlo: ',nArr
        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY
        helpi=-1
        call write_glvec_rl(fnamegenInlo,' ',
     I      xx_gen_lo,nArr,helpi,myThid)

        write(*,*) 'nArr for fnamegenOutlo: ',nArr
        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY
        helpj=-1
        call write_glvec_rl(fnamegenOutlo,' ',
     I   fld_lo,nArr,helpj,myThid)
#endif /* MULTISCALE_COUPLING_TAPES */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('CTRL_MAP_GENARR2D', myThid )
#endif
#endif /* ALLOW_GENARR2D_CONTROL */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C     !ROUTINE: CTRL_MAP_GENARR3D
C     !INTERFACE:
      SUBROUTINE CTRL_MAP_GENARR3D( fld, iarr, myThid )

C     !DESCRIPTION: \bv
C     *=================================================================
C     | SUBROUTINE CTRL_MAP_GENARR3D
C     | Add the generic 3D-arrays of the
C     | control vector to the model state and update the tile halos.
C     | The control vector is defined in the header file "CTRL.h".
C     *=================================================================
C     \ev

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

#include "CTRL_SIZE.h"
#include "CTRL.h"
#include "OPTIMCYCLE.h"
#include "CTRL_GENARR.h"
#include "CTRL_DUMMY.h"
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h" /* for scaling factors */ 
#endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
      _RL     fld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER iarr
      INTEGER myThid

#ifdef ALLOW_OPENAD
      STOP 'ABNORMAL END: CTRL_MAP_GENARR3D is empty'
#else /* ALLOW_OPENAD */
#ifdef ALLOW_GENARR3D_CONTROL
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
      INTEGER bi,bj
      INTEGER i,j,k
      INTEGER numsmo,k2
      LOGICAL dowc01
      LOGICAL dosmooth
      LOGICAL doscaling
      LOGICAL dolog10ctrl
      _RL log10initval
      _RL ln10
      _RL xx_gen(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      CHARACTER*(MAX_LEN_FNAM) fnamegenIn
      CHARACTER*(MAX_LEN_FNAM) fnamegenOut
      CHARACTER*(MAX_LEN_FNAM) fnamebase
#ifdef MULTISCALE_COUPLING_TAPES
      character*(80) fnamegenInlo
      character*(80) fnamegenOutlo
#endif
      INTEGER ilgen
      LOGICAL doglobalread
      LOGICAL ladinit
      _RS mask3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER ilDir
#ifdef MULTISCALE_COUPLING_TAPES
catn: first attempt, hard-coding 365
      _RL     xx_gen_lo(1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
      _RL     fld_lo(1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
      _RL     wgenarr3d_lo(1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
      integer nArr,tmpcount,tmpind1,tmpind2,tmpind3,k,l,m,helpi,helpj,il
#endif
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('CTRL_MAP_GENARR3D', myThid )
#endif
C--   Now, read the control vector.
      doglobalread = .FALSE.
      ladinit      = .FALSE.
#ifndef ALLOW_TAPENADE
      CALL CTRL_ASSIGN(xx_gen, Nr, zeroRL, myThid)
#endif
      dosmooth=.FALSE.
      dowc01  = .FALSE.
      doscaling=.TRUE.
      dolog10ctrl=.FALSE.
      log10initval=0. _d 0
      ln10 = LOG(10. _d 0)

#ifdef MULTISCALE_COUPLING_TAPES
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO k = 1, Nr
         DO j = 1,sNy/MSfacY
          DO i = 1,sNx/MSfacY
           xx_gen_lo(i,j,k,bi,bj)=0. _d 0
           fld_lo(i,j,k,bi,bj)=0. _d 0
           wgenarr3d_lo(i,j,k,bi,bj)=0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* MULTISCALE_COUPLING_TAPES */

      numsmo = 1
      DO k2 = 1, maxCtrlProc
        IF ( xx_genarr3d_preproc(k2,iarr).EQ.'WC01' ) THEN
          dowc01 = .TRUE.
          IF (xx_genarr3d_preproc_i(k2,iarr).NE.0)
     &      numsmo = xx_genarr3d_preproc_i(k2,iarr)
        ENDIF
        IF ( (.NOT.dowc01).AND.
     &       (xx_genarr3d_preproc(k2,iarr).EQ.'smooth') ) THEN
          dosmooth = .TRUE.
          IF (xx_genarr3d_preproc_i(k2,iarr).NE.0)
     &      numsmo = xx_genarr3d_preproc_i(k2,iarr)
        ENDIF
        IF (xx_genarr3d_preproc(k2,iarr).EQ.'noscaling') THEN
          doscaling = .FALSE.
        ENDIF
        IF (xx_genarr3d_preproc_c(k2,iarr).EQ.'log10ctrl') THEN
C             fld := log10(xx_gen)
C               with initial guess for xx_gen set to log10initval
C               passed to data.ctrl, (default = 0.0)
          dolog10ctrl=.TRUE.
          log10initval = xx_genarr3d_preproc_r(k2,iarr)
        ENDIF
      ENDDO

      fnamebase = xx_genarr3d_file(iarr)
      ilgen = ILNBLNK( fnamebase )
C     Find ctrlDir (w/o trailing blanks) length
      ilDir = ilnblnk(ctrlDir)
      WRITE(fnamegenIn,'(2a,i10.10)')
     & ctrlDir(1:ilDir)//fnamebase(1:ilgen),'.',optimcycle
      WRITE(fnamegenOut,'(2a,i10.10)')
     & ctrlDir(1:ilDir)//fnamebase(1:ilgen),'.effective.',optimcycle
      il=ilnblnk( fnamegenIn )
      write(*,'(A)') 'genarr2d fnamegenIn: ', fnamegenIn(1:il)
      il=ilnblnk( fnamegenOut )
      write(*,'(A)') 'genarr2d fnamegenOut: ', fnamegenOut(1:il)
#ifdef MULTISCALE_COUPLING_TAPES
      write(fnamegenInlo(1:80),'(2a,i10.10)')
     &     fnamebase(1:ilgen),'_lo.',optimcycle
      write(fnamegenOutlo(1:80),'(2a,i10.10)')
     &     fnamebase(1:ilgen),'_lo.effective.',optimcycle
      il=ilnblnk( fnamegenInlo )
      write(*,'(A)') 'genarr2d fnamegenInlo: ', fnamegenInlo(1:il)
      il=ilnblnk( fnamegenOutlo )
      write(*,'(A)') 'genarr2d fnamegenOutlo: ', fnamegenOutlo(1:il)
#endif /* MULTISCALE_COUPLING_TAPES */

      CALL READ_REC_3D_RL( xx_genarr3d_weight(iarr), ctrlprec, Nr,
     &          wgenarr3d(1-OLx,1-OLy,1,1,1,iarr), 1, 1, myThid )

#ifdef ALLOW_AUTODIFF
      CALL ACTIVE_READ_XYZ( fnamegenIn, xx_gen, 1, doglobalread,
     & ladinit, optimcycle, myThid, xx_genarr3d_dummy(iarr) )
#else
      CALL READ_REC_XYZ_RL( fnamegenIn, xx_gen, 1, 1, myThid)
#endif

C --- Get appropriate mask for ctrl variable
      CALL CTRL_GET_MASK3D(xx_genarr3d_file(iarr), mask3D, myThid)

#ifdef MULTISCALE_COUPLING_TAPES
      write(*,*) 'genarr3d: i,j,k,xx_gen_lo(i,j,k,bi,bj),tmpcount'
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        do k=1,Nr                     !1:21
         do j=1,sNy/MSfacY            !1:16
          do i=1,sNx/MSfacX           !1:20
           tmpcount=0
           do l=1,MSfacY              !1:4
            do m=1,MSfacX             !1:4
             tmpind1 = (i-1)*MSfacX+m !1:64
             tmpind2 = (j-1)*MSfacY+l !1:80
             if(mask3D(tmpind1,tmpind2,k,bi,bj).ne. 0.) then
              tmpcount=tmpcount+1
              xx_gen_lo(i,j,k,bi,bj)=
     &        xx_gen_lo(i,j,k,bi,bj)+
     &        xx_gen(tmpind1,tmpind2,k,bi,bj)

              fld_lo(i,j,k,bi,bj)=
     &        fld_lo(i,j,k,bi,bj)+
     &        fld(tmpind1,tmpind2,k,bi,bj)
             endif
            enddo   !m
           enddo    !l
           if(tmpcount .gt. 0) then
            xx_gen_lo(i,j,k,bi,bj)=xx_gen_lo(i,j,k,bi,bj)/tmpcount
            fld_lo(i,j,k,bi,bj)=fld_lo(i,j,k,bi,bj)/tmpcount
           endif
           write(*,*) i,j,k,tmpcount,xx_gen_lo(i,j,k,bi,bj),
     &                fld_lo(i,j,k,bi,bj)
          enddo     !i
         enddo      !j
        enddo       !k
       ENDDO        !bi
      ENDDO         !bj
#endif /* MULTISCALE_COUPLING_TAPES */

C --- Do any smoothing
#ifdef ALLOW_SMOOTH
      IF (useSMOOTH) THEN
       IF (dowc01) CALL SMOOTH_CORREL3D( xx_gen, numsmo, myThid )
       IF (dosmooth) CALL SMOOTH3D( xx_gen, numsmo, myThid )
      ENDIF
#endif

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        IF (doscaling) THEN
C     scale param adjustment
         DO k = 1,Nr
          DO j = 1,sNy
           DO i = 1,sNx
            IF ( wgenarr3d(i,j,k,bi,bj,iarr).GT.0. ) THEN
             xx_gen(i,j,k,bi,bj) = xx_gen(i,j,k,bi,bj)
     &                / SQRT( wgenarr3d(i,j,k,bi,bj,iarr) )
             IF (dolog10ctrl) THEN
              xx_gen(i,j,k,bi,bj) = xx_gen(i,j,k,bi,bj) + log10initval
              xx_gen(i,j,k,bi,bj) = EXP(ln10 * xx_gen(i,j,k,bi,bj))
             ENDIF
            ELSE
             xx_gen(i,j,k,bi,bj) = 0. _d 0
            ENDIF
           ENDDO
          ENDDO
         ENDDO
C     doscaling
        ENDIF

#ifdef MULTISCALE_COUPLING_TAPES
        write(*,*) 'genarr3d: ',
     &    'i,j,k,fld_lo(i,j,k,bi,bj),tmpcount'
        DO k = 1,Nr
         DO j=1,sNy/MSfacY            !1:16
          DO i=1,sNx/MSfacX           !1:20
           IF (doscaling) then
            tmpcount=0
            DO l=1,MSfacY              !1:4
             DO m=1,MSfacX             !1:4
              tmpind1 = (i-1)*MSfacX+m   !1:64
              tmpind2 = (j-1)*MSfacY+l   !1:80
              IF ((mask3D(tmpind1,tmpind2,k,bi,bj).NE.0.).AND.
     &             (wgenarr3d(tmpind1,tmpind2,k,bi,bj,iarr).GT.0.))
     &         tmpcount=tmpcount+1
             ENDDO   !m
            ENDDO    !l
            IF(tmpcount.GT.0) THEN   !pick at least 4 pts 
             DO l=1,MSfacY                !1:4
              DO m=1,MSfacX               !1:4
               tmpind1 = (i-1)*MSfacX+m   !1:64
               tmpind2 = (j-1)*MSfacY+l   !1:80
               wgenarr3d_lo(i,j,k,bi,bj)  =
     &           max( wgenarr3d_lo(i,j,k,bi,bj),
     &           wgenarr3d(tmpind1,tmpind2,k,bi,bj,iarr) )
              ENDDO   !m
             ENDDO    !l
             fld_lo(i,j,k,bi,bj)=fld_lo(i,j,k,bi,bj)+
     &         xx_gen_lo(i,j,k,bi,bj)
     &               /sqrt(wgenarr3d_lo(i,j,k,bi,bj))
             write(*,*) i,j,k,tmpcount,fld_lo(i,j,k,bi,bj)
            ELSE
             fld_lo(i,j,k,bi,bj)= 0. _d 0
            ENDIF    !tmpcount gt 0
           ELSE
            fld_lo(i,j,k,bi,bj)=fld_lo(i,j,k,bi,bj)+
     &            xx_gen_lo(i,j,k,bi,bj)
           ENDIF     !doscaling
          ENDDO      !i
         ENDDO       !j
        ENDDO        !k
#endif /* MULTISCALE_COUPLING_TAPES */

C     add to model parameter
C     or in case of log10ctrl, fld = 10^(xx_gen)
        IF ( dolog10ctrl ) THEN
         DO k = 1,Nr
          DO j = 1,sNy
           DO i = 1,sNx
            fld(i,j,k,bi,bj) = xx_gen(i,j,k,bi,bj)*mask3D(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ELSE
         DO k = 1,Nr
          DO j = 1,sNy
           DO i = 1,sNx
            fld(i,j,k,bi,bj) = fld(i,j,k,bi,bj)
     &           + xx_gen(i,j,k,bi,bj)*mask3D(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
C     end bi/bj-loops
       ENDDO
      ENDDO

C avoid param out of [boundsVec(1) boundsVec(4)]
      CALL CTRL_BOUND_3D(fld,mask3D,xx_genarr3d_bounds(1,iarr),myThid)

C The tile exchange for xx_uvel and xx_vvel will be
C  done in CTRL_MAP_INI_GENARR.F when both
C  xx_uvel and xx_vvel are read in.
      IF (xx_genarr3d_file(iarr)(1:7).NE.'xx_uvel'.AND.
     &    xx_genarr3d_file(iarr)(1:7).NE.'xx_vvel')
     &    CALL EXCH_XYZ_RL( fld, myThid )

#ifdef ALLOW_TAPENADE
C$AD DO-NOT-DIFF
#endif
      CALL WRITE_REC_3D_RL( fnamegenOut, ctrlprec, Nr,
     &                      fld, 1, optimcycle, myThid )
#ifdef ALLOW_TAPENADE
C$AD END-DO-NOT-DIFF
#endif

!write out 
#ifdef MULTISCALE_COUPLING_TAPES
        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*Nr
        write(*,*) 'genarr3d, nArr for fnamegenInlo: ',nArr
        helpi=-1
        call write_glvec_rl(fnamegenInlo,' ',
     I      xx_gen_lo,nArr,helpi,myThid)

        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*Nr
        write(*,*) 'genarr3d, nArr for fnamegenOutlo: ',nArr
        helpj=-1
        call write_glvec_rl(fnamegenOutlo,' ',
     I   fld_lo,nArr,helpj,myThid)
#endif /* ALLOW_GENARR2D_CONTROL */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('CTRL_MAP_GENARR3D', myThid )
#endif
#endif /* ALLOW_GENARR3D_CONTROL */
#endif /* ALLOW_OPENAD */

      RETURN
      END
