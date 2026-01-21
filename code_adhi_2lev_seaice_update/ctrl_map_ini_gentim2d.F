#include "CTRL_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: CTRL_MAP_INI_GENTIM2D
C     !INTERFACE:
      SUBROUTINE CTRL_MAP_INI_GENTIM2D( myThid )

C     !DESCRIPTION: \bv
C     *================================================================
C     | SUBROUTINE CTRL_MAP_INI_GENTIM2D
C     | Dimensionalize and preprocess time variable controls.
C     *================================================================
C     \ev

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "CTRL_SIZE.h"
#include "CTRL.h"
#include "OPTIMCYCLE.h"
#include "CTRL_DUMMY.h"
#include "CTRL_GENARR.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_FIELDS.h"
#endif
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h" /* for scaling factors */
#endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == routine arguments ==
      INTEGER myThid

#ifdef ALLOW_GENTIM2D_CONTROL
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
C     == local variables ==
      integer iarr
      integer numsmo
C     Tapenade requires these 4 variables to be declared as arrays of
C     size maxCtrlTim2D. If declared as scalars, it does not store their
C     temporary values in the do loop during reverse pass.
      character*(MAX_LEN_FNAM) fnamegenIn(1:maxCtrlTim2D)
      character*(MAX_LEN_FNAM) fnamegenOut(1:maxCtrlTim2D)
      character*(MAX_LEN_FNAM) fnamegenTmp(1:maxCtrlTim2D)
      character*(MAX_LEN_FNAM) fnamebase(1:maxCtrlTim2D)
#ifdef MULTISCALE_COUPLING_TAPES
      character*(80) fnamegenInlo
      character*(80) fnamegenOutlo
#endif
      integer startrec
      integer endrec
      integer diffrec
      integer iRec, jrec, kRec, lRec
      integer replicated_nrec
      integer replicated_ntimes
      logical doglobalread
      logical ladinit
      logical dowc01
      logical dosmooth
      logical doscaling
      _RL     xx_gen(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS     mask2D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef MULTISCALE_COUPLING_TAPES
catn: first attempt, hard-coding 365
      _RL     xx_gen_lo(1:sNx/MSfacX,1:sNy/MSfacY,365,nSx,nSy)
      _RL     xx_gen_scale_lo(1:sNx/MSfacX,1:sNy/MSfacY,365,nSx,nSy)
      _RL     wgentim2d_lo(1:sNx/MSfacX,1:sNy/MSfacY,365,nSx,nSy)
      integer nArr,tmpcount,tmpind1,tmpind2,tmpind3,k,l,m,helpi,helpj,il

#endif
#ifdef ALLOW_ECCO
      _RL     xx_gen_tmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer nyearsINT
      _RL     recip_nyearsRL
#endif
      integer bi,bj
      integer i,j,k2
      integer ilgen
      integer ilDir
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('CTRL_MAP_INI_GENTIM2D',myThid)
#endif /* ALLOW_DEBUG */

      write(*,*) 'entering ctrl_map_ini_gentim2d'
C--   Now, read the control vector.
      doglobalread = .false.
      ladinit      = .false.

C     Find ctrlDir (w/o trailing blanks) length
      ilDir = ilnblnk(ctrlDir)

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
          xx_gen(i,j,bi,bj)=0. _d 0
#ifdef ALLOW_ECCO
          xx_gen_tmp(i,j,bi,bj)=0. _d 0
#endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#ifdef MULTISCALE_COUPLING_TAPES
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO k = 1, 365
         DO j = 1,sNy/MSfacY
          DO i = 1,sNx/MSfacY
           xx_gen_lo(i,j,k,bi,bj)=0. _d 0
           xx_gen_scale_lo(i,j,k,bi,bj)=0. _d 0
           wgentim2d_lo(i,j,k,bi,bj)=0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* MULTISCALE_COUPLING_TAPES */
                                             !}}} [init]
      !write(*,*) 'maxCtrlTim2D: ',maxCtrlTim2D
C--   generic 2D control variables
      DO iarr = 1, maxCtrlTim2D                   !{{{ [iarr 1-maxCtrlTim2D]

       diffrec=0
       startrec=0
       endrec=0

#ifndef ALLOW_OPENAD
       if (xx_gentim2d_weight(iarr).NE.' ') then !{{{ [xx_gentim2d_weight ne '']
#endif

        ilgen=ilnblnk( xx_gentim2d_file(iarr) )
        fnamebase(iarr) = xx_gentim2d_file(iarr)(1:ilgen)
        !il=ilnblnk( fnamebase )
        !write(*,'(A)') 'fnamebase(1:il): ', fnamebase(1:il)
        call ctrl_init_rec ( fnamebase(iarr),
     I       xx_gentim2d_startdate1(iarr),
     I       xx_gentim2d_startdate2(iarr),
     I       xx_gentim2d_period(iarr),
     I       1,
     O       xx_gentim2d_startdate(1,iarr),
     O       diffrec, startrec, endrec,
     I       myThid )

C From CTRL_INIT_REC one gets [start,end,diff]rec.
C Below, file xx_*.$iter.data is read in from records startrec to endrec,
C followed by a write to file xx_*.effective.$iter.data of record 1 to diffrec
C (see TAF-generated S/R ctrl_map_ini_gentim2dmd).

#ifdef ALLOW_CTRL_DEBUG
C Note: In the "md" (i.e. TAF modified) version of this fwd S/R, DEBUG calls and
C print statements are removed. See "Automatic Differentiation" chap. in the doc
        write(*,*) 'AA: iarr,xx_gentim2d_startdate(1,iarr): ',
     &              iarr,xx_gentim2d_startdate(1,iarr)
#endif
        dosmooth=.false.
        dowc01  = .false.
        doscaling=.true.

        numsmo=1
        do k2 = 1, maxCtrlProc           !{{{ [do[wc01,smooth,scaling]]
         if (xx_gentim2d_preproc(k2,iarr).EQ.'WC01') then
          dowc01=.TRUE.
          if (xx_gentim2d_preproc_i(k2,iarr).NE.0)
     &         numsmo=xx_gentim2d_preproc_i(k2,iarr)
         endif
         if ((.NOT.dowc01).AND.
     &        (xx_gentim2d_preproc(k2,iarr).EQ.'smooth')) then
          dosmooth=.TRUE.
          if (xx_gentim2d_preproc_i(k2,iarr).NE.0)
     &         numsmo=xx_gentim2d_preproc_i(k2,iarr)
         endif
         if (xx_gentim2d_preproc(k2,iarr).EQ.'noscaling') then
          doscaling=.FALSE.
         endif
        enddo                   !}}} [do[wc01,smooth,scaling]]

        ilgen=ilnblnk( xx_gentim2d_file(iarr) )  !{{{ [fname[base,In,Out]] 
        fnamebase(iarr) = xx_gentim2d_file(iarr)(1:ilgen)
        write(fnamegenIn(iarr),'(2a,i10.10)')
     &   ctrlDir(1:ilDir)//fnamebase(iarr)(1:ilgen),'.',optimcycle
        write(fnamegenOut(iarr),'(2a,i10.10)')
     &   ctrlDir(1:ilDir)//fnamebase(iarr)(1:ilgen),
     &   '.effective.',optimcycle
        write(fnamegenTmp(iarr),'(2a,i10.10)')
     &   ctrlDir(1:ilDir)//fnamebase(iarr)(1:ilgen),'.tmp.',optimcycle
        il=ilnblnk( fnamegenIn )
        write(*,'(A)') 'gentim2d fnamegenIn: ', fnamegenIn(1:il)
        il=ilnblnk( fnamegenOut )
        write(*,'(A)') 'gentim2d fnamegenOut: ', fnamegenOut(1:il)
#ifdef MULTISCALE_COUPLING_TAPES
        write(fnamegenInlo(1:80),'(2a,i10.10)')
     &       fnamebase(1:ilgen),'_lo.',optimcycle
        write(fnamegenOutlo(1:80),'(2a,i10.10)')
     &       fnamebase(1:ilgen),'_lo.effective.',optimcycle
        il=ilnblnk( fnamegenInlo )
        write(*,'(A)') 'gentim2d fnamegenInlo: ', fnamegenInlo(1:il)
        il=ilnblnk( fnamegenOutlo )
        write(*,'(A)') 'gentim2d fnamegenOutlo: ', fnamegenOutlo(1:il)
#endif /* MULTISCALE_COUPLING_TAPES */
                                                !}}} [fname[base,In,Out]]

C--   docycle
!{{{ [docycle] 
        replicated_nrec=diffrec
        replicated_ntimes=0
        !write(*,*) 'ctrl_map_ini_gentim2d A: k2,iarr,',
c     &  'replicated_[nrec,ntimes],diffrec,xx_gentim2d_preproc[,i]: '
        do k2 = 1, maxCtrlProc
         if (xx_gentim2d_preproc(k2,iarr).EQ.'docycle') then
          if (xx_gentim2d_preproc_i(k2,iarr).NE.0) then
           replicated_nrec=min(diffrec,xx_gentim2d_preproc_i(k2,iarr))
           replicated_ntimes=
     &          int(float(diffrec)/float(replicated_nrec))
           if (replicated_ntimes*replicated_nrec.LT.diffrec)
     &          replicated_ntimes=replicated_ntimes+1
           if (replicated_ntimes*replicated_nrec.GT.diffrec)
     &          replicated_ntimes=replicated_ntimes-1
          endif
         endif
c         write(*,*) k2,iarr,replicated_nrec,replicated_ntimes,diffrec,
c     &   xx_gentim2d_preproc(k2,iarr),xx_gentim2d_preproc_i(k2,iarr)
        enddo

#ifdef ALLOW_CTRL_DEBUG
        write(*,'(A,1x,5I6)')
     &   'AD:[start,end,diff]rec, replicated_[nrec,ntimes]: ',
     &              startrec, endrec, diffrec,
     &              replicated_nrec,replicated_ntimes
#endif
        DO jrec = 1, replicated_ntimes+1
         DO iRec = 1, replicated_nrec
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE xx_gentim2d_dummy = dummytape, key = 1 , kind = isbyte
#endif
          kRec=replicated_nrec*(jrec-1)+iRec
c          write(*,*) 'ctrl_map_ini_gentim2d B: [i,j,k,end]rec:',
c     &      irec,jrec,krec,endrec
          lRec=startrec+iRec-1
          IF (kRec.LE.endrec) THEN
#ifdef ALLOW_CTRL_DEBUG
           ilgen=ilnblnk( fnamegenIn(iarr) )
           write(*,*) 'AE: iarr,[i,j,k,l]rec: ',iarr,iRec,jrec,lRec
           write(*,*) 'AF: fnamegenIn: ', fnamegenIn(iarr)(1:ilgen)
#endif
catn: note that we read xx_gen and calc xx_gen_lo here only if docycle is true

C fnamegenIn is xx_*.$iter.data, so the required access records starts at
C startrec+iRec-1 instead of iRec in the read call below. In adjoint mode, S/R
C ADACTIVE_READ_XY reads adxx_*.$iter.data in reverse order, with required
C access from endrec back to startrec. This requires creating adxx_*.$iter.data
C file of size endrec, with valid gradients from records startrec to endrec
C but all zeros from 1 to startrec-1. See documentation:
C https://mitgcm.readthedocs.io/en/latest/ocean_state_est/ocean_state_est.html#ctrl-model-parameter-adjustment-capability
#ifdef ALLOW_AUTODIFF
           CALL ACTIVE_READ_XY( fnamegenIn(iarr), xx_gen, lRec,
     &          doglobalread, ladinit, optimcycle,           !.false.,.false.
     &          myThid, xx_gentim2d_dummy(iarr) )
#else /* ALLOW_AUTODIFF */ 
           CALL READ_REC_XY_RL( fnamegenIn(iarr), xx_gen, lRec,
     &                          1, myThid )
#endif /* ALLOW_AUTODIFF */ 
C Here, we write record kRec (from 1 to diffrec) to file xx_*.effective.$iter
#ifdef ALLOW_AUTODIFF
#ifdef ALLOW_CTRL_DEBUG
           ilgen=ilnblnk( fnamegenOut(iarr) )
           write(*,*) 'AG: iarr,[i,j,k]rec: ',iarr,iRec,jrec,kRec
           write(*,*) 'AH: fnamegenOut: ', fnamegenOut(iarr)(1:ilgen)
#endif
           CALL ACTIVE_WRITE_XY( fnamegenOut(iarr), xx_gen, kRec,
     &          optimcycle, myThid, xx_gentim2d_dummy(iarr) )
#else
           CALL WRITE_REC_XY_RL( fnamegenOut(iarr), xx_gen, kRec,
     &          1 , myThid )
#endif
          ENDIF !krec<endrec
         ENDDO  !irec
        ENDDO   !jrec
                                            !}}} [docycle]

C--   rmcycle
!{{{ [rmcycle]
#ifdef ALLOW_ECCO
        replicated_nrec=diffrec
        replicated_ntimes=0
        do k2 = 1, maxCtrlProc
         if (xx_gentim2d_preproc(k2,iarr).EQ.'rmcycle') then
          if (xx_gentim2d_preproc_i(k2,iarr).NE.0) then
           replicated_nrec=min(diffrec,xx_gentim2d_preproc_i(k2,iarr))
           replicated_ntimes=
     &          int(float(diffrec)/float(replicated_nrec))
           if (replicated_ntimes*replicated_nrec.LT.diffrec)
     &          replicated_ntimes=replicated_ntimes+1
           if (replicated_ntimes*replicated_nrec.GT.diffrec)
     &          replicated_ntimes=replicated_ntimes-1
          endif
         endif
        enddo

#ifdef ALLOW_CTRL_DEBUG
        write(*,'(A,1x,5I6)')
     &  'AI:[start,end,diff]rec, replicated_[nrec,ntimes]: ',
     &      startrec,endrec,diffrec,replicated_nrec,replicated_ntimes
#endif

        IF (replicated_ntimes.GT.0) THEN

C     create cyclic average

         nyearsINT=1+int((diffrec-replicated_nrec)/replicated_nrec)
         recip_nyearsRL=1. _d 0/float(nyearsINT)

#ifdef ALLOW_CTRL_DEBUG
         write(*,*) 'AJ: nyears[INT,RL]: ',nyearsINT,recip_nyearsRL
#endif

         DO iRec = 1, replicated_nrec !{{{ [irec 1-replicated_nrec]

          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j = 1-OLy,sNy+OLy
             DO i = 1-OLx,sNx+OLx
              xx_gen(i,j,bi,bj) = zeroRL
             ENDDO
            ENDDO
           ENDDO
          ENDDO

          DO jrec=1,nyearsINT          !{{{ [jrec 1-nyearsINT]
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE xx_gentim2d_dummy = dummytape, key = 1 , kind = isbyte
#endif
           kRec=iRec+(jrec-1)*replicated_nrec
#ifdef ALLOW_CTRL_DEBUG
           write(*,*) 'AK: iarr,[i,j,k]rec: ',iarr,iRec,jrec,kRec
           ilgen=ilnblnk( fnamegenOut(iarr) )
           write(*,*) 'AL: fnamegenOut: ',fnamegenOut(iarr)(1:ilgen)
#endif
#ifdef ALLOW_AUTODIFF
           CALL ACTIVE_READ_XY( fnamegenOut(iarr), xx_gen_tmp, kRec,
     &          doglobalread, ladinit, optimcycle,
     &          myThid, xx_gentim2d_dummy(iarr) )
#else
           CALL READ_REC_XY_RL( fnamegenOut(iarr), xx_gen_tmp, kRec,
     &          1, myThid )
#endif
           DO bj=myByLo(myThid),myByHi(myThid)
            DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1,sNy
              DO i = 1,sNx
               xx_gen(i,j,bi,bj) = xx_gen(i,j,bi,bj)
     &              +xx_gen_tmp(i,j,bi,bj)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
C     end jrec
          ENDDO                       !}}} [jrec 1-nyearsINT]

          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j = 1,sNy
             DO i = 1,sNx
              xx_gen(i,j,bi,bj) = xx_gen(i,j,bi,bj) * recip_nyearsRL
             ENDDO
            ENDDO
           ENDDO
          ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE xx_gentim2d_dummy = dummytape, key = 1 , kind = isbyte
#endif

#ifdef ALLOW_AUTODIFF
          CALL ACTIVE_WRITE_XY( fnamegenTmp(iarr), xx_gen, iRec,
     &         optimcycle, myThid, xx_gentim2d_dummy(iarr) )
#else
          CALL WRITE_REC_XY_RL( fnamegenTmp(iarr), xx_gen, iRec,
     &         1, myThid )
#endif

         ENDDO                      !}}} [irec 1-replicated_nrec]

C     subtract cyclic average
         DO jrec = 1, replicated_ntimes+1 !{{{ [[i,j,k]rec sub cyclic] 
          DO iRec = 1, replicated_nrec
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE xx_gentim2d_dummy = dummytape, key = 1 , kind = isbyte
#endif
           kRec=replicated_nrec*(jrec-1)+iRec
#ifdef ALLOW_CTRL_DEBUG
           write(*,*) 'AM: iarr,[i,j,k]rec: ',iarr,iRec,jrec,kRec
           ilgen=ilnblnk( fnamegenOut(iarr) )
           write(*,*) 'AN: fnamegenOut: ',fnamegenOut(iarr)(1:ilgen)
#endif
           IF (kRec.LE.diffrec) THEN
#ifdef ALLOW_AUTODIFF
            CALL ACTIVE_READ_XY( fnamegenOut(iarr), xx_gen, kRec,
     &           doglobalread, ladinit, optimcycle,
     &           myThid, xx_gentim2d_dummy(iarr) )
#else
            CALL READ_REC_XY_RL( fnamegenOut(iarr), xx_gen, kRec,
     &           1, myThid )
#endif
#ifdef ALLOW_AUTODIFF
            CALL ACTIVE_READ_XY( fnamegenTmp(iarr), xx_gen_tmp, iRec,
     &           doglobalread, ladinit, optimcycle,
     &           myThid, xx_gentim2d_dummy(iarr) )
#else
            CALL READ_REC_XY_RL( fnamegenTmp(iarr), xx_gen_tmp, iRec,
     &           1, myThid )
#endif
            DO bj=myByLo(myThid),myByHi(myThid)
             DO bi=myBxLo(myThid),myBxHi(myThid)
              DO j = 1,sNy
               DO i = 1,sNx
                xx_gen(i,j,bi,bj)=xx_gen(i,j,bi,bj)
     &              -xx_gen_tmp(i,j,bi,bj)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
#ifdef ALLOW_AUTODIFF
            CALL ACTIVE_WRITE_XY( fnamegenOut(iarr), xx_gen, kRec,
     &           optimcycle, myThid, xx_gentim2d_dummy(iarr) )
#else
            CALL WRITE_REC_XY_RL( fnamegenOut(iarr), xx_gen, kRec,
     &           1, myThid )
#endif
           ENDIF
          ENDDO
         ENDDO                            !}}} [[i,j,k]rec sub cyclic] 

        ENDIF
#endif /* ifdef ALLOW_ECCO */
!}}} [rmcycle]

C--   scaling and smoothing

C The access of records startrec to endrec in  xx_*.*iter was already
C done correctly above. From here, we read in xx_*.effective.$iter.data
C of size diffrec, so no more fix of record is needed from here on out.

catn: if we did not enter loop related to docycle above, then we do not
catn  yet have xx_gen or xx_gen_lo.  So below, xx_gen is being read in
catn  either from raw file or from processes from above.  Thus, we now 
catn  need to figure out how to get xx_gen_lo if we didnt set it above.
catn  I think it is actually the best approach to get xx_gen_lo here
catn  instead of way above.

        DO iRec = 1, diffrec
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE xx_gentim2d_dummy = dummytape, key = 1 , kind = isbyte
#endif

#ifdef ALLOW_CTRL_DEBUG
         ilgen=ilnblnk( fnamegenOut(iarr) )
         lRec=startrec+iRec-1
         write(*,'(A,1x,A,1x,2I6)') 'AO: fnamegenOut, iRec,lRec ',
     &       fnamegenOut(iarr)(1:ilgen),iRec,lRec
#endif
#ifdef ALLOW_AUTODIFF
         CALL ACTIVE_READ_XY( fnamegenOut(iarr), xx_gen, iRec,
     &        doglobalread, ladinit, optimcycle,
     &        myThid, xx_gentim2d_dummy(iarr) )
#else
         CALL READ_REC_XY_RL( fnamegenOut(iarr), xx_gen, iRec,
     &        1, myThid )
#endif

#ifndef ALLOW_OPENAD
         jrec=1
C Not clear which record to use for timevariable weights ;
C For now, just assumes records from startrec to endrec are available in file
         do k2 = 1, maxCtrlProc
          if (xx_gentim2d_preproc(k2,iarr).EQ.'variaweight') then
            jrec=startrec+iRec-1
          endif
         enddo
         CALL READ_REC_3D_RL( xx_gentim2d_weight(iarr), ctrlprec, 1,
     &             wgentim2d(1-OLx,1-OLy,1,1,iarr), jrec, 1, myThid )

C--   Get appropriate mask
         call ctrl_get_mask2D(xx_gentim2d_file(iarr), mask2D, myThid)

#ifdef MULTISCALE_COUPLING_TAPES
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           tmpind3 = krec               !
           do j=1,sNy/MSfacY            !1:16
            do i=1,sNx/MSfacX           !1:20
             tmpcount=0
             do l=1,MSfacY              !1:4
              do m=1,MSfacX             !1:4
               tmpind1 = (i-1)*MSfacX+m !1:64
               tmpind2 = (j-1)*MSfacY+l !1:80
               if(mask2D(tmpind1,tmpind2,bi,bj).ne.0.) then
                tmpcount=tmpcount+1
                xx_gen_lo(i,j,tmpind3,bi,bj)=
     &          xx_gen_lo(i,j,tmpind3,bi,bj)+
     &          xx_gen(tmpind1,tmpind2,bi,bj)
               endif
              enddo   !m
             enddo    !l
             if(tmpcount .gt. 0)
     &        xx_gen_lo(i,j,tmpind3,bi,bj)=
     &        xx_gen_lo(i,j,tmpind3,bi,bj)/tmpcount
            enddo     !i
           enddo      !j
          ENDDO        !bi
         ENDDO         !bj
#endif /* MULTISCALE_COUPLING_TAPES */

#ifdef ALLOW_SMOOTH
         IF (useSMOOTH) THEN
          IF (dowc01) call smooth_correl2d(xx_gen,mask2D,numsmo,myThid)
          IF (dosmooth) call smooth2d(xx_gen,mask2D,numsmo,myThid)
         ENDIF
catn: we do not use xx_gen after smooth for lowres , because
catn: the idea of lowres is already mimicking smoothing.  Thus
catn: we do not avg smooth but approx avg raw = smooth.
#endif /* ALLOW_SMOOTH */

         DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
           DO j = 1,sNy
            DO i = 1,sNx
             if ((mask2D(i,j,bi,bj).NE.0.).AND.
     &            (wgentim2d(i,j,bi,bj,iarr).GT.0.)) then
              IF (doscaling) then
               xx_gen(i,j,bi,bj)=xx_gen(i,j,bi,bj)
     &              /sqrt(wgentim2d(i,j,bi,bj,iarr))
              ENDIF             ! IF (doscaling) then
             else
              xx_gen(i,j,bi,bj)=0. _d 0
             endif
            ENDDO
           ENDDO
          ENDDO
         ENDDO

#ifdef MULTISCALE_COUPLING_TAPES
         DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
           DO j=1,sNy/MSfacY            !1:16
            DO i=1,sNx/MSfacX           !1:20
             IF (doscaling) then
              tmpcount=0
              DO l=1,MSfacY              !1:4
               DO m=1,MSfacX             !1:4
                IF ((mask2D(i,j,bi,bj).NE.0.).AND.
     &               (wgentim2d(i,j,bi,bj,iarr).GT.0.))
     &           tmpcount=tmpcount+1
               ENDDO   !m
              ENDDO    !l
              IF(tmpcount.GT.MSfacX) THEN   !pick at least 4 pts
               DO l=1,MSfacY                !1:4
                DO m=1,MSfacX               !1:4
                 tmpind1 = (i-1)*MSfacX+m   !1:64
                 tmpind2 = (j-1)*MSfacY+l   !1:80
                 wgentim2d_lo(i,j,tmpind3,bi,bj)  =
     &             max( wgentim2d_lo(i,j,tmpind3,bi,bj),
     &             wgentim2d(tmpind1,tmpind2,bi,bj,iarr) )
                ENDDO   !m
               ENDDO    !l
               xx_gen_scale_lo(i,j,tmpind3,bi,bj)=
     &           xx_gen_lo(i,j,tmpind3,bi,bj)
     &                 /sqrt(wgentim2d_lo(i,j,tmpind3,bi,bj))
              ELSE
               xx_gen_scale_lo(i,j,tmpind3,bi,bj)= 0. _d 0
              ENDIF    !tmpcount gt 0
             ELSE
              xx_gen_scale_lo(i,j,tmpind3,bi,bj)=
     &              xx_gen_lo(i,j,tmpind3,bi,bj)
             ENDIF     !doscaling
            ENDDO      !i
           ENDDO       !j
          ENDDO        !bi
         ENDDO         !bj
#endif /* MULTISCALE_COUPLING_TAPES */
#endif /* ALLOW_OPENAD */

catn: for now, i am skipping setting bounds...
         CALL CTRL_BOUND_2D(xx_gen,mask2D,
     &        xx_gentim2d_bounds(1,iarr),myThid)

         CALL EXCH_XY_RL ( xx_gen , myThid )

#ifdef ALLOW_CTRL_DEBUG
         ilgen=ilnblnk( fnamegenOut(iarr) )
         write(*,'(A,1x,I6,1x,A)')
     &      'AQ: iRec,fnamegenOut: ',iRec,fnamegenOut(iarr)(1:ilgen)
#endif
#ifdef ALLOW_AUTODIFF
         CALL ACTIVE_WRITE_XY( fnamegenOut(iarr), xx_gen, iRec,
     &        optimcycle, myThid, xx_gentim2d_dummy(iarr) )
#else
         CALL WRITE_REC_XY_RL( fnamegenOut(iarr), xx_gen, iRec,
     &        1, myThid )
#endif
C--   end iRec loop
        ENDDO

!write out after finish nrec loop:
#ifdef MULTISCALE_COUPLING_TAPES
        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*365
        helpi=-1
        call write_glvec_rl(fnamegenInlo,' ',
     I      xx_gen_lo,nArr,helpi,myThid)

        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*365
        helpj=-1
        call write_glvec_rl(fnamegenOutlo,' ',
     I   xx_gen_scale_lo,nArr,helpj,myThid)
#endif /* MULTISCALE_COUPLING_TAPES */

#ifndef ALLOW_OPENAD
       endif                            !}}} [xx_gentim2d_weight ne '']
#endif

C--   end iarr loop
      ENDDO                             !}}} [iarr 1-maxCtrlTim2D]

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('CTRL_MAP_INI_GENTIM2D',myThid)
#endif /* ALLOW_DEBUG */
#endif /* ALLOW_GENTIM2D_CONTROL */

      RETURN
      END
