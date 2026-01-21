#include "ECCO_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h" /* for multiscale use */
#endif

C--  File cost_generic.F: Generic cost function routines
C--   Contents
C--   o cost_generic
C--   o cost_genloop

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: cost_generic
C     !INTERFACE:
      subroutine cost_generic(
     &     nnzbar, localbarfile, dummy,
     &     nnzobs, localobsfile, localerrfile,
     &     mult_local, nrecloc, nrecobs,
     &     localstartdate, localperiod,
     &     ylocmask, spminloc, spmaxloc, spzeroloc,
     &     preproc, preproc_c, preproc_i, preproc_r,
     &     posproc, posproc_c, posproc_i, posproc_r,
     &     outlev, outname,
     &     objf_local, num_local,
     &     myIter, myTime, myThid )

C     !DESCRIPTION: \bv
C     Generic routine for evaluating time-dependent
c       cost function contribution
C     \ev

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_ECCO
# include "ECCO_SIZE.h"
# include "ECCO.h"
#endif
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h" /* for scaling factors */
#endif
#endif


C     !INPUT/OUTPUT PARAMETERS:
      integer myIter
      integer myThid
      integer nnzbar, nnzobs
      integer nrecloc, nrecobs
      integer localstartdate(4)
      integer outlev
      integer preproc_i(NGENPPROC)
      integer posproc_i(NGENPPROC)

      _RL objf_local(nSx,nSy)
      _RL num_local(nSx,nSy)
#ifdef ALLOW_AUTODIFF 
#ifdef MULTISCALE_COUPLING_TAPES
      _RL objf_local_lo(nSx,nSy)
      _RL num_local_lo(nSx,nSy)
#endif
#endif
      _RL dummy
      _RL mult_local
      _RL myTime
      _RL localperiod
      _RL spminloc
      _RL spmaxloc
      _RL spzeroloc
      _RL preproc_r(NGENPPROC)
      _RL posproc_r(NGENPPROC)

      character*(1) ylocmask
      character*(MAX_LEN_FNAM) localbarfile
      character*(MAX_LEN_FNAM) localobsfile
      character*(MAX_LEN_FNAM) localerrfile
      character*(MAX_LEN_FNAM) preproc(NGENPPROC)
      character*(MAX_LEN_FNAM) preproc_c(NGENPPROC)
      character*(MAX_LEN_FNAM) posproc(NGENPPROC)
      character*(MAX_LEN_FNAM) posproc_c(NGENPPROC)
      character*(MAX_LEN_FNAM) outname

#ifdef ALLOW_ECCO
C     !LOCAL VARIABLES:
      integer bi,bj,k2
      logical domean, doanom

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
      integer nxb, nyb
#endif
#endif
      _RL localdifmean1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localdifmean2  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
CEOP

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
      nyb=sNy/MSfacY
      nxb=sNx/MSfacX
#endif
#endif

c--   Initialise local variables.

      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          objf_local(bi,bj) = 0. _d 0
          num_local(bi,bj) = 0. _d 0
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
          objf_local_lo(bi,bj) = 0. _d 0 
          num_local_lo(bi,bj) = 0. _d 0
#endif
#endif
        ENDDO
      ENDDO

      call ecco_zero(localdifmean1,Nr,zeroRL,myThid)
      call ecco_zero(localdifmean2,Nr,zeroRL,myThid)

      domean=.FALSE.
      doanom=.FALSE.
      do k2 = 1, NGENPPROC
          if (preproc(k2).EQ.'mean') domean=.TRUE.
          if (preproc(k2).EQ.'anom') doanom=.TRUE.
      enddo
      write(*,*) 'domean,doanom: ',domean,doanom !F,F

C Extra time loop to compute time-mean fields and costs
      if ( (.NOT. ( localobsfile.EQ.' ' ) )         !{{{ [not enter]
     &   .AND. ( domean .OR. doanom ) ) then
        call cost_genloop(
     &     localdifmean1,localdifmean2,.FALSE.,
     &     nnzbar, localbarfile, dummy,
     &     nnzobs, localobsfile, localerrfile,
     &     mult_local, nrecloc, nrecobs,
     &     localstartdate, localperiod,
     &     ylocmask, spminloc, spmaxloc, spzeroloc,
     &     preproc, preproc_c, preproc_i, preproc_r,
     &     posproc, posproc_c, posproc_i, posproc_r,
     &     outlev, outname,
     &     objf_local, num_local,
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
     &     objf_local_lo, num_local_lo,
#endif
#endif
     &     myIter, myTime, myThid )
      endif                                        !}}}

      call ecco_zero(localdifmean1,Nr,zeroRL,myThid)

      if ((.NOT.(localobsfile.EQ.' ')).AND.(.NOT.domean)) then
        call cost_genloop(
     &     localdifmean2,localdifmean1,.TRUE.,
     &     nnzbar, localbarfile, dummy,
     &     nnzobs, localobsfile, localerrfile,
     &     mult_local, nrecloc, nrecobs,
     &     localstartdate, localperiod,
     &     ylocmask, spminloc, spmaxloc, spzeroloc,
     &     preproc, preproc_c, preproc_i, preproc_r,
     &     posproc, posproc_c, posproc_i, posproc_r,
     &     outlev, outname,
     &     objf_local, num_local,
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
     &     objf_local_lo, num_local_lo,
#endif
#endif
     &     myIter, myTime, myThid )
      endif

#endif /* ALLOW_ECCO */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: cost_genloop
C     !INTERFACE:
      subroutine cost_genloop(
     &     localdifmeanIn,localdifmeanOut, addVariaCost,
     &     nnzbar, localbarfile, dummy,
     &     nnzobs, localobsfile, localerrfile,
     &     mult_local, nrecloc, nrecobs,
     &     localstartdate, localperiod,
     &     ylocmask, spminloc, spmaxloc, spzeroloc,
     &     preproc, preproc_c, preproc_i, preproc_r,
     &     posproc, posproc_c, posproc_i, posproc_r,
     &     outlev, outname,
     &     objf_local, num_local,
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
     &     objf_local_lo, num_local_lo,
#endif
#endif
     &     myIter, myTime, myThid )

C     !DESCRIPTION: \bv
C     Generic routine for evaluating time-dependent
c       cost function contribution
C     \ev

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_ECCO
# include "ECCO_SIZE.h"
# include "ECCO.h"
#endif
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h"        /* for tape level sizes */
#endif
#endif


C     !INPUT/OUTPUT PARAMETERS:
      integer myIter
      integer myThid
      integer nnzbar, nnzobs
      integer nrecloc, nrecobs
      integer localstartdate(4)
      integer outlev
      integer preproc_i(NGENPPROC)
      integer posproc_i(NGENPPROC)

      _RL objf_local(nSx,nSy)
      _RL num_local(nSx,nSy)
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
      _RL objf_local_lo(nSx,nSy)
      _RL num_local_lo(nSx,nSy)
#endif
#endif
      _RL dummy
      _RL mult_local
      _RL myTime
      _RL localperiod
      _RL spminloc
      _RL spmaxloc
      _RL spzeroloc
      _RL preproc_r(NGENPPROC)
      _RL posproc_r(NGENPPROC)

      character*(1) ylocmask
      character*(MAX_LEN_FNAM) localbarfile
      character*(MAX_LEN_FNAM) localobsfile
      character*(MAX_LEN_FNAM) localerrfile
      character*(MAX_LEN_FNAM) preproc(NGENPPROC)
      character*(MAX_LEN_FNAM) preproc_c(NGENPPROC)
      character*(MAX_LEN_FNAM) posproc(NGENPPROC)
      character*(MAX_LEN_FNAM) posproc_c(NGENPPROC)
      character*(MAX_LEN_FNAM) outname

      logical addVariaCost
      _RL localdifmeanIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localdifmeanOut(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

#ifdef ALLOW_ECCO
C     !FUNCTIONS:
      integer  ilnblnk
      external ilnblnk

C     !LOCAL VARIABLES:
      integer bi,bj
      integer i,j,k
      integer irec, jrec
      integer il, k2
      integer localrec, obsrec
      integer nrecloop, nrecclim, k2smooth
      logical domean, doanom, dovarwei, doclim, dosmooth, dosumsq
      logical dooffset, domaskmindepth

      _RL localmask  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL localbar   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localweight(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localtmp   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localobs   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localdif   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL difmask    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL localdifmsk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localdifsum(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL localdifnum(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
catn: fold the nzobs into 3rd dim, for now hardcode to only
catn: allow 365 records, for testing 
catn: The way the barfields are aggregated and written out
catn  is not easy to downscale.  So i think for now the easiest
catn  route may be to simply write it out right here! (though
catn  inappropriate, but actually most sufficient i think.
c     == local variables ==
      _RL localobslo(1:sNx/MSfacX,1:sNy/MSfacY,365*Nr,nSx,nSy)
      _RL localdiflo(1:sNx/MSfacX,1:sNy/MSfacY,365*Nr,nSx,nSy)
      _RL localbarlo(1:sNx/MSfacX,1:sNy/MSfacY,365*Nr,nSx,nSy)
      _RL localweightlo(1:sNx/MSfacX,1:sNy/MSfacY,365*Nr,nSx,nSy)
      _RL difmasklo(1:sNx/MSfacX,1:sNy/MSfacY,365*Nr,nSx,nSy)
      _RL localmasklo(1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
      _RL localdifmeanOutlo(1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
c     _RL localdifmsklo    (1:sNx/MSfacX,1:sNy/MSfacY,Nr,nSx,nSy)
      integer l,m,tmpind1,tmpind2,tmpind3,nArr, krec, tmpcount
c      integer nxb,nyb
c     integer ioUnit
#endif
#endif /* ALLOW_AUTODIFF */

#ifdef ALLOW_SMOOTH
      _RL tmpfld2D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS mask2D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif

      _RL fac
      _RL topomin

      character*(128) fname1, fname2, fname3

      logical exst
CEOP

      call ecco_zero(localbar,Nr,zeroRL,myThid)
      call ecco_zero(localweight,Nr,zeroRL,myThid)
      call ecco_zero(localtmp,Nr,zeroRL,myThid)
      call ecco_zero(localmask,Nr,zeroRL,myThid)

      call ecco_zero(localobs,Nr,zeroRL,myThid)
      call ecco_zero(localdif,Nr,zeroRL,myThid)
      call ecco_zero(difmask,Nr,zeroRL,myThid)

      call ecco_zero(localdifmsk,Nr,zeroRL,myThid)
      call ecco_zero(localdifsum,Nr,zeroRL,myThid)
      call ecco_zero(localdifnum,Nr,zeroRL,myThid)

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
cinitialize:
c      nxb=sNx/MSfacX
c      nyb=sNy/MSfacY
      tmpind1=sNx/MSfacX
      tmpind2=sNy/MSfacY
c      call ecco_zero_lo(localmasklo,nxb,nyb,Nr,
      call ecco_zero_lo(localmasklo,tmpind1,tmpind2,Nr,
     &                  zeroRL,myThid)
catn because we have to do the cost_lo outside the irec
catn loop, we will make the weightlo of same size as
catn the localobslo
cc      call ecco_zero_lo(localweightlo,nxb,nyb,Nr,
c      call ecco_zero_lo(localweightlo,tmpind1,tmpind2,Nr,
c     &                  zeroRL,myThid)
#endif
#endif

      dosumsq=.TRUE.
      domean=.FALSE.
      doanom=.FALSE.
      dovarwei=.FALSE.
      dosmooth=.FALSE.
      dooffset=.FALSE.
      domaskmindepth=.FALSE.
      k2smooth=1
      doclim=.FALSE.
      nrecclim=nrecloc
      fac=oneRL

      do k2 = 1, NGENPPROC
          if (preproc(k2).EQ.'mean') domean=.TRUE.
          if (preproc(k2).EQ.'anom') doanom=.TRUE.
          if (preproc(k2).EQ.'variaweight') dovarwei=.TRUE.
          if (preproc(k2).EQ.'nosumsq') dosumsq=.FALSE.
          if (preproc(k2).EQ.'offset') dooffset=.TRUE.
          if (preproc(k2).EQ.'mindepth') then
           domaskmindepth=.TRUE.
           topomin = preproc_r(k2)
          endif
          if (posproc(k2).EQ.'smooth') then
            dosmooth=.TRUE.
            k2smooth=k2
          endif
          if (preproc(k2).EQ.'clim') then
            doclim=.TRUE.
            nrecclim=preproc_i(k2)
          endif
          if (preproc(k2).EQ.'factor') then
            fac=preproc_r(k2)
          endif
      enddo

      write(*,*) 'inside cost_genloopA ',
     &'do[mean,anom,varwei,sumsq,offset,smooth,clim]: ',
     & domean,doanom,dovarwei,dosumsq,dooffset,dosmooth,doclim
      write(*,*) 'inside cost_genloopB ',
     &'domaskmindepth,topomin,k2smooth,nrecclim,fac: ',
     & domaskmindepth,topomin,k2smooth,nrecclim,fac
      write(*,*) 'inside cost_genloopC ',
     & 'outlev,ylocmask,addVariaCost',
     &  outlev,ylocmask,addVariaCost

C--   Assign mask, write explicitly to include overlaps and
C     avoid exchange later (e.g. in smooth routine), note that
C     the routine ecco_cprsrl does not copy overlaps
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

        IF ( ylocmask .EQ. 'C' .OR. ylocmask .EQ. 'c' ) THEN
         DO k = 1, Nr
          DO j = 1-OLy,sNy+OLy
           DO i = 1-OLx,sNx+OLx
            localmask(i,j,k,bi,bj) = maskC(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO

        ELSEIF ( ylocmask .EQ. 'S' .OR. ylocmask .EQ. 's' ) THEN

         DO k = 1, Nr
          DO j = 1-OLy,sNy+OLy
           DO i = 1-OLx,sNx+OLx
            localmask(i,j,k,bi,bj) = maskS(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO

        ELSEIF ( ylocmask .EQ. 'W' .OR. ylocmask .EQ. 'w' ) THEN

         DO k = 1, Nr
          DO j = 1-OLy,sNy+OLy
           DO i = 1-OLx,sNx+OLx
            localmask(i,j,k,bi,bj) = maskW(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO

        ELSE
         CALL ECCO_ERROR('COST_GENLOOP', 'wrong ylocmask', myThid)
        ENDIF

#ifdef ALLOW_SMOOTH
C       Set 2D Mask for 2D smooth routine
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
          mask2D(i,j,bi,bj) = localmask(i,j,1,bi,bj)
          tmpFld2D(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif
C-- end bi/bj loop for setting masks
       ENDDO
      ENDDO

c-- set nrecloop to nrecloc
      nrecloop=nrecloc

c-- reset nrecloop, if needed, according to preproc
      if ( doclim ) nrecloop=MIN(nrecloop,nrecclim)

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
cinitialize:
c      nxb=sNx/MSfacX
c      nyb=sNy/MSfacY
      tmpind1=sNx/MSfacX
      tmpind2=sNy/MSfacY
c      tmpind3=nrecloop*nnzobs
c the way this is coded, nnzobs is either 1 or Nr
c so best to initialize to the full size for now
c even though the full size, i initiated as 365*Nr
c which, if we have daily data for more than 1 yr
c and if data fully 3d, then we're underestimating.
c That is right now it works with nrecloop<=365
c If nnzobs is 1, then we're right now accomodating 21yrs
c If nnzobs is Nr, then only 1yr of data is accom. now.
      tmpind3=nrecloop*Nr
c      call ecco_zero_lo(localdiflo,nxb,nyb,tmpind3,
      call ecco_zero_lo(localdiflo,tmpind1,tmpind2,tmpind3,
     &                  zeroRL,myThid)
c      call ecco_zero_lo(localbarlo,nxb,nyb,tmpind3,
      call ecco_zero_lo(localbarlo,tmpind1,tmpind2,tmpind3,
     &                  zeroRL,myThid)
c      call ecco_zero_lo(localobslo,nxb,nyb,tmpind3,
      call ecco_zero_lo(localobslo,tmpind1,tmpind2,tmpind3,
     &                  zeroRL,myThid)
c      call ecco_zero_lo(localweightlo,nxb,nyb,tmpind3,
      call ecco_zero_lo(localweightlo,tmpind1,tmpind2,tmpind3,
     &                  zeroRL,myThid)
c the mean field is mean over the nrecloop, so max size is Nr
      if( domean )
c     &call ecco_zero_lo(localdifMeanOutlo,nxb,nyb,Nr,
     &call ecco_zero_lo(localdifMeanOutlo,tmpind1,tmpind2,Nr,
     &                  zeroRL,myThid)
c      do irec = 1,nrecloop     !{{{
c       do k=1,nnzobs
c        tmpind3=(irec-1)*nnzobs+k
c        do j=1,nyb
c         do i=1,nxb
c          localdiflo(i,j,tmpind3,bi,bj)=0. _d 0
c          localbarlo(i,j,tmpind3,bi,bj)=0. _d 0
c          localobslo(i,j,tmpind3,bi,bj)=0. _d 0
c         enddo
c        enddo
c       enddo
c      enddo                    !}}}
c      DO bj = myByLo(myThid), myByHi(myThid)   !{{{
c       DO bi = myBxLo(myThid), myBxHi(myThid)
c        do k=1,nnzobs
c         do j=1,nyb
c          do i=1,nxb
ccatn below: should be inside loop domean only
c            localdifMeanOutlo(i,j,k,bi,bj)=0. _d 0
c          enddo
c         enddo
c        enddo
c       enddo
c      enddo                                    !}}}
catn now set mask: localmask is simply either grid msk[C,W,S]
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        do k=1,Nr
c         do j=1,nxb
c          do i=1,nyb
         do j=1,sNy/MSfacY
          do i=1,sNx/MSfacX
           do l=1,MSfacY
            do m=1,MSfacX
             tmpind1 = (i-1)*MSfacX+m
             tmpind2 = (j-1)*MSfacY+l
catn: very crude binary 0/1
             localmasklo(i,j,k,bi,bj) = max( localmasklo(i,j,k,bi,bj)
     &                            ,localmask(tmpind1,tmpind2,k,bi,bj) )
            enddo
           enddo
          enddo
         enddo
        enddo
       ENDDO
      ENDDO
#endif /* MULTISCALE_COUPLING_TAPES */
#endif /* ALLOW_AUTODIFF */

c--   loop over obsfile records
      do irec = 1, nrecloop

c--     load weights
        exst=.FALSE.
        jrec=1
        if( dovarwei ) jrec = irec
        call cost_gencal(localbarfile, localerrfile,
     &     jrec, localstartdate, localperiod, fname1,
     &     fname3, localrec, obsrec, exst, myThid )
        call ecco_zero(localweight,Nr,zeroRL,myThid)
        if ( (localrec .GT. 0).AND.(obsrec .GT. 0).AND.(exst) )
catn: here we read in weight, either Nr or 1 
     &  call ecco_readwei(fname3,localweight,
     &                    localrec,Nr,nnzobs,dosumsq,myThid)

c--     determine records and file names
        exst=.FALSE.
        call cost_gencal(localbarfile, localobsfile,
     &     irec, localstartdate, localperiod, fname1,
     &     fname2, localrec, obsrec, exst, myThid )

c--     load model average and observed average
        call ecco_zero(localbar,Nr,zeroRL,myThid)
        call cost_genread( fname1, localbar, localtmp, irec, Nr, nnzbar,
     &       nrecloc, preproc, preproc_c, preproc_i, preproc_r,
     &       dummy, myThid )
        call ecco_mult(localbar,fac,Nr,nnzbar,myThid)

        call ecco_zero(localobs,Nr,spzeroloc,myThid)
        if ( (localrec .GT. 0).AND.(obsrec .GT. 0).AND.(exst) ) then
         CALL READ_REC_LEV_RL( fname2, cost_iprec, Nr, 1, nnzobs,
     &                         localobs, localrec, 1, myThid )
        endif

c--     Compute masked model-data difference
        call ecco_diffmsk( localbar, localobs, localmask,
     &     Nr, nnzobs, spminloc, spmaxloc, spzeroloc,
     &     localdif, difmask, myThid )

c#ifdef ALLOW_AUTODIFF          !{{{
c#ifdef MULTISCALE_COUPLING_TAPES
ccatn there is an assumption in ecco_difmask that nnzobs<=Nr
c        DO bj = myByLo(myThid), myByHi(myThid)
c         DO bi = myBxLo(myThid), myBxHi(myThid)
c          do k=1,nnzobs
cc tmpind3 is shifted by nnzobs=Nr or by nnzobs=1
c           tmpind3 = (irec-1)*nnzobs+k
c           do j=1,nyb
c            do i=1,nxb
c             difmasklo(i,j,k,bi,bj)=localmasklo(i,j,k,bi,bj)
c             tmpcount=0
c             do l=1,MSfacY
c              do m=1,MSfacX
c               tmpind1 = (i-1)*MSfacX+m
c               tmpind2 = (j-1)*MSfacY+l
ccatn we have a choice here of how to treat mask:
ccatn whether if there is any valid data we keep, or
ccatn if there is any invalid data we omit all.
ccatn For now, try 2nd method, but leave code for 1st method as well 
c               if(
c     &      localobs(tmpind1,tmpind2,k,bi,bj) .LT. spMinLoc .OR.
c     &      localobs(tmpind1,tmpind2,k,bi,bj) .GT. spMaxLoc .OR.
c     &      localobs(tmpind1,tmpind2,k,bi,bj) .EQ. spzeroLoc ) then
cc                difmasklo(i,j,k,bi,bj) = 0. _d 0
c                tmpcount=tmpcount+1
ccatn crude treatment of weight: take max; also weightlo dim is tmpind3
c                localweightlo(i,j,tmpind3,bi,bj)=max(
c     &            localweigthlo(i,j,tmpind3,bi,bj)
c     &          , localweight(tmpind1,tmpind2,k,bi,bj) )
c                localobslo(i,j,tmpind3,bi,bj)=
c     &          localobslo(i,j,tmpind3,bi,bj)+
c     &          localobs(tmpind1,tmpind2,k,bi,bj)
c               endif
c              enddo
c             enddo
c             if(tmpcount.eq.0) then
c              difmasklo(i,j,k,bi,bj) = 0. _d 0
c             elseif((tmpcount.gt.0.and.tmpcount.le.MSfacX*MSfacY).and.
c     &           difmasklo(i,j,k,bi.bj).gt.0. _d 0) then
c              localobslo(i,j,tmpind3,bi,bj)=
c     &        localobslo(i,j,tmpind3,bi,bj)/tmpcount
c             endif
c            enddo
c           enddo
c          enddo
c         enddo
c        enddo
c#endif /* MULTISCALE_COUPLING_TAPES */
c#endif /* ALLOW_AUTODIFF */             !}}}

        if ( domaskmindepth ) call ecco_maskmindepth(
     &       difmask, Nr, nnzobs, topomin, myThid )

        if ( doanom ) call ecco_subtract( localdif,
     &       localdifmeanIn, Nr, nnzobs, myThid )

        if ( domean.OR.doanom ) call ecco_addmask(
     &       localdif, difmask, localdifsum,
     &       localdifnum, Nr, nnzobs, myThid )

catn do the interp to lores before or after smooth?      !{{{
catn ideally we want to avoid smooth because lores 
catn is simply the binavg of the hires, so what we
catn want is for hires, we compute cost from smooth
catn while for lores, we compute cost from binavg
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
catn        write(*,*) 'inside cost_genloopE '
c       if(useSingleCPUio) then
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          do k=1,nnzobs                 !1
           tmpind3 = (irec-1)*nnzobs+k  !1:nrec for nz=nnzobs=1
c           do j=1,nyb                   !1:16
c            do i=1,nxb                  !1:20
           do j=1,sNy/MSfacY            !1:16
            do i=1,sNx/MSfacX           !1:20
             difmasklo(i,j,tmpind3,bi,bj)=localmasklo(i,j,k,bi,bj)
             tmpcount=0
             do l=1,MSfacY              !1:4
              do m=1,MSfacX             !1:4
               tmpind1 = (i-1)*MSfacX+m !1:64
               tmpind2 = (j-1)*MSfacY+l !1:80

catn we have a choice here of how to treat mask:
catn whether if there is any valid data we keep, or
catn if there is any invalid data we omit all.
catn For now, try 2nd method, but leave code for 1st method as well 
catn crude treatment of weight: take max
               if(
     &    localobs(tmpind1,tmpind2,k,bi,bj) .LT. spMinLoc .OR.
     &    localobs(tmpind1,tmpind2,k,bi,bj) .GT. spMaxLoc .OR.
     &    localobs(tmpind1,tmpind2,k,bi,bj) .EQ. spzeroLoc ) then
               else
catn method1 ---
c               difmasklo(i,j,k,bi,bj) = 0. _d 0
catn method1 ---
                tmpcount=tmpcount+1
cmodel:
                localbarlo(i,j,tmpind3,bi,bj)=
     &          localbarlo(i,j,tmpind3,bi,bj)+
     &          localbar(tmpind1,tmpind2,k,bi,bj)
cobs;
                localobslo(i,j,tmpind3,bi,bj)=
     &          localobslo(i,j,tmpind3,bi,bj)+
     &          localobs(tmpind1,tmpind2,k,bi,bj)
cdiff:
                localdiflo(i,j,tmpind3,bi,bj)=
     &          localdiflo(i,j,tmpind3,bi,bj)+
     &          localdif(tmpind1,tmpind2,k,bi,bj)

catn crude treatment of weight: take max; also weightlo dim is tmpind3
                localweightlo(i,j,tmpind3,bi,bj) = max(
     &             localweightlo(i,j,tmpind3,bi,bj) ,
     &             localweight(tmpind1,tmpind2,k,bi,bj) )
               endif

catn               if(irec.ge.45 .and. irec.le. 48) then
catn                write(*,'(9(I5,1X),4(1x,F12.6))') 
catn     &           irec,i,j,tmpcount,l,m,tmpind1,tmpind2,tmpind3,
catn     &           localdiflo(i,j,tmpind3,bi,bi),
catn     &           localdif(tmpind1,tmpind2,k,bi,bj),
catn     &           localbarlo(i,j,tmpind3,bi,bi),
catn     &           localbar(tmpind1,tmpind2,k,bi,bj)
catn               endif
              enddo!m
             enddo !l

             if( tmpcount.lt.1 ) then
              difmasklo(i,j,tmpind3,bi,bj) = 0. _d 0
             elseif((tmpcount.gt.0 .and. tmpcount.le.MSfacX*MSfacY)
     &         .and. (difmasklo(i,j,tmpind3,bi,bj) .gt. 0. _d 0) ) then
              localobslo(i,j,tmpind3,bi,bj)=
     &        localobslo(i,j,tmpind3,bi,bj)/tmpcount
              localbarlo(i,j,tmpind3,bi,bj)=
     &        localbarlo(i,j,tmpind3,bi,bj)/tmpcount
              localdiflo(i,j,tmpind3,bi,bj)=
     &        localdiflo(i,j,tmpind3,bi,bj)/tmpcount
             endif

            enddo !i
           enddo  !j
          enddo   !k
         enddo    !bi
        enddo     !bj
#endif /* MULTISCALE_COUPLING_TAPES */
#endif /* ALLOW_AUTODIFF */
c!}}}

        
        if (addVariaCost) then

#ifdef ALLOW_SMOOTH
catn note smooth is only done to 2d obs field
         IF ( useSMOOTH.AND.dosmooth.AND.
     &     (nnzbar.EQ.1).AND.(nnzobs.EQ.1) ) THEN
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j = 1-OLy,sNy+OLy
             DO i = 1-OLx,sNx+OLx
              tmpFld2D(i,j,bi,bj) = localdif(i,j,1,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
catn          CALL SMOOTH_HETERO2D(tmpFld2D,mask2D,
catn     &         posproc_c(k2smooth),posproc_i(k2smooth),myThid)
         call smooth2d(tmpFld2D,mask2D,posproc_i(k2smooth),myThid)
         DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j = 1-OLy,sNy+OLy
             DO i = 1-OLx,sNx+OLx
              localdif(i,j,1,bi,bj) = tmpFld2D(i,j,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
#endif

c--     Compute normalized model-obs cost function
         call ecco_addcost(
     I             localdif, localweight, difmask, Nr, nnzobs, dosumsq,
     U             objf_local, num_local,
     I             myThid
     &             )
c--     output model-data difference to disk
         if ( outlev.GT.0 ) then
          il=ilnblnk(outname)
          write(fname3(1:128),'(2a)') 'misfit_', outname(1:il)
          if ( nnzobs.EQ.1 ) CALL
     &    WRITE_REC_XY_RL( fname3, localdif,irec, eccoiter, myThid )
          if ( nnzobs.EQ.Nr ) CALL
     &    WRITE_REC_XYZ_RL( fname3, localdif,irec, eccoiter, myThid )
         endif !outlev.gt.0

         if ( outlev.GT.1 ) then
          if ( dovarwei ) then
           il=ilnblnk(outname)
           write(fname3(1:128),'(2a)') 'weight_', outname(1:il)
           if ( nnzobs.EQ.1 ) CALL WRITE_REC_XY_RL( fname3,
     &                           localweight, irec, eccoiter, myThid )
           if ( nnzobs.EQ.Nr ) CALL WRITE_REC_XYZ_RL( fname3,
     &                           localweight, irec, eccoiter, myThid )
          endif
         endif
        endif  !addvariacost

      enddo     !irec
c--   End of loop over obsfile records.

      call ecco_zero(localdifmeanOut,Nr,zeroRL,myThid)
      call ecco_cp      (localdifsum,localdifmeanOut,Nr,nnzobs,myThid)
      call ecco_divfield(localdifmeanOut,localdifnum,Nr,nnzobs,myThid)
      call ecco_cp      (localdifnum,localdifmsk,Nr,nnzobs,myThid)
      call ecco_divfield(localdifmsk,localdifnum,Nr,nnzobs,myThid)

      if (dooffset) then
       il=ilnblnk(localbarfile)
       call ecco_offset(
     I      localbarfile(1:il),
     U      localdifmeanOut,
     I      localdifmsk,Nr,nnzobs,myThid)
      endif

      if ( domean ) then
c--     Compute normalized model-obs cost function
        call ecco_addcost(
     I                    localdifmeanOut, localweight, localdifmsk,
     I                    Nr, nnzobs, dosumsq,
     U                    objf_local, num_local,
     I                    myThid )

c--     output model-data difference to disk
        if ( outlev.GT.0 ) then
          il=ilnblnk(outname)
          write(fname3(1:128),'(2a)') 'misfit_', outname(1:il)
          if ( nnzobs.EQ.1 ) CALL
     &    WRITE_REC_XY_RL(fname3,localdifmeanOut,1,eccoiter,myThid)
          if ( nnzobs.EQ.Nr ) CALL
     &    WRITE_REC_XYZ_RL(fname3,localdifmeanOut,1,eccoiter,myThid)
        endif
#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
        if ( outlev.GT.0 ) then
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            do k=1,nnzobs
c             do j=1,nyb
c              do i=1,nxb
             do j=1,sNy/MSfacY
              do i=1,sNx/MSfacX
               do l=1,MSfacY
                do m=1,MSfacX
                 tmpind1 = (i-1)*MSfacX+m
                 tmpind2 = (j-1)*MSfacY+l
                 localdifmeanOutlo(i,j,k,bi,bj)=
     &           localdifmeanOutlo(i,j,k,bi,bj)+
     &           localdifmeanOut(tmpind1,tmpind2,k,bi,bj)/MSfacX/MSfacY
                enddo
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo

        endif   !outlev.gt.0
#endif /* MULTISCALE_COUPLING_TAPES */
#endif /* ALLOW_AUTODIFF */
      endif     !domean

#ifdef ALLOW_AUTODIFF
#ifdef MULTISCALE_COUPLING_TAPES
      if (addVariaCost) then
         write(*,*) 'inside cost_genloop D'
         write(*,*) 'localdiflo ' !,localdiflo(15,14)
         write(*,*) 'localweightlo ' !,localweightlo(15,14)
         write(*,*) 'difmasklo ' !,difmasklo(15,14)
         write(*,*) 'objf_local_lo: ',objf_local_lo
         write(*,*) 'num_local_lo: ',num_local_lo
         call ecco_addcost_lo(
     I             localdiflo, localweightlo, difmasklo,
     I             sNx/MSfacX,sNy/MSfacY,nrecloop*nnzobs, dosumsq,
     U             objf_local_lo, num_local_lo,
     I             myThid
     &             )
c     I             nxb,nyb,nrecloop*nnzobs, dosumsq,
      endif
         write(*,*) 'objf_local_lo: ',objf_local_lo
         write(*,*) 'num_local_lo: ',num_local_lo
c      if ( outlev.GT.0 ) then
          il=ilnblnk(outname)
          write(fname3(1:128),'(3a)') 'misfit_', outname(1:il), '_lo'
        if( domean) then
          nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*nnzobs
c          nArr=nxb*nyb*nSx*nSy*nnzobs
          call write_glvec_rl(fname3,' ',
     I          localdifmeanOutlo,nArr,myIter,myThid)
        else
c          nArr=nxb*nyb*nSx*nSy*nrecloop*nnzobs
          nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*nrecloop*nnzobs
          call write_glvec_rl(fname3,' ',
     I          localdiflo,nArr,myIter,myThid)
        endif
        write(*,*) nArr,domean,fname3
        il=ilnblnk(localbarfile)
        write(fname3(1:128),'(2a)') localbarfile(1:il), '_lo'
        nArr=sNx*sNy*nSx*nSy/MSfacX/MSfacY*nrecloop*nnzobs
c        nArr=nxb*nyb*nSx*nSy*nrecloop*nnzobs
        call write_glvec_rl(fname3,' ',
     I        localbarlo,nArr,myIter,myThid)
ccatn: now also attempt to write meta file!
c        dimList(1,1) = nxb
c        dimList(2,1) = 1
c        dimList(3,1) = nxb
c        dimList(1,2) = nyb
c        dimList(2,2) = 1
c        dimList(3,2) = nyb
c        dimList(1,3) = nnzobs
c        dimList(2,3) = 1
c        dimList(3,3) = nnzobs
c        if(nnzobs .eq. 1) then
c          nDims=2
c        else
c          nDims=3
c        endif
c        MSnf = nrecloop
c        map2gl(1)=0
c        map2gl(2)=1
c        CALL MDS_WRITE_META(
c     I              metaFName, dataFName, the_run_name, ' ',
c     I              filePrec, nDims, dimList, map2gl, MSnf, blank8c,
c     I              0, dummyRL, oneRL, iRec, myIter, myThid )
c
c      endif
catn make an attempt to write local cost, not sure if this make sense...
CSR: doesn't work with mpi

c      il=ilnblnk(outname)
c      write(fname3(1:128),'(2a)') 'costfinallo_',outname(1:il)
c      il=ilnblnk(fname3)
c      open(0,file=fname3(1:il),status='new')
c      write(0,'(G18.9,1x,F18.2)') objf_local_lo,num_local_lo
c      close(0)

#endif /* MULTISCALE_COUPLING_TAPES */
#endif /* ALLOW_AUTODIFF */
c      if ( outlev.GT.1 ) then
       if ( .NOT.dovarwei ) then
          il=ilnblnk(outname)
          write(fname3(1:128),'(2a)') 'weight_', outname(1:il)
          if ( nnzobs.EQ.1 ) CALL
     &    WRITE_REC_XY_RL( fname3, localweight, 1, eccoiter, myThid )
          if ( nnzobs.EQ.nr ) CALL
     &    WRITE_REC_XYZ_RL( fname3, localweight, 1, eccoiter, myThid )
      endif
c      endif


#endif /* ALLOW_ECCO */

      RETURN
      END
