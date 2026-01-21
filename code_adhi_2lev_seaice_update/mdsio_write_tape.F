#include "MDSIO_OPTIONS.h"
catn
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h" /* for multiscale use */
#endif
catn
CAV(
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h" /* for multiscale use */
#endif
CAV)

CBOP
C !ROUTINE: MDS_WRITE_TAPE
C !INTERFACE:
      SUBROUTINE MDS_WRITE_TAPE(
     I   fName,
     I   filePrec,
     I   globalfile,
     I   arrType,
     I   nSize,
     I   fldR8, fldR4,
     I   singleCpuIO,
     I   iRec,
     I   myIter,
     I   myThid )

C !DESCRIPTION:
C MDS_WRITE_TAPE: write an array (treated as vector) to a tape-file
C  (renamed from MDSWRITEVECTOR with 2 explicit input array types)
C
C Arguments:
C fName      string  :: base name for file to write
C filePrec   integer :: number of bits per word in file (32 or 64)
C globalFile logical :: selects between writing a global or tiled file
C arrType    char(2) :: which array (fldR8/R4) to write, either "R8" or "R4"
C nSize      integer :: number of elements of input array "fldR8/R4" to write
C fldR8      ( R8 )  :: array to write if arrType="R8", fldR8(nSize)
C fldR4      ( R4 )  :: array to write if arrType="R4", fldR4(nSize)
C bi,bj      integer :: tile indices (if tiled array)
C singleCpuIO ( L )  :: only proc 0 do IO and collect data from other procs
C iRec       integer :: record number to write
C myIter     integer :: time step number
C myThid     integer :: my Thread Id number

C !USES:
      IMPLICIT NONE

C-- Global variables --
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
CAV(
#ifdef MULTISCALE_COUPLING_TAPES
#include "GRID.h" /* for hFacC */
#ifdef ALLOW_SEAICE
#include "SEAICE_OPTIONS.h"
#include "SEAICE_SIZE.h"
#endif
#include "AUTODIFF_PARAMS.h" /* for scaling factors */
#include "AUTODIFF.h"        /* for tape level sizes */
#include "AUTODIFF_STORE.h"        /* for tape level sizes */
#ifdef ALLOW_CTRL
#include "CTRL_OPTIONS.h"
#include "CTRL_SIZE.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h" /* for PTRACERS_PARAMS */
# include "PTRACERS_PARAMS.h" /* for PTRACERS_numInUse */
#endif
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif
#include "EEBUFF_SCPU.h" /* for gather_2d_r8 */
#endif
CAV)

C !INPUT/OUTPUT PARAMETERS:
      CHARACTER*(*) fName
      INTEGER filePrec
      LOGICAL globalfile
      CHARACTER*(2) arrType
      INTEGER nSize
      _R8     fldR8(*)
      _R4     fldR4(*)
      LOGICAL singleCpuIO
      INTEGER iRec
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_AUTODIFF

C !FUNCTIONS:
      INTEGER ILNBLNK
      INTEGER MDS_RECLEN
      EXTERNAL ILNBLNK
      EXTERNAL MDS_RECLEN

C !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_FNAM) dataFName, metaFName, pfName
      INTEGER iG, jG, jRec, dUnit, IL, pIL
      INTEGER dimList(3,3), nDims, map2gl(2)
      INTEGER length_of_rec
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C simple implementation of singleCpuIO without any specific EXCH2
C feature (should work as long as reading and writing match)
      INTEGER j
      INTEGER vec_size
C Note: would be better to use explicit (allocate/delocate) dynamical
C       allocation instead of this implicit form:
      _R8    gl_buffer_r8(nSize*nPx*nPy)
      _R4    gl_buffer_r4(nSize*nPx*nPy)
      _R8    local_r8    (nSize)
      _R4    local_r4    (nSize)
      _RL dummyRL(1)
      CHARACTER*8 blank8c

CAV(
#ifdef MULTISCALE_COUPLING_TAPES
      INTEGER   nxt,nyt,nxb,nyb
      INTEGER   MSvec_size, nbin, mmirec
      INTEGER   bi,bj,MSai,MSaj,MSak,MSaf,MSpx,MSpy,MSnz,MSnf
      INTEGER   tOLx, tOLy, tmpind1, tmpind2, tbi, tbj
      _R8       tmp_bin_r8
      _R4       tmp_bin_r4
#ifdef ALLOW_EXCH2
      _R8    fld_2d_r8(exch2_global_Nx*exch2_global_NY/MSfacX/MSfacY)
      _R4    fld_2d_r4(exch2_global_Nx*exch2_global_NY/MSfacX/MSfacY)
      _R8 fld_3d_r8(exch2_global_Nx*exch2_global_NY*Nr/MSfacX/MSfacY)
      _R4 fld_3d_r4(exch2_global_Nx*exch2_global_NY*Nr/MSfacX/MSfacY)
      _R8    fld_obx2_r8(exch2_global_Nx/MSfacX)
      _R4    fld_obx2_r4(exch2_global_Nx/MSfacX)
      _R8    fld_obx3_r8(exch2_global_Nx*Nr/MSfacX)
      _R4    fld_obx3_r4(exch2_global_Nx*Nr/MSfacX)
      _R8    fld_oby2_r8(exch2_global_Ny/MSfacY)
      _R4    fld_oby2_r4(exch2_global_Ny/MSfacY)
      _R8    fld_oby3_r8(exch2_global_Ny*Nr/MSfacY)
      _R4    fld_oby3_r4(exch2_global_Ny*Nr/MSfacY)
#else
      _R8    fld_2d_r8(NX*NY/MSfacX/MSfacY)
      _R4    fld_2d_r4(NX*NY/MSfacX/MSfacY)
      _R8    fld_3d_r8(NX*NY*Nr/MSfacX/MSfacY)
      _R4    fld_3d_r4(NX*NY*Nr/MSfacX/MSfacY)
      _R8    fld_obx2_r8(Nx/MSfacX)
      _R4    fld_obx2_r4(Nx/MSfacX)
      _R8    fld_obx3_r8(Nx*Nr/MSfacX)
      _R4    fld_obx3_r4(Nx*Nr/MSfacX)
      _R8    fld_oby2_r8(Ny/MSfacY)
      _R4    fld_oby2_r4(Ny/MSfacY)
      _R8    fld_oby3_r8(Ny*Nr/MSfacY)
      _R4    fld_oby3_r4(Ny*Nr/MSfacY)
#endif
      LOGICAL   useExch2ioLayOut, zeroBuff
      INTEGER   xSize, ySize
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      LOGICAL   IsOBCN, IsOBCE,IsOBCS,IsOBCW
#endif
      LOGICAL   Ifill
#endif
CAV)

CEOP

      DATA dummyRL(1) / 0. _d 0 /
      DATA blank8c / '        ' /
      DATA map2gl  / 0, 1 /

CAV(
#ifdef MULTISCALE_COUPLING_TAPES
      zeroBuff = .FALSE.
      useExch2ioLayOut = .FALSE.
      xSize = Nx
      ySize = Ny
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      IsOBCN = .FALSE.
      IsOBCE = .FALSE.
      IsOBCS = .FALSE.
      IsOBCW = .FALSE.
#endif
#ifdef ALLOW_EXCH2
      zeroBuff = .TRUE.
      useExch2ioLayOut = .TRUE.
      xSize = exch2_global_Nx
      ySize = exch2_global_Ny
#endif 
#endif
CAV)

      vec_size = nSize*nPx*nPy

C--   Copy input array to local buffer
        IF ( arrType.EQ.'R4' ) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
            DO j=1,nSize
              local_r4(j) = fldR4(j)
            ENDDO
          ELSE
            DO j=1,nSize
              local_r8(j) = fldR4(j)
            ENDDO
          ENDIF
        ELSEIF ( arrType.EQ.'R8' ) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
            DO j=1,nSize
              local_r4(j) = fldR8(j)
            ENDDO
          ELSE
            DO j=1,nSize
              local_r8(j) = fldR8(j)
            ENDDO
          ENDIF
        ELSE
          WRITE(msgBuf,'(A)')
     &         ' MDS_WRITE_TAPE: illegal value for arrType'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
        ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Only do I/O if I am the master thread
      _BEGIN_MASTER( myThid )

C-    Record number must be >= 1
      IF ( iRec.LT.1 ) THEN
       WRITE(msgBuf,'(A,I10)')
     &   ' MDS_WRITE_TAPE: argument iRec =',iRec
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' MDS_WRITE_TAPE: invalid value for iRec'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
      ENDIF

C-    Assume nothing
      IL  = ILNBLNK( fName )
      pIL = ILNBLNK( mdsioLocalDir )

C-    Assign special directory
      IF ( pIL.EQ.0 ) THEN
        pfName = fName
      ELSE
        WRITE(pfName,'(2A)') mdsioLocalDir(1:pIL), fName(1:IL)
      ENDIF
      pIL = ILNBLNK( pfName )
      IF ( debugLevel.GE.debLevC .AND.
     &     ( .NOT.singleCpuIO .OR. myProcId.EQ.0 ) ) THEN
        WRITE(msgBuf,'(A,I8,2A)')
     &      ' MDS_WRITE_TAPE: iRec=', iRec, ', file=', pfName(1:pIL)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

C-    Assign a free unit number as the I/O channel for this routine
      CALL MDSFINDUNIT( dUnit, myThid )

CAV(
#ifdef MULTISCALE_COUPLING_TAPES

      IF ( singleCpuIO ) THEN

cmm      print*,'CMM, pfname: ', pfName(pIL-4-7:pIL-2-7)

CMM     HERE FIND OUT IF ITS IsOBC
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      IF (       pfName(pIL-4-7:pIL-7).EQ.'obcsn')
     &      IsOBCN = .TRUE.
      IF (       pfName(pIL-4-7:pIL-7).EQ.'obcse')
     &      IsOBCE = .TRUE.
      IF (       pfName(pIL-4-7:pIL-7).EQ.'obcss')
     &      IsOBCS = .TRUE.
      IF (       pfName(pIL-4-7:pIL-7).EQ.'obcsw')
     &      IsOBCW = .TRUE.
CMM
            DO j=1,8
      IF (       pfName(pIL-8-j:pIL-6-j).EQ.'obn') 
     &      IsOBCN = .TRUE.
      IF (       pfName(pIL-8-j:pIL-6-j).EQ.'obe')
     &      IsOBCE = .TRUE.
      IF (       pfName(pIL-8-j:pIL-6-j).EQ.'obs')
     &      IsOBCS = .TRUE.
      IF (       pfName(pIL-8-j:pIL-6-j).EQ.'obw')
     &      IsOBCW = .TRUE.
            ENDDO
#endif
CMM)
cmm        IF (myProcId.EQ.44)   print*,'CMMa: fname ',pfName(pIL-20:pIL)

CAV dimensions of each tile, including overlap
CMM: SET TO DEFAULT SIZE VALUES. WILL OVERWRITE FOR UNIQUE CASES
      MSnz = 1
      MSnf = 1
      nxt = sNx+2*OLx
      nyt = sNy+2*OLy
CAV dimensions of the binned array
      nxb = xSize/MSfacX
      nyb = ySize/MSfacY

CMM(
C  NEED TO KILL THIS IF NR IS EQUAL TO  ANY THING ELSE
C:  NDV2D (~23) NEXF1 (~23) NEXF2 (~24) NCTRL1 (~20) NSI (~19) maxCtrlTim2D (~14)
C: CANNOT BE SAME AS
C: NR or Nr*2 or 
C: Nr*NDV3D (NDV3d is typically 10 or 14)
C: Nr*ptracers_num (ptracers_num is typically about 8)
C: This should not be a problem for large NR
C:  but will be if NR is less than 24. 
c: likely will be a problem for tuturials. PUT IN A STOP HERE 
CMM)
CMM(
      IF ((NR.EQ.NDV2D).OR.(2*NR.EQ.NDV2D).OR.
     &    (NR.EQ.NEXF1).OR.(2*NR.EQ.NEXF1).OR.
     &    (NR.EQ.NEXF2).OR.(2*NR.EQ.NEXF2).OR.
     &    (NR.EQ.NCTRL1).OR.(2*NR.EQ.NCTRL1)
#if (defined (ALLOW_GENARR2D_CONTROL) || defined (ALLOW_GENARR3D_CONTROL) || defined (ALLOW_GENTIM2D_CONTROL))
     &.OR.(NR.EQ.maxCtrlTim2D).OR.(2*NR.EQ.maxCtrlTim2D)
#endif
catn add loop for sea ice
#ifdef ALLOW_SEAICE
     &.OR.(NR.EQ.NSI).OR.(2*NR.EQ.NSI) 
#endif
     &  ) THEN
         print *,'CMM: AMBIGUITY. MULITSCALE WONT WORK'
        WRITE(msgBuf,'(A,I10)')
     &   ' MDS_WRITE_TAPE: argument nr = ', nr
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)')
     &   ' MDS_WRITE_TAPE: invalid value for nr when using multiscale'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
      ENDIF
CMM)

catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      IF (IsOBCN .OR. IsOBCS) THEN
          nyt = 1
CMM: IS IT A 3D OPEN N or S BOUNDARY
        IF (nSize.eq.(sNx+2*OLx)*nSx*nSy*Nr*NOB) THEN
          MSnz = NR
          MSnf = NOB
        ENDIF
      ELSEIF (IsOBCE .OR. IsOBCW) THEN
          nxt = 1
CMM: IS IT A 3D OPEN E or W BOUNDARY
       IF (nSize.eq.(sNy+2*OLy)*nSx*nSy*Nr*NOB) THEN
          MSnz = NR
          MSnf = NOB
       ENDIF
#else
      IF (.FALSE.) THEN
#endif
      ELSE
CMM: NOT AN OBC. IS IT A 3D FIELD
       IF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*Nr) THEN
          MSnz=Nr
#ifdef ALLOW_PTRACERS
CMM: IS IT PTRACERS
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*Nr*
     &                ptracers_num) THEN
          MSnz=Nr
          MSnf=ptracers_num
#endif
CMM: IS IT OTHER STACKED 3D FIELD
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*Nr*NDV3D) THEN
          MSnz=Nr
          MSnf=NDV3D
CMM: loop_g terms
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*Nr*2) THEN
          MSnz=Nr
          MSnf=2
CMM: IS IT A 2D FIELD
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy) THEN
CMM: all default values
CMM: IS IT A STACKED 2D FIELD
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*NDV2D) THEN
           MSnf=NDV2D
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*NEXF1) THEN
           MSnf=NEXF1
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*NEXF2) THEN
           MSnf=NEXF2
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*NCTRL1) THEN
           MSnf=NCTRL1
catn add loop for sea ice
#ifdef ALLOW_SEAICE
       ELSEIF (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*NSI) THEN
           MSnf=NSI
#endif
#if (defined (ALLOW_GENARR2D_CONTROL) || defined (ALLOW_GENARR3D_CONTROL) || defined (ALLOW_GENTIM2D_CONTROL))
       ELSEIF 
     &    (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*maxCtrlTim2D) THEN
           MSnf=maxCtrlTim2D
#endif
       ELSE
          print*,'CAV) issue with input size in mdsio_write_tape ???'
          print*,'CAV) file ', fName
          print*,'CAV) has size ', nSize
        WRITE(msgBuf,'(A,I10)')
     &   ' MDS_WRITE_TAPE: argument nsize =', nSize
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)')
     &   ' MDS_WRITE_TAPE: invalid value for nSize'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
  
       ENDIF
      ENDIF

catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
CMM. GOING TO WRITE BY RECORD EITHER 2D or 3D or OBC FILES.
      IF (IsOBCN .OR. IsOBCS) THEN
        MSvec_size = nxb*MSnz
      ELSEIF (IsOBCE .OR. IsOBCW) THEN
        MSvec_size = nyb*MSnz
#else
      IF (.FALSE.) THEN
#endif
      ELSE
        MSvec_size = nxb*nyb*MSnz
      ENDIF
CMM. SO NEED TO CONVERT IREC 
      mmirec = (iRec-1)*MSnf

      IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, opens a global file

         WRITE(dataFName,'(2a)') fName(1:IL),'.data'
         length_of_rec = MDS_RECLEN( filePrec, MSvec_size, myThid )
         IF (iRec .EQ. 1) THEN
          OPEN( dUnit, file=dataFName, status=_NEW_STATUS,
     &          access='direct', recl=length_of_rec )
         ELSE
          OPEN( dUnit, file=dataFName, status=_OLD_STATUS,
     &          access='direct', recl=length_of_rec )
         ENDIF
      ENDIF

CAV inspired from MITgcm/eesupp/src/global_sum_singlecpu.F
C--   copy (and conversion to real*8) to Shared buffer:
CMM Z BINNING NOT YET IMPLEMENTED
CMM IT IS NOT OBC  
      DO MSaf=1,MSnf
       DO MSak=1,MSnz
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO MSaj=1,sNy
           DO MSai=1,sNx
CMM FLAG IF OUT OF OBC LIMITS SO NOT INCLUDED
            Ifill = .FALSE.
            IF (hFacC(MSai,MSaj,MSak,bi,bj).GT.0)  THEN
             tOLx = OLx
             tOLy = OLy
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
             IF (IsOBCN .AND. (MSaj.EQ.sNy)) THEN
              Ifill = .TRUE.
              tOLy = 1-MSaj
             ELSEIF (IsOBCS.AND.(MSaj.EQ.1)) THEN
              Ifill = .TRUE.
              tOLy = 1-MSaj
             ELSEIF (IsOBCE.AND.(MSai.EQ.sNx)) THEN
              Ifill = .TRUE.
              tOLx = 1-MSai
             ELSEIF (IsOBCW.AND.(MSai.EQ.1)) THEN
              Ifill = .TRUE.
              tOLx = 1-MSai
             ELSEIF (IsOBCN.OR.IsOBCS.OR.IsOBCE.OR.IsOBCW) THEN
              Ifill = .FALSE.
             ELSE
#else
             IF (.TRUE.) THEN
#endif
              Ifill = .TRUE.
             ENDIF
            ENDIF
CMM ALSO FLAG IF LAND
            IF (Ifill) THEN
                 tmpind1 = (MSaf-1)*nxt*nyt*MSnz*nSx*nSy
     &           + (bj-1)*nxt*nyt*MSnz*nSx
     &           + (bi-1)*nxt*nyt*MSnz
     &           + (MSak-1)*nxt*nyt
     &           + (MSaj+tOLy-1)*nxt
     &           + (MSai+tOLx) 

              IF ( filePrec.EQ.precFloat32 ) THEN
                sharedLocBuf_r4(MSai,MSaj,bi,bj) = local_r4(tmpind1)
              ELSEIF ( filePrec.EQ.precFloat64 ) THEN
                sharedLocBuf_r8(MSai,MSaj,bi,bj) = local_r8(tmpind1)
              ENDIF
            ELSE
                 sharedLocBuf_r4(MSai,MSaj,bi,bj) = -9999.9
                 sharedLocBuf_r8(MSai,MSaj,bi,bj) = -9999.9
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C--   Master thread does the communications and the global sum
C--   Master thread cannot start until everyone is ready:
        CALL BAR2( myThid )

C--   Gather local arrays
        IF ( filePrec.EQ.precFloat32 ) THEN
          CALL GATHER_2D_R4(
     O                       xy_buffer_r4,
     I                       sharedLocBuf_r4,
     I                       xSize, ySize,
     I                       useExch2ioLayOut, zeroBuff, myThid )
        ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          CALL GATHER_2D_R8(
     O                       xy_buffer_r8,
     I                       sharedLocBuf_r8,
     I                       xSize, ySize,
     I                       useExch2ioLayOut, zeroBuff, myThid )
        ENDIF

C All threads wait for Master to finish 
      CALL BAR2( myThid )


CAV If anisotropic coupling: downscale before writing
CAV by bin averaging each field from global buffer

CAV for each point of the bin-averaged array
        DO MSai=1,nxb
         DO MSaj=1,nyb

          tmp_bin_r4 = 0.0
          tmp_bin_r8 = 0.0
          nbin = 0

CAV for each point of the original vector considered in the average
          DO MSpx=1,MSfacX
           DO MSpy=1,MSfacY

        IF ( filePrec.EQ.precFloat32 ) THEN
            IF (xy_buffer_r4(
     &                        ((MSaj-1)*MSfacY+MSpy-1)*xSize+
     &                         (MSai-1)*MSfacX+MSpx).gt.-9999.0) THEN
             tmp_bin_r4 = tmp_bin_r4 + xy_buffer_r4(
     &                        ((MSaj-1)*MSfacY+MSpy-1)*xSize+
     &                         (MSai-1)*MSfacX+MSpx)
             nbin = nbin+1
            ENDIF
        ELSEIF ( filePrec.EQ.precFloat64 ) THEN
            IF (xy_buffer_r8(
     &                        ((MSaj-1)*MSfacY+MSpy-1)*xSize+
     &                         (MSai-1)*MSfacX+MSpx).gt.-9999.0) THEN
             tmp_bin_r8 = tmp_bin_r8 + xy_buffer_r8(
     &                        ((MSaj-1)*MSfacY+MSpy-1)*xSize+ 
     &                         (MSai-1)*MSfacX+MSpx)
             nbin = nbin+1
            ENDIF
        ENDIF
           ENDDO
          ENDDO

           tmpind1 = nxb*nyb*(MSak-1)+nxb*(MSaj-1)+MSai

           IF (nbin.gt.0) THEN
            IF (MSnz.eq.1) THEN
             fld_2d_r4(tmpind1) = tmp_bin_r4/nbin
             fld_2d_r8(tmpind1) = tmp_bin_r8/nbin
            ELSEIF (MSnz.eq.Nr) THEN
             fld_3d_r4(tmpind1) = tmp_bin_r4/nbin
             fld_3d_r8(tmpind1) = tmp_bin_r8/nbin
            ENDIF
           ELSE
            IF (MSnz.eq.1) THEN
             fld_2d_r4(tmpind1) = 0
             fld_2d_r8(tmpind1) = 0
            ELSEIF (MSnz.eq.Nr) THEN
             fld_3d_r4(tmpind1) = 0
             fld_3d_r8(tmpind1) = 0
            ENDIF
           ENDIF

         ENDDO
        ENDDO

CMM IF OBCS CROP WHAT YOU NEED
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
        IF (IsOBCN .OR. IsOBCS) THEN
          IF (IsOBCN) MSaj = nyb
          IF (IsOBCS) MSaj = 1
          DO MSai=1,xSize
           tmpind1 = nxb*nyb*(MSak-1)+nxb*(MSaj-1)+MSai
           tmpind2 = nxb*(MSak-1)+MSai
           IF (MSnz.eq.1) THEN
             fld_obx2_r4(tmpind2) = fld_2d_r4(tmpind1)
             fld_obx2_r8(tmpind2) = fld_2d_r8(tmpind1)
           ELSEIF (MSnz.eq.Nr) THEN
             fld_obx3_r4(tmpind2) = fld_3d_r4(tmpind1)
             fld_obx3_r8(tmpind2) = fld_3d_r8(tmpind1)
           ENDIF
          ENDDO
        ENDIF
        IF (IsOBCE .OR. IsOBCW) THEN
          IF (IsOBCE) MSai = nxb
          IF (IsOBCW) MSai = 1
          DO MSaj=1,nyb
           tmpind1 = nxb*nyb*(MSak-1)+nxb*(MSaj-1)+MSai
           tmpind2 = nyb*(MSak-1)+MSaj
           IF (MSnz.eq.1) THEN
             fld_obx2_r4(tmpind2) = fld_2d_r4(tmpind1)
             fld_obx2_r8(tmpind2) = fld_2d_r8(tmpind1)
           ELSEIF (MSnz.eq.Nr) THEN
             fld_obx3_r4(tmpind2) = fld_3d_r4(tmpind1)
             fld_obx3_r8(tmpind2) = fld_3d_r8(tmpind1)
           ENDIF
          ENDDO
        ENDIF
#endif

       ENDDO /* k loop */

      IF ( myProcId .EQ. 0 ) THEN
      mmirec = mmirec + 1

C-    Write global buffer to file:
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
        IF (IsOBCN .OR. IsOBCS) THEN
         IF ( filePrec.EQ.precFloat32 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_obx2_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_obx2_r4
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_obx3_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_obx3_r4
          ENDIF
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_obx2_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_obx2_r8
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_obx3_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_obx3_r8
          ENDIF
         ENDIF
        ELSEIF (IsOBCE .OR. IsOBCW) THEN
         IF ( filePrec.EQ.precFloat32 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_oby2_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_oby2_r4
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_oby3_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_oby3_r4
          ENDIF
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_oby2_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_oby2_r8
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_oby3_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_oby3_r8
          ENDIF
         ENDIF
#else
        IF (.FALSE.) THEN
#endif
        ELSE
         IF ( filePrec.EQ.precFloat32 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_2d_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_2d_r4
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR4( MSvec_size, fld_3d_r4 )
#endif
           WRITE(dUnit,rec=mmirec) fld_3d_r4
          ENDIF
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          IF (MSnz.eq.1) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_2d_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_2d_r8
          ELSEIF (MSnz.eq.Nr) THEN
#ifdef _BYTESWAPIO
          CALL MDS_BYTESWAPR8( MSvec_size, fld_3d_r8 )
#endif
           WRITE(dUnit,rec=mmirec) fld_3d_r8
          ENDIF
         ENDIF
        ENDIF
      ENDIF
      ENDDO /* MSnf loop */

      IF ( myProcId .EQ. 0 ) THEN
C-    Close data-file and create meta-file
         CLOSE( dUnit )

         IL  = ILNBLNK( fName )
         WRITE(metaFName,'(2A)') fName(1:IL),'.meta'
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
        IF (IsOBCE .OR. IsOBCW) THEN
         dimList(1,1) = 1
         dimList(2,1) = 1
         dimList(3,1) = 1
#else
        IF (.FALSE.) THEN
#endif
        ELSE
         dimList(1,1) = nxb
         dimList(2,1) = 1
         dimList(3,1) = nxb
        ENDIF
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
        IF (IsOBCN .OR. IsOBCS) THEN
         dimList(1,2) = 1
         dimList(2,2) = 1
         dimList(3,2) = 1
#else
        IF (.FALSE.) THEN
#endif
        ELSE
         dimList(1,2) = nyb
         dimList(2,2) = 1
         dimList(3,2) = nyb
        ENDIF
         dimList(1,3) = MSnz
         dimList(2,3) = 1
         dimList(3,3) = MSnz
         nDims=3
         IF (MSnz.EQ.1) nDims=2
         map2gl(1) = 0
         map2gl(2) = 1
         CALL MDS_WRITE_META(
     I              metaFName, dataFName, the_run_name, ' ',
     I              filePrec, nDims, dimList, map2gl, MSnf, blank8c,
     I              0, dummyRL, oneRL, iRec, myIter, myThid )

C-     end if myProcId=0
       ENDIF

CAV   else if singleCpuIO
      ELSE

        print*,'CAV) MULTISCALING only works with singleCpuIO'
        STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'

CAV   end if singleCpuIO
      ENDIF 

#else /* MULTISCALE_COUPLING_TAPES */
CAV)

C     If option globalFile is desired but does not work or if
C     globalFile is too slow, then try using single-CPU I/O.
      IF ( singleCpuIO ) THEN

C-    Gather array from all procs
        IF ( filePrec.EQ.precFloat32 ) THEN
          CALL GATHER_VEC_R4( gl_buffer_r4, local_r4, nSize, myThid )
        ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          CALL GATHER_VEC_R8( gl_buffer_r8, local_r8, nSize, myThid )
        ELSE
          WRITE(msgBuf,'(A)')
     &            ' MDS_WRITE_TAPE: illegal value for filePrec'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
        ENDIF

        IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, opens a global file

         WRITE(dataFName,'(2a)') fName(1:IL),'.data'
         length_of_rec = MDS_RECLEN( filePrec, vec_size, myThid )
         IF (iRec .EQ. 1) THEN
          OPEN( dUnit, file=dataFName, status=_NEW_STATUS,
     &          access='direct', recl=length_of_rec )
         ELSE
          OPEN( dUnit, file=dataFName, status=_OLD_STATUS,
     &          access='direct', recl=length_of_rec )
         ENDIF

C-    Write global buffer to file:
         IF ( filePrec.EQ.precFloat32 ) THEN
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR4( vec_size, gl_buffer_r4 )
#endif
           WRITE(dUnit,rec=iRec) gl_buffer_r4
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR8( vec_size, gl_buffer_r8 )
#endif
           WRITE(dUnit,rec=iRec) gl_buffer_r8
         ENDIF

C-    Close data-file and create meta-file
         CLOSE( dUnit )
         WRITE(metaFName,'(2a)') fName(1:IL),'.meta'
         dimList(1,1) = vec_size
         dimList(2,1) = 1
         dimList(3,1) = vec_size
         nDims = 1
         CALL MDS_WRITE_META(
     I              metaFName, dataFName, the_run_name, ' ',
     I              filePrec, nDims, dimList, map2gl, 0, blank8c,
     I              0, dummyRL, oneRL, iRec, myIter, myThid )

C-    end if myProcId=0
        ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     if ( singleCpuIO ), else
      ELSEIF ( .NOT. singleCpuIO ) THEN

        IF ( globalFile ) THEN
C-    If we are writing to a global file then we open it here
         WRITE(dataFName,'(2A)') fName(1:IL),'.data'
         length_of_rec = MDS_RECLEN( filePrec, nSize, myThid )
         IF ( iRec.EQ.1 ) THEN
          OPEN( dUnit, file=dataFName, status=_NEW_STATUS,
     &          access='direct', recl=length_of_rec )
         ELSE
          OPEN( dUnit, file=dataFName, status=_OLD_STATUS,
     &          access='direct', recl=length_of_rec )
         ENDIF

        ELSE
C-    If we are writing to a tiled MDS file then we open each one here
         iG = 1 + (myXGlobalLo-1)/sNx
         jG = 1 + (myYGlobalLo-1)/sNy
         WRITE(dataFName,'(2A,I3.3,A,I3.3,A)')
     &              pfName(1:pIL),'.',iG,'.',jG,'.data'
         length_of_rec = MDS_RECLEN( filePrec, nSize, myThid )
         IF (iRec .EQ. 1) THEN
          OPEN( dUnit, file=dataFName, status=_NEW_STATUS,
     &          access='direct', recl=length_of_rec )
         ELSE
          OPEN( dUnit, file=dataFName, status=_OLD_STATUS,
     &          access='direct', recl=length_of_rec )
         ENDIF
        ENDIF

C-    Write local buffer to file:
        IF (globalFile) THEN
C-- Original: nPy=2, nSx=2 -> produces too large file (1.5 x normal size)
c          iG   = myXGlobalLo-1+(bi-1)*sNx
c          jG   = myYGlobalLo-1+(bj-1)*sNy
c          jRec = 1 + int(iG/sNx) + (jG/sNy)*nSx*nPx +
c    &            (iRec-1)*nSx*nPx*nSy*nPy
C-- Alternative: same layout as in scatter/gather_vector (for singleCpuIO)
C   problem: nPx=2, nSx=2, writing a global (i.e., with bi,bj dim);
C-         2nd proc get iG=3 -> badly placed data over nPx*nPy*nSize range
C                               that will be overwritten by next record
c          iG   = 1 + (myXGlobalLo-1)/sNx
c          jG   = 1 + (myYGlobalLo-1)/sNy
c          jRec = iG + (jG-1)*nPx + (iRec-1)*nPx*nPy
C-- Simpler: should work (but hard to interpret the sequence of data in file)
           jRec = 1 + myProcId + (iRec-1)*nPx*nPy
        ELSE
           jRec = iRec
        ENDIF
        IF ( filePrec.EQ.precFloat32 ) THEN
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR4( nSize, local_r4 )
#endif
           WRITE(dUnit,rec=jRec) local_r4
        ELSEIF ( filePrec.EQ.precFloat64 ) THEN
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR8( nSize, local_r8 )
#endif
           WRITE(dUnit,rec=jRec) local_r8
        ELSE
           WRITE(msgBuf,'(A)')
     &           ' MDS_WRITE_TAPE: illegal value for filePrec'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE'
        ENDIF

C-    Close data-file and create meta-file
        CLOSE( dUnit )
        IF ( globalFile ) THEN
C     meta-file for global file
          WRITE(metaFName,'(2A)') fName(1:IL),'.meta'
          dimList(1,1) = vec_size
          dimList(2,1) = 1
          dimList(3,1) = vec_size
          nDims = 1
        ELSE
C     meta-file for tiled file
          iG = 1 + (myXGlobalLo-1)/sNx
          jG = 1 + (myYGlobalLo-1)/sNy
          WRITE(metaFName,'(2A,I3.3,A,I3.3,A)')
     &             pfName(1:pIL),'.',iG,'.',jG,'.meta'
          dimList(1,1) = nPx*nPy*nSize
          dimList(2,1) = 1 + myProcId*nSize
          dimList(3,1) = (1+myProcId)*nSize
          nDims = 1
        ENDIF
C-    write meta-file
        CALL MDS_WRITE_META(
     I              metaFName, dataFName, the_run_name, ' ',
     I              filePrec, nDims, dimList, map2gl, 0, blank8c,
     I              0, dummyRL, oneRL, iRec, myIter, myThid )
c    I              metaFName, dataFName, the_run_name, titleLine,
c    I              filePrec, nDims, dimList, map2gl, nFlds, fldList,
c    I              nTimRec, timList, misVal, iRec, myIter, myThid )

C     end-if ( .not. singleCpuIO )
      ENDIF

CAV(
#endif /* MULTISCALE_COUPLING_TAPES */
CAV)

      _END_MASTER( myThid )

#else /* ALLOW_AUTODIFF */
      STOP 'ABNORMAL END: S/R MDS_WRITE_TAPE is empty'
#endif /* ALLOW_AUTODIFF */

      RETURN
      END
