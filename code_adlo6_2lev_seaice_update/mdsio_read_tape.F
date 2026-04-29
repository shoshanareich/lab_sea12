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
C !ROUTINE: MDS_READ_TAPE
C !INTERFACE:
      SUBROUTINE MDS_READ_TAPE(
     I   fName,
     I   filePrec,
     I   arrType,
     I   nSize,
     O   fldR8, fldR4,
     I   singleCpuIO,
     I   iRec,
     I   myThid )

C !DESCRIPTION:
C MDS_READ_TAPE: load an array (treated as vector) for a tape-file
C  (renamed from MDSREADVECTOR with 2 explicit output array typest)
C
C Arguments:
C fName      string  :: base name for file to read
C filePrec   integer :: number of bits per word in file (32 or 64)
C arrType    char(2) :: which array (fldR8/R4) to read, either "R8" or "R4"
C nSize      integer :: number of elements of input array "fldR8/R4" to read
C fldR8      ( R8 )  :: array to read if arrType="R8", fldR8(nSize)
C fldR4      ( R4 )  :: array to read if arrType="R4", fldR4(nSize)
C singleCpuIO ( L )  :: only proc 0 do IO and send data to other procs
C iRec       integer :: record number to read
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
      CHARACTER*(2) arrType
      INTEGER nSize
      _R8     fldR8(*)
      _R4     fldR4(*)
      LOGICAL singleCpuIO
      INTEGER iRec
      INTEGER myThid

#ifdef ALLOW_AUTODIFF

C !FUNCTIONS:
      INTEGER ILNBLNK
      INTEGER MDS_RECLEN
      EXTERNAL ILNBLNK
      EXTERNAL MDS_RECLEN

C !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_FNAM) dataFName, pfName
      INTEGER iG, jG, jRec, dUnit, IL, pIL
      LOGICAL exst
      LOGICAL globalFile, fileIsOpen
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
CAV(
#ifdef MULTISCALE_COUPLING_TAPES
      INTEGER   nxt,nyt,tOLx,tOLy,mmirec
      INTEGER   fld_size,tmpind1,tmpind2
      INTEGER   bi,bj,MSai,MSaj,MSak,MSaf,MSnz,MSnf
#ifdef ALLOW_EXCH2
      _R8    fld_2d_r8(exch2_global_Nx*exch2_global_NY)
      _R4    fld_2d_r4(exch2_global_Nx*exch2_global_NY)
      _R8    fld_3d_r8(exch2_global_Nx*exch2_global_NY*Nr)
      _R4    fld_3d_r4(exch2_global_Nx*exch2_global_NY*Nr)
      _R8    fld_obx2_r8(exch2_global_Nx)
      _R4    fld_obx2_r4(exch2_global_Nx)
      _R8    fld_obx3_r8(exch2_global_Nx*Nr)
      _R4    fld_obx3_r4(exch2_global_Nx*Nr)
      _R8    fld_oby2_r8(exch2_global_Ny)
      _R4    fld_oby2_r4(exch2_global_Ny)
      _R8    fld_oby3_r8(exch2_global_Ny*Nr)
      _R4    fld_oby3_r4(exch2_global_Ny*Nr)
#else
      _R8    fld_2d_r8(NX*NY)
      _R4    fld_2d_r4(NX*NY)
      _R8    fld_3d_r8(NX*NY*Nr)
      _R4    fld_3d_r4(NX*NY*Nr)
      _R8    fld_obx2_r8(Nx)
      _R4    fld_obx2_r4(Nx)
      _R8    fld_obx3_r8(Nx*Nr)
      _R4    fld_obx3_r4(Nx*Nr)
      _R8    fld_oby2_r8(Ny)
      _R4    fld_oby2_r4(Ny)
      _R8    fld_oby3_r8(Ny*Nr)
      _R4    fld_oby3_r4(Ny*Nr)
#endif
      LOGICAL useExch2ioLayOut
      INTEGER xSize, ySize
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      LOGICAL   IsOBCN, IsOBCE,IsOBCS,IsOBCW
#endif
      LOGICAL   Ifill
#endif
CAV)

CEOP

      vec_size = nSize*nPx*nPy

C--   Only do I/O if I am the master thread
      _BEGIN_MASTER( myThid )

C-    Record number must be >= 1
      IF ( iRec.LT.1 ) THEN
       WRITE(msgBuf,'(A,I10)')
     &   ' MDS_READ_TAPE: argument iRec =',iRec
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' MDS_READ_TAPE: invalid value for iRec'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
      ENDIF

C-    Assume nothing
      globalFile = .FALSE.
      fileIsOpen = .FALSE.
      IL  = ILNBLNK( fName )
      pIL = ILNBLNK( mdsioLocalDir )

CAV(
#ifdef MULTISCALE_COUPLING_TAPES
      useExch2ioLayOut = .FALSE.
      xSize = Nx
      ySize = Ny
      Ifill = .FALSE.
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
      IsOBCN = .FALSE.
      IsOBCE = .FALSE.
      IsOBCS = .FALSE.
      IsOBCW = .FALSE.
#endif
#ifdef ALLOW_EXCH2
      useExch2ioLayOut = .TRUE.
      xSize = exch2_global_Nx
      ySize = exch2_global_Ny
#endif 
#endif 
CAV)

C-    Assign special directory
      IF ( mdsioLocalDir .NE. ' ' ) THEN
       WRITE(pfName,'(2A)') mdsioLocalDir(1:pIL), fName(1:IL)
      ELSE
       pfName = fName
      ENDIF
      pIL = ILNBLNK( pfName )

C-    Assign a free unit number as the I/O channel for this routine
      CALL MDSFINDUNIT( dUnit, myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      IF ( singleCpuIO ) THEN

cc        IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, opens a global file

CAV If using coarser adjoint, need to read tapes in the order written
CAV in mdsio_write_tape
#ifdef MULTISCALE_COUPLING_TAPES

CMM INITIALIZE LOCAL BUFFER
      DO j=1,nSize
        local_r4(j) = 0.0
        local_r8(j) = 0.0
      ENDDO
CMM ALSO fld_?d_r?
      DO j=1,xSize*ySize
        fld_2d_r4(j) = 0.0
        fld_2d_r8(j) = 0.0
      ENDDO
      DO j=1,xSize*ySize*NR
        fld_3d_r4(j) = 0.0
        fld_3d_r8(j) = 0.0
      ENDDO

CMM DEBUG
cmm        IF (myProcId.EQ.44) print*,'CMMa: ',pfName(pIL-20:pIL),nSize

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
cmm      print*,'CMMb: pIL, fname ',pIL,pfName(pIL-7-7:pIL),IsOBCN
CMM)

CMM(
CAV dimensions of each tile, including overlap
CMM: SET TO DEFAULT SIZES. WILL OVERWRITE FOR UNIQUE CASES
      MSnz = 1
      MSnf = 1
      nxt = sNx+2*OLx
      nyt = sNy+2*OLy

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
     &     (nSize.eq.(sNx+2*OLx)*(sNy+2*OLy)*nSx*nSy*maxCtrlTim2D) THEN
          MSnf=maxCtrlTim2D
#endif
       ELSE
          print*,'CAV) issue with input size in mdsio_write_tape ???'
#ifdef ALLOW_PTRACERS
          print*,'mdsio_write_tape / nsize, ptracers_num',
     &           nSize,ptracers_num
#endif
       WRITE(msgBuf,'(A,I10)')
     &   ' MDS_READ_TAPE: argument nsize =', nSize
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' MDS_READ_TAPE: invalid value for nSize'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R MDS_READ_TAPE'

       ENDIF
      ENDIF

catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
CMM. GOING TO READ BY RECORD EITHER 2D or 3D or OBC FILES.
      IF (IsOBCN .OR. IsOBCS) THEN
        fld_size = xSize*MSnz
      ELSEIF (IsOBCE .OR. IsOBCW) THEN
        fld_size = ySize*MSnz
#else
      IF (.FALSE.) THEN
#endif
      ELSE
        fld_size = xSize*ySize*MSnz
      ENDIF
CMM. SO NEED TO CONVERT IREC = size nsize*npx*npy
         mmirec = (iRec-1)*MSnf

        IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, opens a global file

cmm       IF (IsOBCN.EQ..TRUE.)
cmm     & print*,'CMMc: mmirec myProcId ',iRec, MSnf, mmirec,myProcId

C-    Check first for global file with with MDS name (ie. fName.data)
         WRITE(dataFName,'(2A)') fName(1:IL),'.data'
         INQUIRE( file=dataFName, exist=exst )
         IF (exst) globalFile = .TRUE.

C-    If global file is visible to process 0, then open it here.
         IF ( globalFile ) THEN

          IF ( debugLevel .GE. debLevB ) THEN
           WRITE(msgBuf,'(A,A)')
     &     ' MDS_READ_TAPE: opening global file: ',dataFName(1:IL+5)
           CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
          length_of_rec = MDS_RECLEN( filePrec, fld_size, myThid )

          OPEN( dUnit, file=dataFName, status='old',
     &          access='direct', recl=length_of_rec )
         ELSE
C     Otherwise stop program.
          WRITE(msgBuf,'(2A)')
     &      ' MDS_READ_TAPE: filename: ',dataFName(1:IL)
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(A)')
     &      ' MDS_READ_TAPE: File does not exist'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
         ENDIF

         ENDIF /* ( myProcId .EQ. 0 ) */

C-    Read into global buffer:

      DO MSaf=1,MSnf

        mmirec = mmirec + 1

        IF ( myProcId .EQ. 0 ) THEN

cmm       IF (IsOBCN)
cmm     & print*,'CMMd: MSaf, mmirec, MSnz ',MSaf, mmirec,fld_size,MSnz,filePrec

catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
C--   Master thread of process 0, only, reads a global file
         IF (IsOBCN .OR. IsOBCS) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_obx2_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_obx2_r4 )
#endif
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_obx3_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_obx3_r4 )
#endif
           ENDIF
          ELSEIF ( filePrec.EQ.precFloat64 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_obx2_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_obx2_r8 )
#endif
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_obx3_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_obx3_r8 )
#endif
           ENDIF
          ENDIF
         ELSEIF (IsOBCE .OR. IsOBCW) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_oby2_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_oby2_r4 )
#endif
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_oby3_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_oby3_r4 )
#endif
           ENDIF
          ELSEIF ( filePrec.EQ.precFloat64 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_oby2_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_oby2_r8 )
#endif
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_oby3_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_oby3_r8 )
#endif
           ENDIF
          ENDIF
#else
         IF (.FALSE.) THEN
#endif
         ELSE
          IF ( filePrec.EQ.precFloat32 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_2d_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_2d_r4 )
#endif
cmm      print*,'CMM RdChk2: ',fld_2d_r4(216500),fld_2d_r4(216650)
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_3d_r4
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_3d_r4 )
#endif
cmm      print*,'CMM RdChk3: ',fld_3d_r4(216500),fld_3d_r4(216650)
           ENDIF
          ELSEIF ( filePrec.EQ.precFloat64 ) THEN
           IF (MSnz.eq.1) THEN
            READ(dUnit,rec=mmirec) fld_2d_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_2d_r8 )
#endif
           ELSEIF (MSnz.eq.Nr) THEN
            READ(dUnit,rec=mmirec) fld_3d_r8
#ifdef _BYTESWAPIO
            CALL MDS_BYTESWAPR4( fld_size, fld_3d_r8 )
#endif
           ENDIF
          ENDIF
         ENDIF

cmm       IF (IsOBCN)
cmm     & print*,'CMMe: MSaf, mmirec, fld_size, MSnz ',MSaf, mmirec, fld_size,MSnz

        ENDIF /* ( myProcId .EQ. 0 ) */


CAV must reorder according to what adjoint would have written 
CAV given SIZE.h in low res coMSnfiguration
CMM. IF OBCS REPOPULATE fld_?d_r?
CMM(
       DO MSak=1,MSnz
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
        IF (IsOBCN .OR. IsOBCS) THEN
          IF (IsOBCN) MSaj = ysize
          IF (IsOBCS) MSaj = 1
          DO MSai=1,xSize
           tmpind1 = xSize*ySize*(MSak-1)+xSize*(MSaj-1)+MSai
           tmpind2 = xSize*(MSak-1)+MSai
           IF (MSnz.eq.1) THEN
             fld_2d_r4(tmpind1) = fld_obx2_r4(tmpind2)
             fld_2d_r8(tmpind1) = fld_obx2_r8(tmpind2)
           ELSEIF (MSnz.eq.Nr) THEN
             fld_3d_r4(tmpind1) = fld_obx3_r4(tmpind2)
             fld_3d_r8(tmpind1) = fld_obx3_r8(tmpind2)
           ENDIF
          ENDDO
        ENDIF
        IF (IsOBCE .OR. IsOBCW) THEN
          IF (IsOBCE) MSai = xSize
          IF (IsOBCW) MSai = 1
          DO MSaj=1,ySize
           tmpind1 = xSize*ySize*(MSak-1)+xSize*(MSaj-1)+MSai
           tmpind2 = ySize*(MSak-1)+MSaj
           IF (MSnz.eq.1) THEN
             fld_2d_r4(tmpind1) = fld_obx2_r4(tmpind2)
             fld_2d_r8(tmpind1) = fld_obx2_r8(tmpind2)
           ELSEIF (MSnz.eq.Nr) THEN
             fld_3d_r4(tmpind1) = fld_obx3_r4(tmpind2)
             fld_3d_r8(tmpind1) = fld_obx3_r8(tmpind2)
           ENDIF
          ENDDO
        ENDIF
#endif

            DO MSaj=1,ySize
             DO MSai=1,xSize
         IF ( filePrec.EQ.precFloat32 ) THEN
              IF (MSnz.eq.1) THEN
               xy_buffer_r4(MSai+(MSaj-1)*xSize)
     &           = fld_2d_r4(MSai+(MSaj-1)*xSize)
              ELSEIF (MSnz.eq.Nr) THEN
               xy_buffer_r4(MSai+(MSaj-1)*xSize)
     &           = fld_3d_r4(MSai+(MSaj-1)*xSize + (MSak-1)*xSize*ySize)
              ENDIF
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
              IF (MSnz.eq.1) THEN
               xy_buffer_r8(MSai+(MSaj-1)*xSize)
     &           = fld_2d_r8(MSai+(MSaj-1)*xSize)
              ELSEIF (MSnz.eq.Nr) THEN
               xy_buffer_r8(MSai+(MSaj-1)*xSize)
     &          = fld_3d_r8(MSai+(MSaj-1)*xSize + (MSak-1)*xSize*ySize)
              ENDIF
         ENDIF
             ENDDO
            ENDDO

cmm              IF (MSak.eq.1) THEN
cmm      print*,'CMM RdChk3: ',xy_buffer_r4(216500),xy_buffer_r4(216650)
cmm              ENDIF

cmm       IF (IsOBCN)
cmm     & print*,'CMMg: MSaj, MSak, MSaf, tmpind1, tmpind2 ', MSaj, MSak, MSaf, xSize, ySize

C All threads wait for Master to finish
            CALL BAR2( myThid )
         IF ( filePrec.EQ.precFloat32 ) THEN
            CALL SCATTER_2D_R4(
     U                        xy_buffer_r4,
     O                        sharedLocBuf_r4,
     I                        xSize, ySize,
     I                        useExch2ioLayOut, .FALSE., myThid )
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
            CALL SCATTER_2D_R8(
     U                        xy_buffer_r8,
     O                        sharedLocBuf_r8,
     I                        xSize, ySize,
     I                        useExch2ioLayOut, .FALSE., myThid )
         ENDIF             
C All threads wait for Master to finish 
            CALL BAR2( myThid )

            DO bj=1,nSy
             DO bi=1,nSx
              DO MSaj=1,sNy
               DO MSai=1,sNx
                Ifill = .FALSE.
                IF (hFacC(MSai,MSaj,MSak,bi,bj).GT.0)  THEN
                 tOLX = OLx
                 tOLy = OLy
catn add loop checking for obcs before using nob
#if ( defined ALLOW_OBCS && defined AUTODIFF_USE_STORE_RESTORE_OBCS )
                 IF (IsOBCN .AND. (MSaj.EQ.sNy)) THEN
                  Ifill = .TRUE.
                  tOLy = 1-MSaj
                 ELSEIF (IsOBCS .AND. (MSaj.EQ.1)) THEN 
                  Ifill = .TRUE.
                  tOLy = 1-MSaj
                 ELSEIF (IsOBCE .AND. (MSai.EQ.sNx)) THEN
                  Ifill = .TRUE.
                  tOLx = 1-MSai
                 ELSEIF (IsOBCW .AND. (MSai.EQ.1) )THEN
                  Ifill = .TRUE.
                  tOLx = 1-MSai
                 ELSEIF (IsOBCN.OR.IsOBCS.OR.IsOBCE.OR.IsOBCW) THEN
                  Ifill = .FALSE.
#else
                 IF (.FALSE.) THEN
#endif
                 ELSE
                  Ifill = .TRUE.
                 ENDIF
                ENDIF
                IF (Ifill) THEN
                 tmpind1 = (MSai+tOLx)
     &            + (MSaj+tOLy-1)*nxt
     &            + (MSak-1)*nxt*nyt
     &            + (bi-1)*nxt*nyt*MSnz
     &            + (bj-1)*nxt*nyt*MSnz*nSx
     &            + (MSaf-1)*nxt*nyt*MSnz*nSx*nSy

cmm       IF (IsOBCN.EQ..TRUE.)
cmm     & print*,'CMMi: sNy, nxt,nyt,tmpind1 ', sNy, nxt,nyt,MSnz,tmpind1,tOLx,tOLy

                 IF ( filePrec.EQ.precFloat32 ) THEN

                  local_r4(tmpind1) = sharedLocBuf_r4(MSai,MSaj,bi,bj)

CMM: DEBUG
cmm        IF ((local_r4(tmpind1).NE.0).AND.(myProcId.EQ.44)
cmm     &        .AND.(MSai*MSaj*MSak*MSaf.EQ.1) )
cmm     &       print*,'CMM RdChk: ',nSize,MSnz,MSnf,irec,local_r4(tmpind1)

                 ELSEIF ( filePrec.EQ.precFloat64 ) THEN
                   local_r8(tmpind1) = sharedLocBuf_r8(MSai,MSaj,bi,bj)
                 ENDIF
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDDO

CMM( AK LOOP 
           ENDDO
CMM( AF LOOP 
          ENDDO
 
        IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, closes a global file
C-    Close data-file
         CLOSE( dUnit )

        ENDIF /* ( myProcId .EQ. 0 ) */



#else /* if not defined MULTISCALE_COUPLING_TAPES */

        IF ( myProcId .EQ. 0 ) THEN
C--   Master thread of process 0, only, opens a global file

C-    Check first for global file with with MDS name (ie. fName.data)
         WRITE(dataFName,'(2A)') fName(1:IL),'.data'
         INQUIRE( file=dataFName, exist=exst )
         IF (exst) globalFile = .TRUE.

C-    If global file is visible to process 0, then open it here.
         IF ( globalFile ) THEN
          IF ( debugLevel .GE. debLevB ) THEN
           WRITE(msgBuf,'(A,A)')
     &     ' MDS_READ_TAPE: opening global file: ',dataFName(1:IL+5)
           CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
          length_of_rec = MDS_RECLEN( filePrec, vec_size, myThid )
          OPEN( dUnit, file=dataFName, status='old', _READONLY_ACTION
     &          access='direct', recl=length_of_rec )
         ELSE
C     Otherwise stop program.
          WRITE(msgBuf,'(2A)')
     &      ' MDS_READ_TAPE: filename: ',dataFName(1:IL)
C-jmc: why double print (stdout + stderr) ?
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(A)')
     &      ' MDS_READ_TAPE: File does not exist'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
         ENDIF

C-    Read into global buffer:
         IF ( filePrec.EQ.precFloat32 ) THEN
           READ(dUnit,rec=iRec) gl_buffer_r4
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR4( vec_size, gl_buffer_r4 )
#endif
         ELSEIF ( filePrec.EQ.precFloat64 ) THEN
           READ(dUnit,rec=iRec) gl_buffer_r8
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR8( vec_size, gl_buffer_r8 )
#endif
         ENDIF

C-    Close data-file
         CLOSE( dUnit )

C--   end if myProcId=0
        ENDIF

        IF ( filePrec.EQ.precFloat32 ) THEN
          CALL SCATTER_VEC_R4( gl_buffer_r4, local_r4, nSize, myThid )
        ELSEIF ( filePrec.EQ.precFloat64 ) THEN
          CALL SCATTER_VEC_R8( gl_buffer_r8, local_r8, nSize, myThid )
        ELSE
          WRITE(msgBuf,'(A)')
     &            ' MDS_READ_TAPE: illegal value for filePrec'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
        ENDIF

#endif /* if not defined MULTISCALE_COUPLING_TAPES */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     if ( singleCpuIO ), else
      ELSEIF ( .NOT. singleCpuIO ) THEN

C-    Check first for global file with with MDS name (ie. fName.data)
        WRITE(dataFName,'(2A)') fName(1:IL),'.data'
        INQUIRE( file=dataFName, exist=exst )
        IF (exst) THEN
          IF ( debugLevel .GE. debLevB ) THEN
           WRITE(msgBuf,'(A,A)')
     &     ' MDS_READ_TAPE: opening global file: ',dataFName(1:IL+5)
           CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
          globalFile = .TRUE.
C-    And open it here
          length_of_rec = MDS_RECLEN( filePrec, nSize, myThid )
          OPEN( dUnit, file=dataFName, status='old',_READONLY_ACTION 
     &          access='direct', recl=length_of_rec )
          fileIsOpen=.TRUE.
        ENDIF

C-    If we are reading from a tiled MDS file then we open each one here
        IF ( .NOT.globalFile ) THEN
         iG = 1 + (myXGlobalLo-1)/sNx
         jG = 1 + (myYGlobalLo-1)/sNy
         WRITE(dataFName,'(2A,I3.3,A,I3.3,A)')
     &              pfName(1:pIL),'.',iG,'.',jG,'.data'
         INQUIRE( file=dataFName, exist=exst )
         IF (exst) THEN
          IF ( debugLevel .GE. debLevB ) THEN
           WRITE(msgBuf,'(A,A)')
     &      ' MDS_READ_TAPE: opening file: ',dataFName(1:pIL+13)
           CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
          length_of_rec = MDS_RECLEN( filePrec, nSize, myThid )
          OPEN( dUnit, file=dataFName, status='old',_READONLY_ACTION 
     &          access='direct', recl=length_of_rec )
          fileIsOpen=.TRUE.
         ELSE
          fileIsOpen=.FALSE.
          WRITE(msgBuf,'(4A)')
     &      ' MDS_READ_TAPE: missing file: ',fName(1:IL),
     &                                 ' , ',dataFName(1:pIL+13)
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
         ENDIF
        ENDIF

        IF ( fileIsOpen ) THEN
          IF ( globalFile ) THEN
C-    read the same way it was written:
            jRec = 1 + myProcId + (iRec-1)*nPx*nPy
          ELSE
            jRec = iRec
          ENDIF
          IF ( filePrec.EQ.precFloat32 ) THEN
           READ(dUnit,rec=jRec) local_r4
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR4( nSize, local_r4 )
#endif
          ELSEIF ( filePrec.EQ.precFloat64 ) THEN
           READ(dUnit,rec=jRec) local_r8
#ifdef _BYTESWAPIO
           CALL MDS_BYTESWAPR8( nSize, local_r8 )
#endif
          ELSE
           WRITE(msgBuf,'(A)')
     &            ' MDS_READ_TAPE: illegal value for filePrec'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
          ENDIF
C--   If file was opened then close it
          CLOSE( dUnit )
          fileIsOpen = .FALSE.
        ENDIF

C     end-if ( .not. singleCpuIO )
      ENDIF

      _END_MASTER( myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Copy local buffer into output array
        IF ( arrType.EQ.'R4' ) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
            DO j=1,nSize
              fldR4(j) = local_r4(j)
            ENDDO
          ELSE
            DO j=1,nSize
              fldR4(j) = local_r8(j)
            ENDDO
          ENDIF
        ELSEIF ( arrType.EQ.'R8' ) THEN
          IF ( filePrec.EQ.precFloat32 ) THEN
            DO j=1,nSize
              fldR8(j) = local_r4(j)
            ENDDO
          ELSE
            DO j=1,nSize
              fldR8(j) = local_r8(j)
            ENDDO
          ENDIF
        ELSE
          WRITE(msgBuf,'(A)')
     &          ' MDS_READ_TAPE: illegal value for arrType'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R MDS_READ_TAPE'
        ENDIF

#else /* ALLOW_AUTODIFF */
      STOP 'ABNORMAL END: S/R MDS_READ_TAPE is empty'
#endif /* ALLOW_AUTODIFF */

      RETURN
      END
