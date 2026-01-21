#ifndef AUTODIFF_OPTIONS_H
#define AUTODIFF_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: AUTODIFF_OPTIONS.h
C !INTERFACE:
C #include "AUTODIFF_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for AutoDiff (autodiff) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_AUTODIFF
#ifdef ECCO_CPPOPTIONS_H

C-- When multi-package option-file ECCO_CPPOPTIONS.h is used (directly included
C    in CPP_OPTIONS.h), this option file is left empty since all options that
C   are specific to this package are assumed to be set in ECCO_CPPOPTIONS.h

#else /* ndef ECCO_CPPOPTIONS_H */
C   ==================================================================
C-- Package-specific Options & Macros go here

CAV when using coarser adjoint, interp model fields to adjoint grid
CAV before writing to tapes
#define MULTISCALE_COUPLING_TAPES

C o Include/exclude code in order to be able to automatically
C   differentiate the MITgcmUV by using the Tangent Linear and
C   Adjoint Model Compiler (TAMC).
#define ALLOW_AUTODIFF_TAMC

C o Checkpointing as handled by TAMC
#define ALLOW_TAMC_CHECKPOINTING
#define AUTODIFF_2_LEVEL_CHECKPOINT

C o Extract adjoint state
#define ALLOW_AUTODIFF_MONITOR

C o use divided adjoint to split adjoint computations
#define ALLOW_DIVIDED_ADJOINT
#define ALLOW_DIVIDED_ADJOINT_MPI

C o tape settings
#undef ALLOW_AUTODIFF_WHTAPEIO
#undef AUTODIFF_USE_MDSFINDUNITS
#define AUTODIFF_USE_OLDSTORE_2D
#define AUTODIFF_USE_OLDSTORE_3D
#undef EXCLUDE_WHIO_GLOBUFF_2D
#undef ALLOW_INIT_WHTAPEIO

catnC this is for Writing adj rho
catn#define ALLOW_AUTODIFF_MONITOR_DIAG

C o write separate tape files for each ptracer
#define AUTODIFF_PTRACERS_SPLIT_FILES

C o use the deprecated autodiff_store/restore method where multiple fields
C   are collected in a single buffer field array before storing to tape.
C   This functionality has been replaced by WHTAPEIO method (see above).
C   Might still be used for OBCS since WHTAPEIO does not support OBCS fields.
#define AUTODIFF_USE_STORE_RESTORE
#define AUTODIFF_USE_STORE_RESTORE_OBCS

C o allow using viscFacInAd to recompute viscosities in AD
#define AUTODIFF_ALLOW_VISCFACADJ

C   ==================================================================
#endif /* ndef ECCO_CPPOPTIONS_H */
#endif /* ALLOW_AUTODIFF */
#endif /* AUTODIFF_OPTIONS_H */
