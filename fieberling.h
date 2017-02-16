/*
** svn $Id: seamount_b.h 778 2015-08-22 01:20:04Z arango $
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Seamount Test.
**
** Application flag:   SEAMOUNT_B
** Input script:       ocean_seamount.in
*/

#define UV_ADV
#define UV_U3ADV_SPLIT
#define UV_COR
#define SOLVE3D
#define DJ_GRADPS /*Splines density Jacobian*/
#define ANA_DIAG /*diagnostics*/

/*Horizontal drag and viscosity*/
#define UV_QDRAG /*Quadratic drag coefficient (See Arbic 2008). Other options: UV_LDRAG,UV_LOGDRAG*/
!#define UV_DRAG_GRID /*spatially-varying bottom drag*/
#define UV_VIS4  /*Harmonic viscosity*/
#define MIX_S_UV /*Momentum mixing along lines of const density */
#define SPLINES_VDIFF
#define SPLINES_VVISC

/*Tracer advection/diffusion*/
#define TS_U3ADV_SPLIT /*Contains correction for terrain-following coordinates*/
#define MIX_GEO_TS     /*Needed with TS_U3ADV_SPLIT*/
#define TS_DIF4        /*Needed with TS_U3ADV_SPLIT*/
!#define TS_A4HADVECTION
!#define TS_A4VADVECTION


/*Defines mean flow and tidal components*/
#define ANA_M2OBC
#define ANA_FSOBC

/*restart parameters*/
#define NO_LBC_ATT /*use this to tell ROMS to ignore a restart bug*/
!#define PERFECT_RESTART

/****************************************************************
LMD Turbulence Closure Scheme*/
#define LMD_MIXING     /*Large/McWilliams/Doney scheme*/

#ifdef LMD_MIXING
#define LMD_BKPP       /*Bottom boundary layer from local KPP*/
#define LMD_SKPP       /*needed to add this or the code wouldn't compile*/
#define LMD_CONVEC     /*Convecting mixing from shear instabilities*/
#define LMD_RIMIX      /*Diffusivity due to shear instabilities*/
#define RI_SPLINES     /*Spline reconstruction of vertical shear*/
#define LMD_SHAPIRO    /*eliminates grid-scale noise in boundary layer thickness*/
#endif
/*
!*****************************************************************

!Bottom Boundary Layer Closure?
*/

/*GLS Turbulence Closure Scheme*/
!#define GLS_MIXING /*Generic length scale mixing*/
#ifdef GLS_MIXING
#define K_C4ADVECTION /*Centered 4th order advection*/
#define N2S2_HORAVG /*Horizontal smoothing of buoyancy/shear */
#endif

/*Fluxes on top and bottom of domain, all set to zero. These do not correspond to the analytical files precisely, but all of these are necessary for the code to run.*/
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX

#define ANA_M3OBC
!#define ANA_GRID   /*grid and initial condition files are generated by python*/
!#define ANA_INITIAL
!#define ANA_NUDGE   /*Climatology fields are activated in the .in file*/
!#define ANA_M2CLIMA
!#define ANA_M3CLIMA
