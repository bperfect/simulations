      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id: ana_m2obc.h 751 2015-01-07 22:56:36Z arango $
!!======================================================================
!! Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m2obc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(12)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
     &                           angler, h, pm, pn, on_u,               &
#ifdef MASKING
     &                           umask,                                 &
#endif
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: angle, cff, fac, major, minor, omega, phase, val
      real(r8) :: ramp
#if defined ESTUARY_TEST || defined INLET_TEST
      real(r8) :: my_area, my_flux, tid_flow, riv_flow, cff1, cff2,     &
     &            model_flux
#endif
#if defined TEST_CHAN
      real(r8) :: my_area, my_width
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined ESTUARY_TEST
      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        cff1=0.40_r8                                          ! west end
        cff2=0.08_r8
        riv_flow=cff2*300.0_r8*5.0_r8
        tid_flow=cff1*300.0_r8*10.0_r8
        my_area=0.0_r8
        my_flux=0.0_r8
        DO j=JstrP,JendP
          cff=0.5_r8*(zeta(Istr  ,j,knew)+h(Istr  ,j)+                  &
     &                zeta(Istr-1,j,knew)+h(Istr-1,j))/pn(Istr,j)
          my_area=my_area+cff
        END DO
        my_flux=-tid_flow*SIN(2.0_r8*pi*time(ng)/                       &
     &          (12.0_r8*3600.0_r8))-riv_flow
        DO j=JstrP,JendP
          BOUNDARY(ng)%ubar_west(j)=my_flux/my_area
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        cff2=0.08_r8                                          ! east end
        riv_flow=cff2*300.0_r8*5.0_r8
        my_area=0.0_r8
        my_flux=0.0_r8
        DO j=JstrP,JendP
          cff=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                  &
     &                zeta(Iend+1,j,knew)+h(Iend+1,j))/pn(Iend,j)
          my_area=my_area+cff
        END DO
        my_flux=-riv_flow
        DO j=JstrP,JendP
          BOUNDARY(ng)%ubar_east(j)=my_flux/my_area
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF

#elif defined KELVIN
      fac=1.0_r8                                ! zeta0
      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
      val=fac*SIN(omega*time(ng))
      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          cff=SQRT(g*GRID(ng)%h(Istr-1,j))
          BOUNDARY(ng)%ubar_west(j)=(val*cff/GRID(ng)%h(Istr-1,j))*     &
     &                              EXP(-GRID(ng)%f(Istr-1,j)*          &
     &                                   GRID(ng)%yp(Istr-1,j)/cff)
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          cff=SQRT(g*GRID(ng)%h(Iend,j))
          val=fac*EXP(-GRID(ng)%f(Iend,j)*GRID(ng)%yp(Istr-1,j)/cff)
          BOUNDARY(ng)%ubar_east(j)=(val*cff/GRID(ng)%h(Iend,j))*       &
     &                              SIN(omega*GRID(ng)%xp(Iend,j)/cff-  &
     &                                  omega*time(ng))
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
#elif defined FIEBERLING
        fac = TANH((tdays(ng)-dstart)/10.0_r8)
!Western Edge
        IF (LBC(iwest,isUbar,ng)%acquire.and. &
      & DOMAIN(ng)%Western_Edge(tile)) THEN
                DO j=JstrP,JendP
                    BOUNDARY(ng)%ubar_west(j)=0.1_r8*fac
                END DO
      END IF
        IF (LBC(iwest,isVbar,ng)%acquire.and. &
      & DOMAIN(ng)%Western_Edge(tile)) THEN
                DO j=JstrP,JendP
                    BOUNDARY(ng)%vbar_west(j)=0.0_r8
                END DO
      END IF
!Eastern Edge
        IF (LBC(ieast,isUbar,ng)%acquire.and. &
      & DOMAIN(ng)%Eastern_Edge(tile)) THEN
                DO j=JstrP,JendP
                    BOUNDARY(ng)%ubar_east(j)=0.1_r8*fac
                END DO
      END IF
        IF (LBC(ieast,isVbar,ng)%acquire.and. &
      & DOMAIN(ng)%Eastern_Edge(tile)) THEN
                DO j=JstrP,JendP
                    BOUNDARY(ng)%vbar_east(j)=0.0_r8
                END DO
      END IF
!Northern Edge
        IF (LBC(inorth,isUbar,ng)%acquire.and. &
      & DOMAIN(ng)%Northern_Edge(tile)) THEN
                DO j=IstrP,IendT
                    BOUNDARY(ng)%ubar_north(j)=0.1_r8*fac
                END DO
      END IF
        IF (LBC(inorth,isVbar,ng)%acquire.and. &
      & DOMAIN(ng)%Northern_Edge(tile)) THEN
                DO j=IstrT,IendT
                    BOUNDARY(ng)%vbar_north(j)=0.0_r8
                END DO
      END IF
!Southern Edge
        IF (LBC(isouth,isUbar,ng)%acquire.and. &
      & DOMAIN(ng)%Southern_Edge(tile)) THEN
                DO j=IstrP,IendT
                    BOUNDARY(ng)%ubar_south(j)=0.1_r8*fac
                END DO
      END IF

        IF (LBC(isouth,isVbar,ng)%acquire.and. &
      & DOMAIN(ng)%Southern_Edge(tile)) THEN
                DO j=IstrT,IendT
                    BOUNDARY(ng)%vbar_south(j)=0.0_r8
                END DO
      END IF



!        IF (LBC(ieast,isUbar,ng)%acquire.and. &
!      & LBC(ieast,isVbar,ng)%acquire.and. &
!      & DOMAIN(ng)%Eastern_Edge(tile)) THEN
!                DO j=JstrT,JendT
!                    BOUNDARY(ng)%ubar_east(j)=0.1_r8*fac
!                END DO
!                DO j=JstrP,JendT
!                    BOUNDARY(ng)%vbar_east(j)=0.0_r8
!                END DO
!      END IF



#else
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF

      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_west(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(isouth,isUbar,ng)%acquire.and.                            &
     &    LBC(isouth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_south(i)=0.0_r8
        END DO
      END IF

      IF (LBC(inorth,isUbar,ng)%acquire.and.                            &
     &    LBC(inorth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE ana_m2obc_tile
