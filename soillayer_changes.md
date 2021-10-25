## Process for adding additional soil layers to Noah-MP 
Author: Carolina Bieri, bieri2@illinois.edu

1.	Add extra variables to module_NoahMP_hrldas_driver.F. 
    - New soil variables with additional layers:
      - NEWSMOIS
      - NEWSH2O 
      - NEWTSLB
      - NEWSMOISEQ
    - Other variables necessary for adding additional layers:
      - NEWDZS
      - NEWZSNSOXY
      - ADDL_SOIL_DZ
      - ADDL_SOIL_LAYERS

Add to variable declarations:

```fortran
REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  NEWSMOIS        ! CB
REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  NEWTSLB         ! CB
REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  NEWSH2O         ! CB
REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  NEWSMOISEQ      ! CB
REAL,    ALLOCATABLE, DIMENSION(:)        ::  NEWDZS          ! CB
REAL,    ALLOCATABLE, DIMENSION(:)        ::  ADDL_SOIL_DZ    ! CB
REAL,    ALLOCATABLE, DIMENSION(:,:,:)    ::  NEWZSNSOXY      ! CB
integer                                   :: addl_soil_layers = 8  ! CB
```

Add to array allocations:
```fortran
ALLOCATE ( NEWSMOIS   (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
ALLOCATE ( NEWSH2O    (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
ALLOCATE ( NEWSMOISEQ (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
ALLOCATE ( NEWTSLB    (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
ALLOCATE ( NEWDZS     (1:NSOIL+ADDL_SOIL_LAYERS) )  ! CB
ALLOCATE ( NEWZSNSOXY (XSTART:XEND,-NSNOW+1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) )  ! CB
ALLOCATE ( ADDL_SOIL_DZ(1:ADDL_SOIL_LAYERS) )       ! CB
```

2.	Comment out array initializations to missing values (in order to run with MMF).
3.	Make changes to NOAHMP_INIT arguments in module_sf_noahmpdrv.F:

```fortran
SUBROUTINE NOAHMP_INIT ( MMINLU, SNOW , SNOWH , CANWAT , ISLTYP ,   IVGTYP, XLAT, &
       ADDL_SOIL_LAYERS, & ! CB
       TSLB , SMOIS , SH2O , DZS , FNDSOILW , FNDSNOWH ,             &
       TSK, isnowxy , tvxy     ,tgxy     ,canicexy ,         TMN,     XICE,   &
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy    ,                                                   &
       croptype ,cropcat  ,                      &
       irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
       irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
       irmivol  ,irfivol  ,irrsplh  ,            &
!jref:start
       t2mvxy   ,t2mbxy   ,chstarxy,             &
!jref:end
       NSOIL, restart,                 &
       allowed_to_read , iopt_run,  iopt_crop, iopt_irr, iopt_irrm,           &
       sf_urban_physics,                         &  ! urban scheme
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte,                &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&     ! Optional groundwater
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &      ! Optional groundwater
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &      ! Optional groundwater
       rechclim,                                                             &      ! Optional groundwater
       gecros_state,                                                         &       ! Optional gecros crop
       NEWTSLB, ADDL_SOIL_DZ, NEWDZS, NEWZSNSOXY)  ! CB

```
4. Add new variables to declaration section of NOAHMP_INIT:

```fortran
INTEGER, INTENT(IN) :: ADDL_SOIL_LAYERS ! CB

REAL, DIMENSION(1:ADDL_SOIL_LAYERS), INTENT(IN), OPTIONAL :: ADDL_SOIL_DZ ! CB
REAL, DIMENSION(ims:ime,1:nsoil+addl_soil_layers,jms:jme), INTENT(INOUT), OPTIONAL :: NEWTSLB ! CB
REAL, DIMENSION(1:nsoil+addl_soil_layers), INTENT(INOUT), OPTIONAL :: NEWDZS ! CB
REAL, DIMENSION(ims:ime,-2:nsoil+addl_soil_layers,jms:jme), INTENT(INOUT), OPTIONAL :: NEWZSNSOXY  ! CB

INTEGER                :: NEWNSOIL ! CB
REAL                      :: DEPTH  ! CB
REAL, DIMENSION(1:NSOIL+ADDL_SOIL_LAYERS) :: NEWZSOIL ! CB
```
5. Add soil temp interpolation code near end of NOAHMP_INIT:

```fortran

IF (PRESENT(NEWTSLB) .AND. &     ! CB soil temp initialization for additional layers
   PRESENT(NEWDZS) .AND. &
   PRESENT(NEWZSNSOXY) .AND. &
   PRESENT(ADDL_SOIL_DZ)) THEN

 NEWTSLB(:,1:NSOIL,:) = TSLB

 NEWDZS(1:NSOIL)  = DZS
 NEWDZS(NSOIL+1:) = ADDL_SOIL_DZ

 NEWZSNSOXY(:,-2:NSOIL,:) = ZSNSOXY

 NEWNSOIL = NSOIL+ADDL_SOIL_LAYERS

 DO J = jts, jtf
  DO I = its, itf
   DO NS = NSOIL+1, NEWNSOIL
     IF (NS .EQ. NSOIL+1) THEN
       DEPTH = ZSOIL(NS-1) - NEWDZS(NS)
     ELSE
       DEPTH = DEPTH - NEWDZS(NS)
     ENDIF

     NEWZSNSOXY(I,NS,J) = DEPTH

     NEWTSLB(I,NS,J) = NEWTSLB(I,4,J)+ &
                       ((DEPTH-ZSOIL(4))*((TMN(I,J)-TSLB(I,4,J))/(20.0-ZSOIL(4))))

   ENDDO
  ENDDO
 ENDDO

 PRINT *, NEWZSNSOXY(100,:,30)

ENDIF

```

7. Make changes accordingly to calls to NOAHMP_INIT in module_NoahMP_hrldas_driver.F.

First call to NOAHMP_INIT (only executed if restart file used):
```fortran
CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                ADDL_SOIL_LAYERS, & ! CB
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                 GRAINXY,    GDDXY,                                                   &
                CROPTYPE,  CROPCAT,                                                   &
                irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
                irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
                irmivol  ,irfivol  ,irrsplh  ,            &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .true.,                                                   &
                  .true.,runoff_option, crop_option, irrigation_option, irrigation_method,   &
                  sf_urban_physics,                         &  ! urban scheme
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                     rechclim ,gecros_state                 &
                     )
```
Second call to NOAHMP_INIT (always executed, regardless if restart run or not):

```fortran
CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, XLAT, &   ! call from WRF phys_init
                    ADDL_SOIL_LAYERS, & ! CB
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, QRAINXY,  WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                 GRAINXY,    GDDXY,                                                   &
                CROPTYPE,  CROPCAT,                                                   &
                irnumsi  ,irnummi  ,irnumfi  ,irwatsi,    &
                irwatmi  ,irwatfi  ,ireloss  ,irsivol,    &
                irmivol  ,irfivol  ,irrsplh  ,            &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .false.,                                                   &
                  .true.,runoff_option, crop_option, irrigation_option, irrigation_method,  &
                  sf_urban_physics,                         &  ! urban scheme
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                     rechclim ,gecros_state, &
                     NEWTSLB, ADDL_SOIL_DZ, NEWDZS, NEWZSNSOXY      )     ! CB

```
7. Define addl_soil_dz before second call to NOAHMP_INIT in module_NoahMP_hrldas_driver.F

```fortran
addl_soil_dz = (/1.0,1.0,2.0,2.0,2.0,2.0,3.0,5.0/) ! CB
```

8.	Redefine variables after second call to NOAHMP_INIT in module_NoahMP_hrldas_driver.F:

```fortran
call move_alloc(newzsnsoxy, zsnsoxy) ! CB
call move_alloc(newdzs, dzs)   ! CB
```

9. Make changes to GROUNDWATER_INIT.

Changes to GROUNDWATER_INIT arguments:

```fortran
SUBROUTINE GROUNDWATER_INIT (   &
            &            GRID, NSOIL , DZS, ISLTYP, IVGTYP, WTDDT , &
            &            FDEPTH, TOPO, RIVERBED, EQWTD, RIVERCOND, PEXP , AREA ,WTD ,  &
            &            SMOIS,SH2O, SMOISEQ, SMCWTDXY, DEEPRECHXY, RECHXY ,  &
            &            QSLATXY, QRFSXY, QSPRINGSXY,                  &
            &            rechclim  ,                                   &
            &            NEWSMOIS, NEWSH2O, NEWSMOISEQ, TSLB, ADDL_SOIL_LAYERS, & ! CB
            &            ids,ide, jds,jde, kds,kde,                    &
            &            ims,ime, jms,jme, kms,kme,                    &
            &            ips,ipe, jps,jpe, kps,kpe,                    &
            &            its,ite, jts,jte, kts,kte                     )

```
Additions to variable declarations in GROUNDWATER_INIT:

```fortran
INTEGER,    INTENT(IN)                           :: ADDL_SOIL_LAYERS ! CB
REAL,    INTENT(IN), DIMENSION(1:NSOIL+ADDL_SOIL_LAYERS)  :: DZS     ! CB
REAL,    INTENT(INOUT), DIMENSION(ims:ime, 1:NSOIL+ADDL_SOIL_LAYERS,jms:jme) :: NEWSMOIS   ! CB
REAL,    INTENT(INOUT), DIMENSION(ims:ime, 1:NSOIL+ADDL_SOIL_LAYERS,jms:jme) :: NEWSH2O    ! CB
REAL,    INTENT(INOUT), DIMENSION(ims:ime, 1:NSOIL+ADDL_SOIL_LAYERS,jms:jme) :: NEWSMOISEQ ! CB
REAL,    INTENT(IN), DIMENSION( ims:ime, 1:NSOIL+ADDL_SOIL_LAYERS, jms:jme) :: TSLB        ! CB

REAL :: SMCBELOW, FK ! CB
REAL, PARAMETER           :: HLICE = 3.335E5 ! CB
REAL, PARAMETER           :: GRAV = 9.81     ! CB
REAL, PARAMETER           :: T0 = 273.15     ! CB
LOGICAL :: FILLLAYERS ! CB

REAL, DIMENSION(1:NSOIL+ADDL_SOIL_LAYERS) :: SMCEQ,ZSOIL ! CB
```

10.	Add code to soil moisture initialization in GROUNDWATER_INIT:

```fortran
 ! Fill in SM where water table is below resolved soil layers - CB
 IF (WTD(I,J) < ZSOIL(NSOIL)) THEN
     DO K = NSOIL,5,-1

        ! Define DDZ
        IF ( K < NSOIL ) THEN
           DDZ = ( ZSOIL(K-1) - ZSOIL(K+1) ) * 0.5
        ELSE
           DDZ = ZSOIL(K-1) - ZSOIL(K)
        ENDIF

        ! Use Newton-Raphson method to define SM for each layer
        EXPON = BEXP +1.
        AA = DWSAT/DDZ
        BB = DKSAT / SMCMAX ** EXPON

        IF ( K .EQ. NSOIL ) THEN
          SMCBELOW = SMCWTDXY(I,J)
        ELSE
          SMCBELOW = SMC
        ENDIF

       ! First guess for SMC
        SMC = 0.5 * SMCMAX

       ! Newton-Raphson iteration assuming vertical flux = 0
        DO ITER = 1, 100
          FUNC = (SMC - SMCBELOW) * AA +  BB * SMC ** EXPON
          DFUNC = AA + BB * EXPON * SMC ** BEXP

          DX = FUNC/DFUNC
          SMC = SMC - DX
          IF ( ABS (DX) < 1.E-6)EXIT
        ENDDO
        NEWSMOIS(I,K,J) = MAX(SMC,1.E-4)

     ENDDO

 ELSEIF (WTD(I,J) .GE. ZSOIL(NSOIL)) THEN

    DO K=NSOIL,5,-1
       FILLLAYERS = .FALSE.

      ! Define SM for a given layer if the WTD is in that layer or at
      ! the bottom of that layer
       IF ((WTD(I,J) .LT. ZSOIL(K-1)) .AND. (WTD(I,J) .GT. ZSOIL(K))) THEN
          NEWSMOIS(I,K,J) = SMCMAX !* ( WTD(I,J) -  ZSOIL(K)) + &
                           !NEWSMOISEQ(I,K,J) * (ZSOIL(K-1) - WTD(I,J))
          FILLLAYERS = .TRUE.
       ELSEIF (WTD(I,J) .EQ. ZSOIL(K)) THEN
          NEWSMOIS(I,K,J) = NEWSMOISEQ(I,K,J)
          FILLLAYERS = .TRUE.
       ENDIF

 !     ! Fill in other layers
       IF ((K .GT. 5) .AND. (FILLLAYERS .EQ. .TRUE.)) THEN
         DO NS = K-1,5,-1
            DDZ = ( ZSOIL(NS-1) - ZSOIL(NS+1) ) * 0.5

            EXPON = BEXP +1.
            AA = DWSAT/DDZ
            BB = DKSAT / SMCMAX ** EXPON

            SMCBELOW = NEWSMOIS(I,NS+1,J)

            SMC = 0.5 * SMCMAX
            DO ITER = 1, 100
              FUNC = (SMC - SMCBELOW) * AA +  BB * SMC ** EXPON
              DFUNC = AA + BB * EXPON * SMC ** BEXP

              DX = FUNC/DFUNC
              SMC = SMC - DX
              IF ( ABS (DX) < 1.E-6)EXIT
            ENDDO

            NEWSMOIS(I,NS,J) = MAX(SMC,1E-4)
         ENDDO
       ENDIF
    ENDDO
 ENDIF

 ! Fill in values for SH2O based on SMC values
 DO K = NSOIL,2,-1
      IF(K .LE. 4)THEN
         IF(NEWSMOIS(I,K,J) .EQ. SMCMAX)THEN
           FRLIQ = NEWSH2O(I,K,J) / NEWSMOIS(I,K,J)
           NEWSH2O(I,K,J) = SMCMAX * FRLIQ
         ENDIF
      ELSE
         IF ( TSLB(I,K,J) < 273.149 ) THEN    ! Use explicit as initial soil ice
           FK=(( (HLICE/(GRAV*(-PSISAT))) *                              &
              ((TSLB(I,K,J)-T0)/TSLB(I,K,J)) )**(-1/BEXP) )*SMCMAX
           FK = MAX(FK, 0.02)
           NEWSH2O(I,K,J) = MIN( FK, NEWSMOIS(I,K,J) )
         ELSE
           NEWSH2O(I,K,J)=NEWSMOIS(I,K,J)
         ENDIF
      ENDIF
 END DO

```

11.	Add code to module_NoahMP_hrldas_driver.F after call to GROUNDWATER_INIT:

```fortran
call move_alloc(newsmois, smois)     ! CB
call move_alloc(newsh2o, sh2o)       ! CB
call move_alloc(newsmoiseq, smoiseq) ! CB

NSOIL = num_soil_layers
```

12.	Change/add lines in TRANSFER_MP_PARAMETERS:
```fortran
INTEGER, INTENT(IN)    :: SOILTYPE(NSOIL)
INTEGER, INTENT(IN)    ::  NSOIL     ! number of soil layers
```

Also add NSOIL to arguments in TRANSFER_MP_PARAMETERS. 

13.	Change NSOIL in module_noahmplsm.F to new value.

