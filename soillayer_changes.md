## Process for adding additional soil layers to Noah-MP 
Author: Carolina Bieri, bieri2@illinois.edu

### 1.	Add extra variables to module_NoahMP_hrldas_driver.F. 
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

### 2.	Comment out array initializations to missing values (in order to run with MMF)
### 3.	Make changes to NOAHMP_INIT arguments in module_sf_noahmpdrv.F:

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
### 4. Add new variables to declaration section of NOAHMP_INIT:

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
5. Make changes accordingly to calls to NOAHMP_INIT in module_NoahMP_hrldas_driver.F.

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
