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

2.	Comment out array initializations to missing values (in order to run with MMF)
3.	Make changes to NoahMP_INIT arguments:

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

