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

Added to variable declarations:

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

Added to array allocations:
```fortran
  ALLOCATE ( NEWSMOIS   (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
  ALLOCATE ( NEWSH2O    (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
  ALLOCATE ( NEWSMOISEQ (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
  ALLOCATE ( NEWTSLB    (XSTART:XEND,1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) ) ! CB
  ALLOCATE ( NEWDZS     (1:NSOIL+ADDL_SOIL_LAYERS) )  ! CB
  ALLOCATE ( NEWZSNSOXY (XSTART:XEND,-NSNOW+1:NSOIL+ADDL_SOIL_LAYERS,YSTART:YEND) )  ! CB
  ALLOCATE ( ADDL_SOIL_DZ(1:ADDL_SOIL_LAYERS) )       ! CB
```
