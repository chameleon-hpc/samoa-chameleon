! This header file maps the chosen precision to ISO and HDF5 data sizes

# if defined(_SINGLE_PRECISION)
#   define XDMF_ISO_P REAL32
#   define XDMF_HDF_P H5T_NATIVE_REAL
#   define XDMF_XMF_P 4
# elif defined(_DOUBLE_PRECISION) || defined(_QUAD_PRECISION)
#   define XDMF_ISO_P REAL64
#   define XDMF_HDF_P H5T_NATIVE_DOUBLE
#   define XDMF_XMF_P 8
# endif

! This is needed to prevent module dependency conflicts
! Note that ADER-DG still uses the old patch format, as opposed to master, which uses the new one (Tools_patch)

# if defined(_SWE_PATCH)
#   define _XDMF_PATCH
#   define _XDMF_PATCH_ORDER_SQUARE _SWE_PATCH_ORDER_SQUARE 
# endif