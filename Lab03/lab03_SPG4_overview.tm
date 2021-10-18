KPL/MK

KERNEL CONTENT:

 - LSK kernel: 
    naif0012.tls -> Standard leap second kernel

 - SPK kernels:
    de425s.bsp -> Planetary ephemeris
    estrack_v03.bsp -> ESA tracking station kernels. Defines the stations
                       KIRUNA1, KIRUNA2, SVALBARD, MALARGUE, SANTIAGO,
                       MALINDI, CEBREROS, NEW_NORCIA, NNORCIA2, V_FRANCA,
                       STA_MARIA, REDU, PERTH, MASPALOMAS, KOUROU

 - PCK kernels: 
    pck00010.tpc -> planetary constant kernel (TEXT kernel)
    earth_latest_high_prec.bpc -> high precision Earth planetary constants
                                 (nutation, polar motion). Binary kernel.
                                 Updated daily.

 - FK kernels:
    earth_fixed.tf -> Defines an alias for a EARTH_FIXED frame, which can 
                      be called for transformations.
                      EARTH_FIXED can map either IAU_EARTH or ITRF93.

\begindata

    PATH_VALUES = ( 'kernels' )
    PATH_SYMBOLS = ( 'KERNELS' )
    KERNELS_TO_LOAD = (
                       '$KERNELS\naif0012.tls',
                       '$KERNELS\pck00010.tpc',
                       '$KERNELS\earth_latest_high_prec.bpc',
					   '$KERNELS\earth_fixed.tf',
					   '$KERNELS\estrack_v03.tf',
					   '$KERNELS\estrack_v03.bsp',
					   '$KERNELS\de425s.bsp'
                      )

\begintext

KERNELS: 
(downloaded from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/)