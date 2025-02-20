RELEM, IZR, DELEM, IZD, LVLD, ....... /LPARMS/  1a2,1x,1a2,5x,1a2,1x,1a2,2x,1a1
do until NENER.eq.-1 
    NENER 
    NMIN 
    NMAX 
    (ENER(IE), IE = 1, NENER) 
    (BETA(IE), I = 1, NENER) 
    if LPARMS is set (LTYP(IE), I = 1, NENER) (XLCR(IE), I = 1, NENER) (PL2A(IE), I = 1, NENER) (PL3A(IE), I = 1, NENER) endif (OMTOT(IE), IE = 1, NENER) for IN = NMIN to NMAX (OMN(IN,IE), IE = 1, NENER) for IL = 0 to IN-1 (OMNL(IN,IL,IE), IE=1,NENER) [ for IM= 0 to IL repeat repeat repeat repeat-1 1x,1i4 1x,1i4 1x,1i4 10x,9f9.2 10x,9f9.2 10,9i9 10x,9f9.2 10x,9f9.2 10x,9f9.2 10x,1p,9d9.2 10x,1p,9d9.2 10x,1p,9d9.2 (OMNLM(IN,IL,IM,IE,IE=1,NENER) ]     -1 10x,1p,9d