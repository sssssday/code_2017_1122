
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... RESPIR.F

      subroutine respir(co2los,nlr,lyr,tcstva,resp,&
                       estatv,minerl,gromin,netmnr,newminrl)

      implicit none
      !include 'const.inc'
      !include 'param.inc'
      !include 'parfx.inc'
      !include 'zztim.inc'

! ... Argument declarations
      integer   nlr, lyr
      real      co2los, tcstva(nlr), estatv(nlr,3), minerl(mlyr,3), gromin(3), netmnr(mlyr,3)
      double precision newminrl

! ... Compute flows associated with microbial respiration.

! ... Input:
! ...   co2los      = CO2 loss associated with decomposition
! ...   nlr         = total number of layers modeled for Box A;
! ...                 (=2 for structural, metabolic, som1, som2;
! ...                  =1 for som3, wood compartments)
! ...   lyr         = layer for which respiration is being computed
! ...   tcstva(nlr) = total C (unlabeled + labeled) in layer 'lyr' of 
! ...                 Box A.  For components with only 1 layer,
! ...                 tcstva will be dimensioned (1).
! ...   time        = passed in /zztim/
! ...
! ... Transput:
! ...   cstatv          = C state variables for Box A;
! ...                     cstatv(lyr,iso) represents C in layer 'lyr',
! ...                     isotope 'iso' where iso=1 for unlabeled C
! ...                     and iso=2 for labeled C.  For components with
! ...                     only 1 layer, the 1st dimension of cstatv will
! ...                     be 1.
! ...   csrsnk          = C source/sink 
! ...   resp            = output variable
! ...   estatv          = N, P, and S variables for Box A.  In 'estatv(lyr,iel)',
! ...                     lyr represents layer and iel indicates N, P, or S.  For
! ...                     components with only 1 layer, the 1st dimension of
! ...                     estatv will be 1.
! ...   minerl(1,iel)   = labile N, P, or S in layer 1.
! ...   gromin(iel)     = gross mineralization -- stored in /comput/
! ...   netmnr(lyr,iel) = net mineralization for layer lyr (N, P, or S)
! ...                     For components with only 1 layer, the 1st
! ...                     dimension of netmnr will be 1.

! ... Local variables
      integer   iel
      real      mnrflo

! ... C flow from cstatv to CO2
      if (labtyp .eq. 2) then
        call csched(co2los,tcstva(lyr),dresp,resp)
      else
        call csched(co2los,tcstva(lyr),1.0,resp)
      endif
 
! ... Mineralization associated with respiration
      do 10 iel = 1, nelem
        mnrflo = co2los*estatv(lyr,iel)/tcstva(lyr)
        call flow(estatv(lyr,iel),minerl(SRFC,iel),time,mnrflo)

! ..... Update gross mineralization
        gromin(iel) = gromin(iel) + mnrflo

! ..... Update net mineralization
        netmnr(lyr,iel) = netmnr(lyr,iel) + mnrflo
! ..... newminrl should be updated only for nitrogen
! ..... akm via cak 07/31/01
        if (iel .eq. N) then
          newminrl = newminrl + mnrflo
        endif

10    continue
 
      return
      end
