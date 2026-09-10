!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module directsum
!
! This module computes self-gravity by direct summation
!  over the particles - used to test the treecode.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, kernel, options, part
!
 implicit none
 public :: directsum_grav

 private

contains

!----------------------------------------------------------------------------
! Calculates the solution to the gravitational Poisson equation
!
! \nabla^2 \phi = 4 *pi*G \rho
!
! With force softening from the kernel
! by a direct summation over the particles.
!
! Use this to check the accuracy of the tree code
!
! Input:
!
!   xyzh(ndim,ntot)  : coordinates and smoothing length of the particles
!   gradh(ngradh,ntot) : 1/OmegaTilde, hydro zeta, and gradsoft (with gravity)
!
! Output:
!
!   phitot          : total potential phi
!   fgrav(ndim,ntot) : gravitational force
!----------------------------------------------------------------------------

subroutine directsum_grav(xyzh,gradh,fgrav,phitot,ntot)
 use kernel,    only:grkern,kernel_softening,radkern2,cnormk,cnormk_tilde, &
                     get_kernel_tilde
 use part,      only:igas,iamtype,maxphase,maxp,iphase, &
                     iactive,isdead_or_accreted,massoftype,maxgradh, &
                     apr_level,aprmassoftype
 use dim,       only:maxvxyzu,maxp,use_apr,igradsoft,igradomega
 use options,   only:two_kernel
 use io,        only:error
 integer,      intent(in)    :: ntot
 real,         intent(in)    :: xyzh(4,ntot)
 real(kind=4), intent(in)    :: gradh(:,:)
 real,         intent(inout) :: fgrav(maxvxyzu,ntot)
 real,         intent(out)   :: phitot
 integer :: i,j,iamtypei,iamtypej
 real :: dx(3),dr(3),fgravi(3),xi(3)
 real :: rij,rij1,rij2,rij21,pmassi,pmassj
 real :: gradhi,gradsofti,grkerni,grkernj,dsofti,dsoftj
 real :: grkern_tildei,grkern_tildej,wtilde,softomegaj
 real :: phii,phij,phiterm,fmi,fmj,phitemp,potensoft0,qi,qj
 real :: hi,hj,hi1,hj1,hi21,hj21,hi41,hj41,q2i,q2j
 real :: fgrav_pair
 logical :: iactivei
!
!--reset potential (but not force) initially
!
 phitot = 0.

 iamtypei = igas
 iamtypej = igas
 pmassi = massoftype(iamtypei)
 pmassj = massoftype(iamtypej)

 call kernel_softening(0.,0.,potensoft0,fmi)
 if (size(gradh(:,1)) < 3) then
    call error('directsum','cannot do direct sum with ngradh < 3')
    return
 endif
 if (maxgradh /= maxp) then
    call error('directsum','gradh not stored (maxgradh /= maxp in part.F90)')
    return
 endif
!
!--one-sided N x N sum: force on i from all j /= i (OpenMP over i)
!
!$omp parallel do default(none) schedule(static) &
!$omp shared(ntot,xyzh,gradh,fgrav,iphase,massoftype,apr_level,aprmassoftype) &
!$omp shared(maxphase,maxp,two_kernel) &
!$omp firstprivate(potensoft0) &
!$omp private(i,j,xi,hi,hi1,hi21,hi41,iamtypei,iactivei,pmassi,gradhi,gradsofti) &
!$omp private(fgravi,phitemp,dx,dr,hj,hj1,hj21,hj41,rij2,rij,rij1,rij21) &
!$omp private(iamtypej,pmassj,q2i,q2j,dsofti,dsoftj,qi,qj) &
!$omp private(grkerni,grkernj,grkern_tildei,grkern_tildej,wtilde) &
!$omp private(phii,phij,fmi,fmj,softomegaj,phiterm,fgrav_pair) &
!$omp reduction(+:phitot)
 overi: do i=1,ntot
    xi(1:3) = xyzh(1:3,i)
    hi      = xyzh(4,i)
    if (isdead_or_accreted(hi)) cycle overi

    iamtypei = igas
    iactivei = .true.
    pmassi = massoftype(igas)
    if (maxphase==maxp) then
       iamtypei = iamtype(iphase(i))
       iactivei = iactive(iphase(i))
       if (use_apr) then
          pmassi = aprmassoftype(iamtypei,apr_level(i))
       else
          pmassi = massoftype(iamtypei)
       endif
    elseif (use_apr) then
       pmassi = aprmassoftype(igas,apr_level(i))
    endif
    hi1  = 1./hi
    hi21 = hi1*hi1
    hi41 = hi21*hi21
    gradhi    = real(gradh(igradomega,i))
    gradsofti = real(gradh(igradsoft,i))
    fgravi(:) = 0.
    phitemp   = 0.

    overj: do j=1,ntot
       if (j==i) cycle overj
       dx(1) = xi(1) - xyzh(1,j)
       dx(2) = xi(2) - xyzh(2,j)
       dx(3) = xi(3) - xyzh(3,j)
       hj    = xyzh(4,j)
       if (isdead_or_accreted(hj)) cycle overj
       hj1   = 1./hj
       rij2  = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
       rij   = sqrt(rij2)
       rij1  = 1./rij
       rij21 = rij1*rij1
       dr(1) = dx(1)*rij1
       dr(2) = dx(2)*rij1
       dr(3) = dx(3)*rij1
       hj21  = hj1*hj1
       hj41  = hj21*hj21
       pmassj = massoftype(igas)
       if (maxphase==maxp) then
          iamtypej = iamtype(iphase(j))
          if (use_apr) then
             pmassj = aprmassoftype(iamtypej,apr_level(j))
          else
             pmassj = massoftype(iamtypej)
          endif
       elseif (use_apr) then
          pmassj = aprmassoftype(igas,apr_level(j))
       endif
       q2i = rij2*hi21
       q2j = rij2*hj21
       dsofti = 0.
       dsoftj = 0.
       if (q2i < radkern2) then
          qi = sqrt(q2i)
          grkerni = cnormk*grkern(q2i,qi)*hi41
          if (two_kernel) then
             call get_kernel_tilde(q2i,qi,wtilde,grkern_tildei)
             grkern_tildei = grkern_tildei*hi41*cnormk_tilde
          else
             grkern_tildei = grkerni
          endif
          call kernel_softening(q2i,qi,phii,fmi)
          phii   = phii*hi1
          fmi    = fmi*hi21
          dsofti = gradsofti*grkern_tildei*gradhi
       else
          phii = -rij1
          fmi  = rij21
       endif
       if (q2j < radkern2) then
          qj = sqrt(q2j)
          grkernj = cnormk*grkern(q2j,qj)*hj41
          if (two_kernel) then
             call get_kernel_tilde(q2j,qj,wtilde,grkern_tildej)
             grkern_tildej = grkern_tildej*hj41*cnormk_tilde
          else
             grkern_tildej = grkernj
          endif
          call kernel_softening(q2j,qj,phij,fmj)
          phij = phij*hj1
          fmj  = fmj*hj21
          softomegaj = real(gradh(igradsoft,j))*real(gradh(igradomega,j))
          dsoftj = softomegaj*grkern_tildej
       else
          phij = -rij1
          fmj  = rij21
       endif

       phiterm = 0.5*(phii + phij)
       phitemp = phitemp + pmassj*phiterm

       ! force on i from j (same pair assembly as force.F90 / previous i<j loop)
       if (iactivei) then
          fgrav_pair = 0.5*pmassj*(fmi + fmj) + 0.5*(dsofti + dsoftj*(pmassj/pmassi))
          fgravi(1) = fgravi(1) - fgrav_pair*dr(1)
          fgravi(2) = fgravi(2) - fgrav_pair*dr(2)
          fgravi(3) = fgravi(3) - fgrav_pair*dr(3)
       endif
    enddo overj

    if (iactivei) then
       fgrav(1,i) = fgrav(1,i) + fgravi(1)
       fgrav(2,i) = fgrav(2,i) + fgravi(2)
       fgrav(3,i) = fgrav(3,i) + fgravi(3)
    endif
    ! one-sided pair sum / 2 matches original unordered-pair potential;
    ! self term also carries the conventional 1/2
    phitot = phitot + 0.5*pmassi*phitemp + 0.5*pmassi*pmassi*potensoft0*hi1
 enddo overi
!$omp end parallel do

end subroutine directsum_grav

end module directsum
