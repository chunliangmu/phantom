!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testapr
!
! Unit test for adaptive particle refinement
!
! :References:
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, apr, boundary, deriv, dim, dust_formation,
!   energies, eos, io, mpidomain, mpiutils, part, random, testutils,
!   unifdis, units, utils_apr
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master,fatal
 implicit none
 public :: test_apr,setup_apr_region_for_test

 private

contains

!--------------------------------------------
!+
!  Various tests of the apr module
!+
!--------------------------------------------
subroutine test_apr(ntests,npass)
 use unifdis,      only:set_unifdis
 use boundary,     only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use part,         only:npart,npartoftype,hfact,xyzh,init_part,massoftype,radprop
 use part,         only:isetphase,igas,iphase,vxyzu,fxyzu,apr_level,maxvxyzu,init_rho_from_h
 use part,         only:rho,igasP,eos_vars
 use mpidomain,    only:i_belong
 use mpiutils,     only:reduceall_mpi
 use dim,          only:periodic,use_apr,do_radiation
 use apr,          only:update_apr
 use utils_apr,    only:apr_centre,apr_max
 use energies,     only:compute_energies,angtot,etot,totmom,ekin,etherm
 use random,       only:ran2
 use eos,          only:gamma
 use deriv,        only:get_density_global
 integer, intent(inout) :: ntests,npass
 real :: psep,rhozero,time,totmass, etotin, totmomin
 real :: angtotin, ekinin, ethermin
 real :: tolang, tolen, tolmom, hcoarse,hfine,tot_entropy, rhoi, tolent, entropin
 integer :: original_npart,splitted,nfailed(13),i,iseed,ncoarse,nfine,tolpart

 if (use_apr) then
    if (id==master) write(*,"(/,a)") '--> TESTING APR MODULE'
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING APR TEST (REQUIRES -DAPR)'
    return
 endif

 ! Tolerances
 tolmom = 2.e-15
 tolang = 3.0e-14
 tolen  = 2.e-15
 tolent = 2.0E-15
 nfailed(:) = 0
 iseed = -92757
 gamma = 5./3.

 ! Set up a uniform box of particles
 call init_part()
 psep = dxbound/20.
 time = 0.
 npartoftype(:) = 0
 npart = 0
 rhozero = 1.0
 totmass = rhozero*dxbound*dybound*dzbound
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,mask=i_belong)

 original_npart = npart
 npartoftype(igas) = npart
 massoftype(igas) = totmass/reduceall_mpi('+',npart)
 iphase(1:npart) = isetphase(igas,iactive=.true.)
 call init_rho_from_h()

 ! this is to prevent a (reasonable) problem when running this test with DEBUG=yes and radiation
 if (do_radiation) then
    radprop(4,:) = 23.0421 ! just some inconsequential number
 endif

 ! Set some random velocities
 do i=1,npart
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
    if (maxvxyzu > 3) vxyzu(4,i) = ran2(iseed)**2
 enddo

 ! Here we need to initialse the pressure (for the entropy calculation later)
 do i=1,npart
    rhoi = rho(i)
    eos_vars(igasP,i) = (gamma - 1.)*rhoi*vxyzu(4,i)
 enddo

 ! Initialise APR
 call setup_apr_region_for_test()
 apr_centre(:,1:2) = 20. ! just moves the APR region away from the box so you don't have any split or merge
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 !
 ! mixed masses on a uniform lattice: neighbour counts (hence h) should match,
 ! unlike h(rho) which would split by m^{1/3}
 !
 write(*,"(/,a)") '--> checking h from number density at mixed APR levels'
 do i=1,npart
    if (mod(i,2)==0) apr_level(i) = 2
 enddo
 call get_density_global(1)
 ncoarse = 0
 nfine   = 0
 hcoarse = 0.
 hfine   = 0.
 do i=1,npart
    if (apr_level(i) == 1) then
       ncoarse = ncoarse + 1
       hcoarse = hcoarse + xyzh(4,i)
    else
       nfine = nfine + 1
       hfine = hfine + xyzh(4,i)
    endif
 enddo
 ncoarse = int(reduceall_mpi('+',ncoarse))
 nfine   = int(reduceall_mpi('+',nfine))
 hcoarse = reduceall_mpi('+',hcoarse)
 hfine   = reduceall_mpi('+',hfine)
 if (ncoarse > 0) hcoarse = hcoarse/ncoarse
 if (nfine > 0) hfine = hfine/nfine
 if (ncoarse > 0 .and. nfine > 0) then
    call checkval(hfine,hcoarse,5.e-2,nfailed(8),'h fine vs coarse (number density)')
 else
    nfailed(8) = 1
    if (id==master) write(*,*) 'ERROR: mixed APR levels not present for h(n) test'
 endif
 call update_test_scores(ntests,nfailed(8:8),npass)
 apr_level(1:npart) = 1
 call init_rho_from_h()

 ! Initialise the energies values
 call compute_energies(0.)
 call calc_entropy(npart,eos_vars,xyzh,apr_level,tot_entropy)
 etotin   = etot
 totmomin = totmom
 angtotin = angtot
 ekinin = ekin
 ethermin = etherm
 entropin = tot_entropy

 ! Now set for a split
 write(*,"(/,a)") '--> conducting a split'
 apr_centre(1:2,1) = 0.25    ! this puts a sphere centred at (0.25,0.25)
 apr_centre(1:2,2) = -0.25   ! and a second sphere at (-0.25,-0.25)
 apr_centre(3,1:2) = 0.      ! and ensures they are in the plane
 ! move the second region just for testing
 apr_centre(1:3,2) = 20.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 ! Check the new conserved values
 call compute_energies(0.)
 call calc_entropy(npart,eos_vars,xyzh,apr_level,tot_entropy)
 call checkval(angtot,angtotin,tolang,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tolmom,nfailed(2),'linear momentum')
 call checkval(etot,etotin,tolen,nfailed(3),'total energy')
 call checkval(ekin,ekinin,tolen,nfailed(4),'kinetic energy')
 call checkval(etherm,ethermin,tolen,nfailed(5),'thermal energy')
 call checkval(tot_entropy,entropin,tolent,nfailed(6),'total entropy')
 call update_test_scores(ntests,nfailed(1:6),npass)

 ! after splitting, the total number of particles should have been updated
 splitted = npart

 ! Move the apr zone out of the box and update again to merge
 write(*,"(/,a)") '--> conducting a merge'
 apr_centre(:,1:2) = 20. ! move the APR zones away again
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)

 ! Check the new conserved values
 call compute_energies(0.)
 call calc_entropy(npart,eos_vars,xyzh,apr_level,tot_entropy)
 call checkval(angtot,angtotin,tolang,nfailed(7),'angular momentum')
 call checkval(totmom,totmomin,tolmom,nfailed(8),'linear momentum')
 call checkval(etot,etotin,tolen,nfailed(9),'total energy')
 call checkval(ekin,ekinin,tolen,nfailed(10),'kinetic energy')
 call checkval(etherm,ethermin,tolen,nfailed(11),'thermal energy')
 call checkval(tot_entropy,entropin,tolent,nfailed(12),'total entropy')

 ! Check that the original particle number returns
 tolpart = (apr_max-1)*6 ! this corresponds to a maximum of 11 "spare" particles per level
 call checkval(npart,original_npart,tolpart,nfailed(13),'number of particles conserved')
 call update_test_scores(ntests,nfailed,npass)

 ! test the dust moments for dust nucleations
 call test_apr_dust(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- APR TEST COMPLETE'

end subroutine test_apr

!--------------------------------------------
!+
!  Check that APR preserves the dust moments.
!
!  The nucleation moments are number densities per hydrogen atom,
!  so the conserved quantities for particle i are
!     m_i * nucleation(idK0:idK3,i) / mass_per_H.
!  Since mass_per_H is constant, summing m_i*K_i is sufficient.
!  Splitting must duplicate K_i and merging equal-mass particles must average it.
!  The same average also conserves the mass-weighted opacity exactly.
!+
!--------------------------------------------
subroutine test_apr_dust(ntests,npass)
 use part,         only:npart,xyzh,vxyzu,fxyzu,apr_level,nucleation,n_nucleation,&
                        idJstar,idK0,idK1,idK2,idK3,idkappa,dust_temp,rho,get_pmass
 use dust_formation, only:init_nucleation,calc_kappa_dust,evolve_chem
 use apr,          only:update_apr
 use dim,          only:nucleation_supported=>nucleation,do_nucleation,inucleation,&
                        maxp_nucleation,maxTdust,store_dust_temperature
 use utils_apr,    only:apr_centre
 use allocutils,   only:allocate_array
 use units,        only:unit_density
 use io,           only:id,master
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ndust_tests = 14
 integer :: i,nfailed(ndust_tests),ndust
 real :: moments_in(0:3),moments(0:3),moments_expected(0:3)
 real :: opacity_in,opacity,pmassi,Tchem
 real :: JKmuS(n_nucleation)
 real, parameter :: dt_chem = 1.e3 ! seconds

 if (.not.nucleation_supported) return

 if (id==master) write(*,"(/,a)") '--> checking dust moment conservation across APR'

 ! allocate only what this test needs; phantomtest normally has no nucleation arrays
 if (allocated(nucleation)) deallocate(nucleation)
 maxp_nucleation = 2*npart + 16
 call allocate_array('nucleation',nucleation,n_nucleation,maxp_nucleation)
 inucleation = 1
 do_nucleation = .true.
 if (allocated(dust_temp)) deallocate(dust_temp)
 maxTdust = maxp_nucleation
 call allocate_array('dust_temp',dust_temp,maxTdust)
 store_dust_temperature = .true.

 call init_nucleation()
 ! use a deterministic, spatially varying state
 ! so that averaging and the subsequent chemistry are genuinely tested
 do i=1,npart
    nucleation(idK0,i) = 1.e-8*(1. + mod(i,5))
    nucleation(idK1,i) = 10.*nucleation(idK0,i)*(1. + 0.1*mod(i,3))
    nucleation(idK2,i) = 100.*nucleation(idK0,i)*(1. + 0.2*mod(i,7))
    nucleation(idK3,i) = 1.e-5*(1. + mod(i,11))
    nucleation(idJstar,i) = 1.e-12*(1. + mod(i,13))
    dust_temp(i) = 800. + 10.*mod(i,17)
    nucleation(idkappa,i) = calc_kappa_dust(nucleation(idK3,i),dust_temp(i),rho(i)*unit_density)
 enddo
 call get_dust_moment_sums(npart,moments_in)

 ! reference evolution: each original particle is evolved once,
 ! independent of how many APR children later represent it
 moments_expected = 0.
 do i=1,npart
    JKmuS = nucleation(:,i)
    Tchem = dust_temp(i)
    call evolve_chem(dt_chem,Tchem,rho(i)*unit_density,JKmuS)
    pmassi = get_pmass(i,.true.)
    moments_expected(0) = moments_expected(0) + pmassi*JKmuS(idK0)
    moments_expected(1) = moments_expected(1) + pmassi*JKmuS(idK1)
    moments_expected(2) = moments_expected(2) + pmassi*JKmuS(idK2)
    moments_expected(3) = moments_expected(3) + pmassi*JKmuS(idK3)
 enddo

 ! split the same region used by the main APR test
 apr_centre(1:2,1) = 0.25
 apr_centre(3,1) = 0.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 call get_dust_moment_sums(npart,moments)
 call checkval(moments(0),moments_in(0),1.e-13,nfailed(1),'dust K0 conserved across split')
 call checkval(moments(1),moments_in(1),1.e-13,nfailed(2),'dust K1 conserved across split')
 call checkval(moments(2),moments_in(2),1.e-13,nfailed(3),'dust K2 conserved across split')
 call checkval(moments(3),moments_in(3),1.e-13,nfailed(4),'dust K3 conserved across split')

 ! chemistry depends on the normalised moments, not on particle mass:
 ! evolving all split children must reproduce the single-parent reference
 do i=1,npart
    Tchem = dust_temp(i)
    call evolve_chem(dt_chem,Tchem,rho(i)*unit_density,nucleation(:,i))
    nucleation(idkappa,i) = calc_kappa_dust(nucleation(idK3,i),dust_temp(i),rho(i)*unit_density)
 enddo
 call get_dust_moment_sums(npart,moments)
 call checkval(moments(0),moments_expected(0),1.e-13,nfailed(5),'dust K0 evolved consistently after split')
 call checkval(moments(1),moments_expected(1),1.e-13,nfailed(6),'dust K1 evolved consistently after split')
 call checkval(moments(2),moments_expected(2),1.e-13,nfailed(7),'dust K2 evolved consistently after split')
 call checkval(moments(3),moments_expected(3),1.e-13,nfailed(8),'dust K3 evolved consistently after split')
 moments_in = moments
 opacity_in = 0.
 do i=1,npart
    opacity_in = opacity_in + get_pmass(i,.true.)*nucleation(idkappa,i)
 enddo

 ! move the region away again and merge
 apr_centre(1:2,1) = 20.
 call update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
 call get_dust_moment_sums(npart,moments)
 call checkval(moments(0),moments_in(0),1.e-13,nfailed(9),'dust K0 conserved across merge')
 call checkval(moments(1),moments_in(1),1.e-13,nfailed(10),'dust K1 conserved across merge')
 call checkval(moments(2),moments_in(2),1.e-13,nfailed(11),'dust K2 conserved across merge')
 call checkval(moments(3),moments_in(3),1.e-13,nfailed(12),'dust K3/dust mass conserved across merge')

 ! equal-mass averaging conserves the mass-weighted opacity exactly;
 ! recomputing kappa(<K3>,<Tdust>) would instead discard the covariance
 opacity = 0.
 do i=1,npart
    opacity = opacity + get_pmass(i,.true.)*nucleation(idkappa,i)
 enddo
 call checkval(opacity,opacity_in,1.e-13,nfailed(13),'mass-weighted opacity conserved across merge')

 ndust = 0
 do i=1,npart
    if (.not.(all(ieee_is_finite(nucleation(idK0:idK3,i))) .and. &
              ieee_is_finite(nucleation(idkappa,i)) .and. ieee_is_finite(dust_temp(i)))) ndust = ndust + 1
 enddo
 call checkval(ndust,0,0,nfailed(14),'finite dust state after APR')

 call update_test_scores(ntests,nfailed,npass)

 do_nucleation = .false.
 inucleation = 0
 store_dust_temperature = .false.

end subroutine test_apr_dust

!--------------------------------------------
!+
!  H-weighted (up to the constant 1/mass_per_H) dust moment sums
!+
!--------------------------------------------
subroutine get_dust_moment_sums(npart,moments)
 use part, only:nucleation,idK0,idK1,idK2,idK3,get_pmass
 integer, intent(in)  :: npart
 real,    intent(out) :: moments(0:3)
 integer :: i
 real :: pmassi

 moments = 0.
 do i=1,npart
    pmassi = get_pmass(i,.true.)
    moments(0) = moments(0) + pmassi*nucleation(idK0,i)
    moments(1) = moments(1) + pmassi*nucleation(idK1,i)
    moments(2) = moments(2) + pmassi*nucleation(idK2,i)
    moments(3) = moments(3) + pmassi*nucleation(idK3,i)
 enddo

end subroutine get_dust_moment_sums

!--------------------------------------------
!+
!  Set up an APR region that is used in other tests
!+
!--------------------------------------------
subroutine setup_apr_region_for_test()
 use apr,        only:init_apr,update_apr
 use utils_apr,  only:apr_type,apr_rad,apr_max_in,ref_dir,ntrack,split_dir
 use part,       only:apr_level
 integer :: ierr

 if (id==master) write(*,"(/,a)") '--> adding an apr region'

 ! set parameters for the region
 apr_max_in  =   1     ! number of additional refinement levels (3 -> 2x resolution)
 ref_dir     =   1     ! increase (1) or decrease (-1) resolution
 apr_type    =  -1     ! choose this so you get the default option which is reserved for the test suite
 apr_rad     =   0.2   ! radius of innermost region
 ntrack      =   2     ! number of regions to track
 split_dir   =   1     ! we only test random split directions as we use multiple regions

 ! initialise
 call init_apr(apr_level,ierr)

end subroutine setup_apr_region_for_test

subroutine calc_entropy(npart,eos_vars,xyzh,apr_level,tot_entropy)
 use part, only:igasP,igas,aprmassoftype,rho
 use eos,  only:gamma
 real, intent(in) :: eos_vars(:,:),xyzh(:,:)
 integer(kind=1), intent(in) :: apr_level(:)
 integer, intent(in) :: npart
 real, intent(out) :: tot_entropy
 integer :: ii
 real :: pmassi,rhoi

 tot_entropy = 0.
 do ii = 1,npart
    pmassi = aprmassoftype(igas,apr_level(ii))
    rhoi = rho(ii)
    tot_entropy = tot_entropy + (eos_vars(igasP,ii)*rhoi**(-gamma)*pmassi)
    write(10,*) ii, pmassi, eos_vars(igasP,ii),(eos_vars(igasP,ii)*rhoi**(-gamma)*pmassi)
 enddo

end subroutine calc_entropy

end module testapr
