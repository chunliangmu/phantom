!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module initialises the wave damping test, as per
!   Choi et al. 2009 (has been generalised for additional studies)
!
! :References:
!   Wurster, Price & Ayliffe (2014), MNRAS 444, 1104 (Section 4.1)
!   Wurster, Price & Bate (2016), MNRAS 457, 1037
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Bxin      : *Initial x-magnetic field*
!   - amplitude : *Initial wave amplitude*
!   - geo_cp    : *Using close-packed grid (F: cubic).*
!   - isowave   : *Modelling a sound wave (F: Alfven wave)*
!   - kwave     : *Wavenumber (k/pi)*
!   - kx_kxy    : *Using wavenumber in x only (F: initialise in x,y)*
!   - nx        : *Particles in the x-direction*
!   - polyk     : *Initial polyk*
!   - realvals  : *Using physical units (F: arbitrary units)*
!   - rect      : *Using rectangular cp grid (F: cubic cp grid)*
!   - rhoin     : *Initial density*
!   - use_ambi  : *Test ambipolar diffusion?*
!   - use_hall  : *Test the Hall effect?*
!   - use_ohm   : *Test Ohmic resistivity?*
!   - vx_vz     : *Using velocity in x (F: initialise in z)*
!
! :Dependencies: boundary, dim, infile_utils, io, mpidomain, nicil,
!   options, part, physcon, prompting, setup_params, timestep, unifdis,
!   units
!
 use part,  only:mhd
 use nicil, only:use_ohm,use_hall,use_ambi
 use nicil, only:eta_constant,eta_const_type,C_OR,C_HE,C_AD,icnstphys,icnstsemi,icnst
 use nicil, only:n_e_cnst,rho_i_cnst,rho_n_cnst,gamma_AD,alpha_AD,hall_lt_zero
 implicit none
 integer, private :: nx
 real,    private :: kwave,amplitude,polykin,rhoin0,Bxin0
 logical, private :: realvals,geo_cp,rect,oned
 logical, private :: isowave,kx_kxy,vx_vz

 public :: setpart

 private

contains
!----------------------------------------------------------------
!+
!  Setup for wave damping test
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master,fatal
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:set_particle_type,igas,Bxyz,periodic
 use timestep,     only:tmax,dtmax
 use options,      only:nfulldump
 use physcon,      only:pi,fourpi,solarm,c,qe
 use units,        only:set_units,unit_density,unit_Bfield,umass,udist
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use prompting,    only:prompt
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=100)               :: inname
 integer                          :: i,idir,ierr
 real                             :: totmass,deltax,deltay,deltaz
 real                             :: x_min,x_max,y_min,y_max,z_min,z_max
 real                             :: length,rhoin,Bxin
 real                             :: uuzero,kx,ky,xi,yi
 real                             :: Bxini,Byini,Bzini,vA,vcoef
 logical                          :: iexist
 !
 !--in-file parameters (only if not already done so)
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 time        = 0.0
 if (.not. iexist) then
    tmax      = 5.0
    dtmax     = 0.01
    nfulldump = 10
    !--Turn on constant, uncalculated resistivities
    eta_constant   = .true.
    eta_const_type = icnstsemi
 endif
 !
 !--Turn off non-ideal MHD by default; each will have the option to be
 !  modified either by the interactive script or the .setup file
 !
 use_ohm  = .false.
 use_hall = .false.
 use_ambi = .false.

 !
 !--Default runtime parameters
 !
 realvals  = .false.
 oned      = .true.
 isowave   = .false.
 kx_kxy    = .true.
 vx_vz     = .false.
 geo_cp    = .true.
 rect      = .false.
 nx        = 64
 amplitude = 0.01
 kwave     = 2.0
 use_ambi  = .true.
 if (.not. mhd) isowave = .true.

 if (oned) then
    rhoin0  = 1.0
    Bxin0   = 1.0
    polykin = 1.0
 else
    rhoin0  = 3.00
    Bxin0   = 1.30
    polykin = 0.97
 endif

 !
 !--Read runtime parameters from setup file
 !
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Wave damping test'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 if (kwave <=0) call fatal('setup','k > 0 is required')
 !
 !--Convert to real values
 !
 if (realvals) then
    !--matches dimensions of sphereinbox
    udist = 1.d16
    umass = solarm
    call set_units(dist=udist,mass=umass,G=1.d0)
    rhoin = rhoin0 * 1.0d-18
    Bxin  = Bxin0  * 1.0d-4
 else
    rhoin = rhoin0
    Bxin  = Bxin0
 endif
 if (isowave) Bxin = 0.0
 rhoin    = rhoin             / unit_density
 Bxin     = Bxin*sqrt(fourpi) / unit_Bfield
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma  = 1.
 endif
 npart          = 0
 npart_total    = 0
 npartoftype(:) = 0
 call set_boundary()
 deltax = dxbound/nx
 !
 !--Change box ratio if closepacked geometry
 !
 if ( geo_cP ) then
    deltay = deltax*sqrt(3.0)/2.0
    deltaz = deltax*sqrt(6.0)/3.0
    if ( rect ) then
       x_min = -(0.5*nx)*deltax
       x_max =  (0.5*nx)*deltax
       y_min =  -6.0*deltay
       y_max =   6.0*deltay
       z_min =  -6.0*deltaz
       z_max =   6.0*deltaz
    else
       x_min = -0.5
       x_max =  0.5
       y_min = x_min*sqrt(3.0)/2.0
       y_max = x_max*sqrt(3.0)/2.0
       z_min = x_min*sqrt(6.0)/3.0
       z_max = x_max*sqrt(6.0)/3.0
    endif
    call set_boundary(x_min,x_max,y_min,y_max,z_min,z_max)
 endif
 !
 !--Set remaining values
 !
 rhozero  = rhoin
 if (maxvxyzu < 4) then
    polyk = polykin**2
 else
    polyk = 0.
 endif
 if (maxvxyzu >= 4) then
    if (gamma > 1.) then
       uuzero = polykin/sqrt(gamma*(gamma-1.))
    else
       uuzero = 3./2.*polykin
    endif
 endif
 if ( geo_cp ) then
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                     hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 else
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                     hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 endif
 npartoftype(igas) = npart
 totmass           = rhozero*dxbound*dybound*dzbound
 massoftype(igas)  = totmass/npartoftype(igas)

 length = xmax-xmin
 kx     = kwave*pi/(xmax-xmin)
 ky     = kwave*pi/(ymax-ymin)
 !
 !--Choi et al. 2009 damping test, Bx<>0, By=Bz=0.
 !
 Bxini = Bxin
 Byini = 0.0
 Bzini = 0.0
 vA    = sqrt((Bxini**2 + Byini**2 + Bzini**2)/rhozero)
 !
 !--Determine velocity coefficient
 !
 if ( isowave ) then
    vcoef = amplitude*polykin
 else
    vcoef = amplitude*vA
 endif
 !--Determine direction of wave
 if ( vx_vz ) then
    idir = 1
 else
    idir = 3
 endif
 !
 !--Set remaining particle properties
 do i=1,npart
    call set_particle_type(i,igas)
    xi = xyzh(1,i)-xmin
    yi = xyzh(2,i)-ymin
    vxyzu(:,i) = 0.
    if ( kx_kxy ) then
       vxyzu(idir,i) = vcoef*sin(kx*xi)
    else
       vxyzu(idir,i) = vcoef*sin(kx*xi + ky*yi)
    endif
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(1,i) = Bxini
       Bxyz(2,i) = Byini
       Bxyz(3,i) = Bzini
    endif
 enddo
 if (mhd) ihavesetupB = .true.
 !
 !--Print statements
 !
 write(*,*) "setup: rho_0         = ",rhozero
 write(*,*) "setup: B_(x,0)       = ",Bxin
 write(*,*) "setup: polyk         = ",polyk
 write(*,*) "setup: k/pi          = ",kwave
 write(*,*) "setup: npart         = ",npart
 write(*,*) "setup: particle mass = ",massoftype(igas)
 write(*,*) "setup: total mass    = ",totmass
 write(*,*) "setup: total volume  = ",dxbound*dybound*dzbound

end subroutine setpart

!------------------------------------------------------------------------
!+
!  write options to .setup file
!+
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,        only:  write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wave damping setup routines'

 write(iunit,"(/,a)") '# units and orientation'
 call write_inopt(realvals,'realvals','Using physical units (F: arbitrary units)',iunit)
 if (mhd) call write_inopt(isowave,'isowave','Modelling a sound wave (F: Alfven wave)',iunit)
 call write_inopt(kx_kxy,'kx_kxy','Using wavenumber in x only (F: initialise in x,y)',iunit)
 call write_inopt(vx_vz,'vx_vz','Using velocity in x (F: initialise in z)',iunit)

 write(iunit,"(/,a)") '# Grid setup'
 call write_inopt(geo_cp,'geo_cp','Using close-packed grid (F: cubic).',iunit)
 if (geo_cp) call write_inopt(rect,'rect','Using rectangular cp grid (F: cubic cp grid)',iunit)
 call write_inopt(nx,'nx','Particles in the x-direction',iunit)
 call write_inopt(rhoin0,'rhoin','Initial density',iunit)
 call write_inopt(polykin,'polyk','Initial polyk',iunit)
 if (mhd  .and. .not. isowave) call write_inopt(Bxin0,'Bxin','Initial x-magnetic field',iunit)
 call write_inopt(amplitude,'amplitude','Initial wave amplitude',iunit)
 call write_inopt(kwave,'kwave','Wavenumber (k/pi)',iunit)

 write(iunit,"(/,a)") '# Test problem and values'
 if (mhd) then
    call write_inopt(use_ambi,'use_ambi','Test ambipolar diffusion?',iunit)
    call write_inopt(use_hall,'use_hall','Test the Hall effect?',iunit)
    call write_inopt(use_ohm, 'use_ohm', 'Test Ohmic resistivity?',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!+
!  read options from .setup file
!+
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer,          parameter   :: iunit = 21
 type(inopts),     allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(realvals,'realvals',db,ierr)
 if (mhd) call read_inopt(isowave,'isowave',db,ierr)
 call read_inopt(kx_kxy,'kx_kxy',db,ierr)
 call read_inopt(vx_vz,'vx_vz',db,ierr)
 call read_inopt(geo_cp,'geo_cp',db,ierr)
 if (geo_cp) call read_inopt(rect,'rect',db,ierr)
 call read_inopt(nx,'nx',db,ierr)
 call read_inopt(rhoin0,'rhoin',db,ierr)
 call read_inopt(polykin,'polyk',db,ierr)
 if (mhd .and. .not. isowave) call read_inopt(Bxin0,'Bxin',db,ierr)
 call read_inopt(amplitude,'amplitude',db,ierr)
 call read_inopt(kwave,'kwave',db,ierr)
 if (mhd) then
    call read_inopt(use_ambi,'use_ambi',db,ierr)
    call read_inopt(use_hall,'use_hall',db,ierr)
    call read_inopt(use_ohm, 'use_ohm', db,ierr)
 endif
 call close_db(db)

end subroutine read_setupfile

!
!---Interactive setup-----------------------------------------------------
!
subroutine setup_interactive()
 use prompting, only:prompt
 use part,      only:maxp

 call prompt('Initialise units to physical values [no: arbitrary values].',realvals)
 call prompt('Initialise values to unity [no: prime numbers].',oned)
 if (.not. oned) then
    rhoin0  = 3.00
    Bxin0   = 1.30
    polykin = 0.97
 else
    rhoin0  = 1.0
    Bxin0   = 1.0
    polykin = 1.0
 endif
 if (mhd) call prompt('Model a sound wave? [no: Alfven wave].',isowave)
 call prompt('Initialise wavenumber in x only [no: initialise in x,y].',kx_kxy)
 call prompt('Initialise velocity in x [no: initialise in z].',vx_vz)
 call prompt('Initialise the grid as close-packed [no: cubic].',geo_cp)
 if (geo_cp) call prompt('Initialise the close-packed as rectangular [no: cubic].',rect)
 call prompt('Set number of particles in the x-direction.',nx,1,int((maxp**(1.0/3.0))))
 call prompt('Set initial density.',rhoin0)
 call prompt('Set polyk.',polykin)
 if (.not. isowave) call prompt('Set initial x-magnetic field.',Bxin0)
 call prompt('Set the initial wave amplitude.',amplitude)
 call prompt('Set wavenumber (k/pi).',kwave)
 if (mhd) then
    call prompt('Pick Coefficient type: 1:phys.cnst+B+rho, 2:C_NI+B+rho, 3:C_NI',eta_const_type,1,3)
    call prompt('Test Ohmic resistivity?',use_ohm)
    call prompt('Test the Hall effect?',use_hall)
    call prompt('Test ambipolar diffusion?',use_ambi)
    if (eta_const_type==icnstsemi .or. eta_const_type==icnst) then
       if (use_ohm ) call prompt('Set C_OR.',C_OR)
       if (use_hall) call prompt('Set C_HE.',C_HE)
       if (use_ambi) call prompt('Set C_AD.',C_AD)
    elseif (eta_const_type==icnstphys) then
       if (use_ohm .or. use_hall) then
          call prompt('Set the electron number density (cgs).',n_e_cnst)
       endif
       if (use_hall) then
          call prompt('Is eta_hall < 0?',hall_lt_zero )
       endif
       if (use_ambi) then
          call prompt('Set the collisional coupling coefficient (cgs) ',gamma_AD)
          call prompt('Set ionised gas density (cgs).',rho_i_cnst)
          call prompt('Set neutral gas density (cgs).',rho_n_cnst)
          call prompt('Set power-law exponent.',alpha_AD)
       endif
    else
       print*, "setup_wavedamp: this option should not be possible"
    endif
 endif

end subroutine setup_interactive

end module setup
