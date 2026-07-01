!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr
!
! Everything needed for live adaptive particle refinement
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr_region, dim, get_apr_level, io, io_summary, kdtree,
!   mpiforce, neighkdtree, part, physcon, quitdump, relaxem, utils_apr,
!   vectorutils
!
 use dim, only:gr,use_apr
 use apr_region
 use utils_apr

 implicit none

 public :: init_apr,update_apr
 public :: use_apr

 private
 real    :: sep_factor = 0.2
 logical :: apr_verbose = .false.
 logical :: do_relax = .false.
 logical :: adjusted_split = .true.

contains

!-----------------------------------------------------------------------
!+
!  Initialising all the apr arrays and properties
!+
!-----------------------------------------------------------------------
subroutine init_apr(apr_level,ierr)
 use part,          only:npart,massoftype,aprmassoftype
 use apr_region,    only:set_apr_centre,set_apr_regions
 use utils_apr,     only:ntrack_max
 use get_apr_level, only:set_get_apr
 use io_summary,    only:print_apr,iosum_apr
 use io,            only:warning,fatal
 use dim,           only:maxvxyzu
 integer,         intent(inout) :: ierr
 integer(kind=1), intent(inout) :: apr_level(:)
 logical :: previously_set
 integer :: i

 ! the resolution levels are in addition to the base resolution
 apr_max = apr_max_in + 1
 if (split_dir == 2) do_relax = .true.

 ! if we're reading in a file that already has the levels set,
 ! don't override these
 previously_set = .false.
 if (sum(int(apr_level(1:npart))) > npart) then
    previously_set = .true.
    if (split_dir /= 2) do_relax = .false.
 endif

 if (.not.previously_set) then
    ! initialise the base resolution level
    if (ref_dir == 1) then
       apr_level(1:npart) = int(1,kind=1)
    else
       apr_level(1:npart) = int(apr_max,kind=1)
    endif

    ! also set the massoftype array
    ! if we are derefining we make sure that
    ! massoftype(igas) is associated with the
    ! largest particle (don't do it twice accidentally!)
    if (ref_dir == -1) then
       massoftype(:) = massoftype(:) * 2.**(apr_max -1)
       top_level = 1
    else
       top_level = apr_max
    endif
 endif

 ! now set the aprmassoftype array, this stores all the masses for the different resolution levels
 do i = 1,apr_max
    aprmassoftype(:,i) = massoftype(:)/(2.**(i-1))
 enddo

 ! how many regions do we need
 if (apr_type == 3) then
    ntrack_max = 999
    ntrack = 0 ! to start with
 elseif (apr_type == -1) then
    ntrack_max = 2
 else
    ntrack_max = 1
 endif

 if ((ntrack_max > 1) .and. (split_dir /= 3)) then
    split_dir = 3 ! no directional splitting for creating/multiple regions
    call warning('init_apr','resetting split_dir=3 because using multiple regions')
 endif

 allocate(apr_centre(3,ntrack_max),track_part(ntrack_max))
 apr_centre(:,:) = 0.

 ! initialise the shape of the region
 call set_get_apr()

 ! initiliase the regions
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 if (.not.allocated(apr_regions)) allocate(apr_regions(apr_max),npart_regions(apr_max))
 call set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
 npart_regions = 0
 icentre = 1 ! to initialise

 ! certain splitdir need certain things
 if (maxvxyzu < 4 .and. split_dir == 2) then
    call fatal('init_apr','split_dir == 2 not compatible with choice of eos')
 endif

 ierr = 0

 ! print summary please
 print_apr = .true.
 iosum_apr(1) = ntrack
 iosum_apr(2) = apr_max

 if (apr_verbose) print*,'initialised apr'

end subroutine init_apr

!-----------------------------------------------------------------------
!+
!  Subroutine to check if particles need to be split or merged
!+
!-----------------------------------------------------------------------
subroutine update_apr(npart,xyzh,vxyzu,fxyzu,apr_level)
!$ use omp_lib
 use dim,        only:maxp,ind_timesteps,maxvxyzu
 use part,       only:ntot,isdead_or_accreted,igas,aprmassoftype,&
                    shuffle_part,iphase,iactive,maxp,npartoftype
 use part,       only:igasP,rhoh,eos_vars,iorig
 use quitdump,   only:quit
 use relaxem,    only:relax_particles
 use utils_apr,  only:find_closest_region,icentre
 use apr_region, only:set_apr_centre
 use io,         only:fatal
 use get_apr_level, only:get_apr,create_or_update_apr_clump
 use io_summary, only:iosum_apr,print_apr
 use eos,        only:gamma
 real,    intent(inout)         :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:)
 integer, intent(inout)         :: npart
 integer(kind=1), intent(inout) :: apr_level(:)
 integer :: ii,jj,kk,npartnew,nsplit_total,apri,npartold,ll,idx_len,j,apr_last
 integer :: n_ref,nrelax,nmerge,nkilled,nmerge_total,mm,n_to_split,iclosest,localtmp
 real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:)
 real, allocatable :: xyzh_merge(:,:),vxyzu_merge(:,:), rneighs(:)
 integer, allocatable :: relaxlist(:),mergelist(:),should_split(:)
 integer, allocatable :: idx_merge(:),should_merge(:),scan_array(:),idx_split(:)
 real :: get_apr_in(3),ientropy,P_i,pmassi,rhoi,xi,yi,zi,dx,dy,dz,rmin_local
 logical :: relax_in_loop

 ! if this routine doesn't need to be used, just skip it
 if (apr_max == 1) return

 if (npart >= 0.9*maxp) then
    call fatal('apr','maxp is not large enough; set --maxp on the command line to something larger than ',var='maxp',ival=maxp)
 endif
 ! if the centre of the region can move, update it
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)

 ! if we don't have any regions, skip routine
 if (ntrack == 0) return

 ! Just a metric
 if (apr_verbose) print*,'original npart is',npart

 ! initialise for the entropy storage
 if (allocated(entropy_list)) deallocate(entropy_list,entropy_stored)
 allocate(entropy_list(maxp*3),entropy_stored(maxp*3))
 entropy_count = 0
 entropy_list(:) = 0
 entropy_stored = 0.

 ! Before adjusting the particles, if we're going to
 ! relax them then let's save the reference particles
 if (do_relax) then
    allocate(xyzh_ref(4,maxp),force_ref(3,maxp),pmass_ref(maxp),relaxlist(maxp))
    relaxlist = -1

    n_ref = 0
    xyzh_ref = 0.
    force_ref = 0.
    pmass_ref = 0.

    do ii = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,ii))) then ! ignore dead particles
          n_ref = n_ref + 1
          xyzh_ref(1:4,n_ref) = xyzh(1:4,ii)
          pmass_ref(n_ref) = aprmassoftype(igas,apr_level(ii))
          force_ref(1:3,n_ref) = fxyzu(1:3,ii)*pmass_ref(n_ref)
       endif
    enddo
 else
    allocate(relaxlist(1))  ! it is passed but not used in merge
 endif

 ! Do any particles need to be split?
 npartnew = npart
 npartold = npart
 nsplit_total = 0
 nrelax = 0
 apri = 0 ! to avoid compiler errors
 apr_last = 0
 ! generally a safe guess, gets checked later
 allocate(scan_array(npart*apr_max),rneighs(npart*apr_max),idx_split(npart*apr_max),should_split(maxp))

 if (apr_verbose) print*,'started splitting'

 do jj = 1,apr_max-1
    do ll = 1,ntrack ! for multiple regions
       icentre = ll
       npartold = npartnew ! to account for new particles as they are being made
       should_split(:) = 0 ! reset
       rneighs(:) = 0.
       idx_split(:) = 0
       n_to_split = 0

       !$omp parallel default(none) &
       !$omp shared(npartold,iphase,apr_level,xyzh,should_split,get_apr,icentre) &
       !$omp shared(idx_split) &
       !$omp private(ii,get_apr_in,apri) &
       !$omp reduction(+:nsplit_total,n_to_split) reduction(max:apr_last)
       !$omp do
       split_over_active: do ii = 1,npartold
          ! only do this on active particles
          if (ind_timesteps) then
             if (.not.iactive(iphase(ii))) cycle split_over_active
          endif

          get_apr_in(1:3) = xyzh(1:3,ii)
          ! this is the refinement level it *should* have based
          ! on it's current position
          call get_apr(get_apr_in,icentre,apri)
          ! if the level it should have is greater than the
          ! level it does have, increment it up one
          if (apri > apr_level(ii)) then
             should_split(ii) = 1 ! record that this should be split
             nsplit_total = nsplit_total + 1
             n_to_split = n_to_split + 1
             apr_last = apri
          endif
       enddo split_over_active
       !$omp enddo
       !$omp end parallel

       ! reallocate if required; if this happens even once just use the biggest possible
       if (n_to_split > size(scan_array)) then
          deallocate(scan_array,rneighs,idx_split)
          allocate(scan_array(maxp),rneighs(maxp),idx_split(maxp))
       endif

       ! create the scan array - this loop should *not* be parallelised
       scan_array(:) = 0
       do ii = 2,npartold
          scan_array(ii) = scan_array(ii-1) + should_split(ii-1)
       enddo

       ! make the particle list
       idx_len = n_to_split
       npartnew = npartnew + idx_len ! total number of particles (for now)
       npartoftype(igas) = npartoftype(igas) + n_to_split ! add to npartoftype
       npart = npartnew ! for splitpart

       ! exit here if there's nothing more to do
       if (n_to_split == 0) cycle

       !$omp parallel default(none) &
       !$omp shared(npartold,should_split,idx_split,scan_array,idx_len) &
       !$omp shared(rneighs,xyzh,adjusted_split) &
       !$omp private(ii,mm,rmin_local,j,xi,yi,zi,dx,dy,dz)
       !$omp do
       do ii = 1,npartold
          if (should_split(ii) == 1) then
             idx_split(scan_array(ii) + 1) = ii
          endif
       enddo
       !$omp enddo

       if (adjusted_split) then
          !$omp do schedule(dynamic)
          do ii = 1,idx_len
             mm = idx_split(ii) ! original particle that should be split
             xi = xyzh(1,mm)
             yi = xyzh(2,mm)
             zi = xyzh(3,mm)

             rmin_local = huge(1.0)

             do j = 1,npartold
                if (j == mm) cycle
                dx = xi - xyzh(1,j)
                dy = yi - xyzh(2,j)
                dz = zi - xyzh(3,j)
                rmin_local = min(rmin_local,dx*dx + dy*dy + dz*dz)
             enddo
             rneighs(ii) = sqrt(rmin_local)
          enddo
          !$omp enddo
       endif
       !$omp end parallel

       ! if relaxing, make some adjustments here:
       ! just use the first particle that has been marked to split
       ! to establish if we should be relaxing at all
       relax_in_loop = (do_relax .and. (gr .or. apr_last == top_level))

       ! now go through and actually split them - this should *probably* not be parallelised
       ! due to the content of the nested functions, idx_len probably isn't that long either
       do ii = 1,idx_len
          mm = idx_split(ii) ! original particle that should be split
          kk = npartold + ii ! location in array for new particle
          if (adjusted_split) then
             call splitpart(mm,kk,rneigh=rneighs(ii))
          else
             call splitpart(mm,kk)
          endif
          if (relax_in_loop) then
             relaxlist(nrelax + ii) = mm
             relaxlist(nrelax + n_to_split + ii) = kk
          endif
          pmassi = aprmassoftype(igas,apr_level(ii))
          P_i = eos_vars(igasP,ii)
          rhoi = rhoh(xyzh(4,ii),pmassi)
          ientropy = pmassi*(P_i*rhoi**(-gamma))
          entropy_count = entropy_count + 2
          entropy_stored(entropy_count - 1:entropy_count) = 0.5*ientropy ! because we share it across both evenly
          entropy_list(entropy_count - 1) = iorig(ii)
          entropy_list(entropy_count) = iorig(kk)
       enddo

       ! if relaxing, update the total number that will be relaxed
       if (relax_in_loop) nrelax = nrelax + 2*n_to_split
    enddo
 enddo

 ! Take into account all the added particles
 npart = npartnew
 ntot = npartnew
 if (apr_verbose) then
    print*,'split: ',nsplit_total
    print*,'npart: ',npart
 endif

 ! Do any particles need to be merged?
 deallocate(scan_array)
 allocate(mergelist(npart),xyzh_merge(4,npart),vxyzu_merge(maxvxyzu,npart))
 allocate(idx_merge(npart),should_merge(npart),scan_array(npart))
 npart_regions = 0
 nmerge_total = 0
 iclosest = 1
 do jj = 1,apr_max-1
    do ll = 1, ntrack
       icentre = ll
       kk = apr_max - jj + 1             ! to go from apr_max -> 2
       mergelist = -1 ! initialise
       nmerge = 0
       nkilled = 0
       xyzh_merge = 0.
       vxyzu_merge = 0.

       should_merge(:) = 0
       scan_array(:) = 0
       idx_merge(:) = 0

       ! identify what should be merged
       !$omp parallel do default(none) &
       !$omp shared(npart,apr_level,kk,xyzh,vxyzu,ntrack,ll,should_merge,iphase,apr_centre) &
       !$omp private(ii,iclosest) &
       !$omp reduction(+:nmerge)
       merge_over_active: do ii = 1,npart
          if ((apr_level(ii) == kk) .and. (.not.isdead_or_accreted(xyzh(4,ii)))) then ! avoid already dead particles
             if (ind_timesteps) then
                if (.not.iactive(iphase(ii))) cycle merge_over_active
             endif
             if (ntrack > 1) call find_closest_region(xyzh(1:3,ii),ntrack,apr_centre,iclosest)

             if ((ntrack == 1) .or. (iclosest == ll)) then
                should_merge(ii) = 1
                nmerge = nmerge + 1
             endif
          endif
       enddo merge_over_active
       !$omp end parallel do

       ! create the scan array - this loop should *not* be parallelised
       scan_array(:) = 0
       do ii = 2,npart
          scan_array(ii) = scan_array(ii-1) + should_merge(ii-1)
       enddo

       !$omp parallel do default(none) &
       !$omp shared(should_merge,idx_merge,scan_array,xyzh_merge,vxyzu_merge,kk) &
       !$omp shared(xyzh,vxyzu,npart) &
       !$omp private(ii,mm) &
       !$omp reduction(+:npart_regions)
       do ii = 1,npart
          if (should_merge(ii) == 1) then
             mm = scan_array(ii) + 1
             idx_merge(mm) = ii
             xyzh_merge(1:4,mm) = xyzh(1:4,ii)
             vxyzu_merge(1:3,mm) = vxyzu(1:3,ii)
             npart_regions(kk) = npart_regions(kk) + 1
          endif
       enddo
       !$omp end parallel do

       if (apr_verbose) print*,nmerge,'particles selected for merge'
       ! Now send them to be merged
       if (nmerge > 11) call merge_with_special_tree(nmerge,idx_merge,xyzh_merge(:,1:nmerge),&
                                            vxyzu_merge(:,1:nmerge),kk,xyzh,vxyzu,apr_level,nkilled,&
                                            nrelax,relaxlist,npartnew,entropy_list,entropy_count,entropy_stored)
       nmerge_total = nmerge_total + nkilled ! actually merged
       if (apr_verbose) then
          print*,'merged: ',nkilled,kk
          print*,'npart: ',npartnew - nkilled
       endif
       npart_regions(kk) = npart_regions(kk) - nkilled
    enddo
 enddo
 ! update npart as required
 npart = npartnew
 npart_regions(1) = npartnew - sum(npart_regions(2:apr_max))
 if (apr_verbose) print*,'particles at each level:',npart_regions(:)


 ! If we need to relax, do it here
 if (nrelax > 0 .and. do_relax) call relax_particles(npart,n_ref,xyzh_ref,force_ref,nrelax,relaxlist)
 ! Turn it off now because we only want to do this on first splits
 if (.not. gr) do_relax = .false.

 ! As we may have killed particles, time to do an array shuffle
 call shuffle_part(npart)

 ! Tidy up
 if (do_relax) then
    deallocate(xyzh_ref,force_ref,pmass_ref)
 endif
 deallocate(relaxlist,should_merge,idx_merge,scan_array,rneighs,idx_split,should_split)

 if (apr_verbose) print*,'total particles at end of apr: ',npart

 ! summary variables
 print_apr = .true.
 iosum_apr(1) = ntrack
 iosum_apr(2) = apr_max
 iosum_apr(3) = iosum_apr(3) + nsplit_total
 iosum_apr(4) = iosum_apr(4) + nmerge_total
 do ii = 1,apr_max
    iosum_apr(ii+4) = count(apr_level(1:npart) == ii)
 enddo

end subroutine update_apr

!-----------------------------------------------------------------------
!+
!  routine to split one particle into two
!+
!-----------------------------------------------------------------------
subroutine splitpart(i,i_new,rneigh)
 use part,         only:xyzh
 use physcon,      only:pi
 use vectorutils, only:cross_product3D,rotatevec
 use get_apr_level, only:split_dir_func
 use dim, only:ind_timesteps
 integer, intent(in) :: i,i_new
 real, optional :: rneigh
 real :: sep

 if (adjusted_split) then
    sep = min(sep_factor*xyzh(4,i),0.35*rneigh)
    sep = sep/xyzh(4,i)  ! for consistency later on
 else
    sep = sep_factor
 endif

 call split_dir_func(i,i_new,sep)

end subroutine splitpart

!-----------------------------------------------------------------------
!+
!  Take in all particles that *might* be merged at this apr_level
!  and use our special tree to merge what has left the region
!+
!-----------------------------------------------------------------------
subroutine merge_with_special_tree(nmerge,mergelist,xyzh_merge,vxyzu_merge,current_apr,&
                                     xyzh,vxyzu,apr_level,nkilled,nrelax,relaxlist,npartnew,&
                                     entropy_list,entropy_count,entropy_stored)
 use neighkdtree,   only:build_tree,ncells,leaf_is_active,get_cell_location
 use mpiforce,      only:cellforce
 use kdtree,        only:inodeparts,inoderange
 use part,          only:kill_particle,igas,igasP,igamma,eos_vars,rhoh
 use part,          only:combine_two_particles,aprmassoftype,iorig
 use dim,           only:ind_timesteps,maxvxyzu
 use get_apr_level, only:get_apr,put_in_smallest_bin
 use physcon,       only:pi
 use utils_apr,     only:apr_centre
 use vectorutils, only:cross_product3D,matrixinvert3D
 use eos,           only:gamma
 integer,         intent(inout) :: nmerge,nkilled,nrelax,relaxlist(:),npartnew,entropy_count
 integer(kind=8), intent(inout) :: entropy_list(:)
 integer(kind=1), intent(inout) :: apr_level(:)
 integer,         intent(in)    :: current_apr,mergelist(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),entropy_stored(:)
 real,            intent(inout) :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer :: remainder,icell,n_cell,apri,m,i,ierr,k,already_stored,localtmp
 integer :: eldest,tuther,testp,testpp,n,child_list(12),parent_list(6)
 integer,         allocatable :: apri_at_cells_com(:)
 real,            allocatable :: cells_com(:,:)
 real    :: com(3),pmassi,xyzh_fromicentre(3)
 real    :: r_ave,phi_ave,theta_ave,r_part,phi_part,ekin
 real    :: pos_com(3),vel_com(3),am(3),ogen,ogam(3),am_term(3)
 real    :: Q(3,3),pdash,qdash,det,phi,lamb(3),es(3,3),sum_temp,s_min,S,dist(3),inv_iner(3,3)
 real    :: test_a,test_b,test_c,vec_a(3),vec_b(3),vec_c(3),u(3),v(3),w(3)
 real    :: lm(3),iner(3,3),lm_ave(3),term(3),omega(3),delta_ekin
 real    :: alpha,alpha1,alpha2,discriminant,A,B,C,un(3)
 real    :: ientropy_tuther, rho_eldest, rho_tuther, P_eldest, P_tuther, gammai
 logical :: spherical
 type(cellforce)        :: cell

 ! First ensure that we're only sending in groups of 12 to the tree
 remainder = modulo(nmerge,12)
 nmerge = nmerge - remainder

 call build_tree(nmerge,nmerge,xyzh_merge(:,1:nmerge),vxyzu_merge(:,1:nmerge),&
                      for_apr=.true.)

 allocate(cells_com(3,ncells),apri_at_cells_com(ncells))

 ! get the center of the cell
 !$omp parallel do default(none) &
 !$omp shared(ncells,leaf_is_active,inoderange,inodeparts,spherical) &
 !$omp shared(xyzh_merge,apr_centre,icentre,cells_com) &
 !$omp private(icell,n_cell,com,m,i) &
 !$omp private(cell,r_ave,theta_ave,phi_ave,r_part,phi_part,xyzh_fromicentre)
 over_cells_part0: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells_part0 !--skip empty cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1

    com = 0.
    if (.not.spherical) then
       ! if not using spherical coordinates to check the cell location, just use existing info
       call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
       com(1:3) = cell%xpos(1:3)
    else
       ! if spherical chosen, calculated the com in spherical coordinates and check
       ! if that is within the boundary or not (convert back to cartesian com later on)
       r_ave = 0.
       theta_ave = 0.
       phi_ave = 0.
       ! spherically average the position of the particles around the current APR region
       do m = 1,n_cell
          i = inodeparts(inoderange(1,icell) + m - 1)
          xyzh_fromicentre(1:3) = xyzh_merge(1:3,i) - apr_centre(1:3,icentre)
          !print*,i,xyzh_merge(1:3,i)
          r_part = sqrt(dot_product(xyzh_fromicentre(1:3),xyzh_fromicentre(1:3)))
          r_ave = r_ave + r_part
          theta_ave = theta_ave + acos(xyzh_fromicentre(3)/r_part)
          phi_part = atan2(xyzh_fromicentre(2),xyzh_fromicentre(1))
          !if (phi_ave < 0.) phi_ave = phi_ave + 2.*pi
          phi_ave = phi_ave + phi_part
       enddo
       r_ave = r_ave/real(n_cell)
       theta_ave = theta_ave/real(n_cell)
       phi_ave = phi_ave/real(n_cell)

       ! now convert back to cartesian equivalents
       com(1) = r_ave*sin(theta_ave)*cos(phi_ave)
       com(2) = r_ave*sin(theta_ave)*sin(phi_ave)
       com(3) = r_ave*cos(theta_ave)
       com(:) = com(:) + apr_centre(1:3,icentre) ! for sending back into get_apr
    endif
    cells_com(:,icell) = com
 enddo over_cells_part0
 !$omp end parallel do
 
 ! not sure how to parallelize this, so I am just gonna run it separately
 over_cells_part1: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells_part1 !--skip empty cells
    call get_apr(cells_com(1:3,icell),icentre,apri)
    apri_at_cells_com(i) = apri
 enddo over_cells_part1

 ! Now use the centre of mass of each cell to check whether it should
 ! be merged or not
 spherical = .true.
 !$omp parallel do default(none) &
 !$omp shared(xyzh,vxyzu,iorig,ncells,leaf_is_active,inoderange,inodeparts,spherical) &
 !$omp shared(cells_com,apri_at_cells_com,do_relax,nrelax,relaxlist) &
 !$omp shared(apr_centre,current_apr,aprmassoftype,mergelist,eos_vars,gamma) &
 !$omp shared(apr_level,xyzh_merge,vxyzu_merge,entropy_count,entropy_list,entropy_stored) &
 !$omp private(icell,n_cell,i,m,u,v,w,vec_a,vec_b,vec_c,test_a,test_b,test_c,testp,testpp,ierr) &
 !$omp private(pos_com,vel_com,am,am_term,lm,lm_ave,ekin,delta_ekin,dist,child_list) &
 !$omp private(apri,pmassi,ogen,ogam,Q,pdash,qdash,det,phi,lamb,es,un,iner,inv_iner,omega) &
 !$omp private(r_part,sum_temp,s_min,S,gammai,parent_list,already_stored,localtmp,term) &
 !$omp private(A,B,C,discriminant,alpha,alpha1,alpha2) &
 !$omp private(eldest,rho_eldest,P_eldest) &
 !$omp private(tuther,rho_tuther,P_tuther,ientropy_tuther) &
 !$omp firstprivate(com) &
 !$omp reduction(+:nkilled)
 over_cells: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells !--skip empty cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1

    com = cells_com(:,icell)
    apri = apri_at_cells_com(icell)

    ! If the apr level based on the com is lower than the current level,
    ! we merge!
    if (apri < current_apr) then
       ! here we take 12 particles from each leaf in the tree and combine these into six new particles
       ! the new particles are constructed to conserve the average properties of the children

       pmassi = aprmassoftype(igas,apr_level(inodeparts(inoderange(1,icell)))) ! this *current* mass is correct
       ! because only particles to merge are sent in

       ! start by calculating (or using) the average properties of the 12 children
       pos_com = 0.
       vel_com(:) = 0.
       am(:) = 0.
       lm(:) = 0.
       i = inodeparts(inoderange(1,icell))
       ekin = 0.
       do m = 1,n_cell
          i = inodeparts(inoderange(1,icell) + m - 1)
          child_list(m) = i ! save these for later
          vel_com(:) = vel_com(:) + vxyzu_merge(1:3,i)
          pos_com(:) = pos_com(:) + xyzh_merge(1:3,i)
          ekin = ekin + 0.5*pmassi*(vxyzu_merge(1,i)**2 &
                 + vxyzu_merge(2,i)**2 + vxyzu_merge(3,i)**2)
          lm(:) = lm(:) + pmassi*vxyzu_merge(1:3,i)
       enddo

       vel_com(:) = vel_com(:)/real(n_cell)
       pos_com(:) = pos_com(:)/real(n_cell)
       ogen = ekin
       lm_ave(:) = lm(:)/(n_cell*pmassi)

       ! adjust the particle positions to the com frame
       do m = 1,n_cell
          i = inodeparts(inoderange(1,icell) + m - 1)
          xyzh_merge(1:3,i) = xyzh_merge(1:3,i) - pos_com(1:3)
       enddo

       ! calculate the quadrupole mass moment and the adjusted angular momentum
       Q = 0.
       do m = 1,n_cell
          i = child_list(m)
          r_part = dot_product(xyzh_merge(1:3,i),xyzh_merge(1:3,i))
          Q(1,1) = Q(1,1) + pmassi*(3.*xyzh_merge(1,i)**2 - r_part)
          Q(2,2) = Q(2,2) + pmassi*(3.*xyzh_merge(2,i)**2 - r_part)
          Q(3,3) = Q(3,3) + pmassi*(3.*xyzh_merge(3,i)**2 - r_part)
          Q(1,2) = Q(1,2) + pmassi*(3.*xyzh_merge(1,i)*xyzh_merge(2,i))
          Q(1,3) = Q(1,3) + pmassi*(3.*xyzh_merge(1,i)*xyzh_merge(3,i))
          Q(2,3) = Q(2,3) + pmassi*(3.*xyzh_merge(2,i)*xyzh_merge(3,i))
          call cross_product3D(xyzh_merge(1:3,i),vxyzu_merge(1:3,i),am_term(:))
          am(:) = am(:) + pmassi*am_term(:)
       enddo
       ! by definition
       Q(2,1) = Q(1,2)
       Q(3,1) = Q(1,3)
       Q(3,2) = Q(2,3)
       ogam(:) = am(:)

       ! calculate the terms we need for the eigenvectors
       pdash = -(Q(1,1)*Q(1,1) + Q(2,2)*Q(2,2) + Q(3,3)*Q(3,3) &
           + 2.*Q(1,2)*Q(1,2) + 2.*Q(1,3)*Q(1,3) + 2.*Q(2,3)*Q(2,3))/2.
       det =  Q(1,1)*(Q(2,2)*Q(3,3) - Q(2,3)*Q(3,2)) &
            - Q(1,2)*(Q(2,1)*Q(3,3) - Q(2,3)*Q(3,1)) &
            + Q(1,3)*(Q(2,1)*Q(3,2) - Q(2,2)*Q(3,1))
       qdash = -det
       phi = 1./3. * acos(3.*qdash * sqrt(-3/pdash)/(2*pdash))

       ! I think we don't need to sort these as we take the shortest later
       do m = 1,3
          lamb(m) = 2.*sqrt(-pdash/3) * cos(phi - ((2.*pi * (m - 1))/3))
       enddo

       ! now construct the eigenvectors
       do i = 1,3
          u = (/Q(1,1) - lamb(i), Q(1,2), Q(1,3)/)
          v = (/Q(2,1), Q(2,2) - lamb(i), Q(2,3)/)
          w = (/Q(3,1), Q(3,2), Q(3,3) - lamb(i)/)

          call cross_product3D(u,v,vec_a)
          call cross_product3D(u,w,vec_b)
          call cross_product3D(v,w,vec_c)

          test_a = dot_product(vec_a,vec_a)
          test_b = dot_product(vec_b,vec_b)
          test_c = dot_product(vec_c,vec_c)

          ! take the one that is longest (safest choice)
          if ((test_a > test_b) .and. (test_a > test_c)) then
             es(1:3,i) = vec_a/sqrt(test_a)
          else if ((test_b > test_a) .and. (test_b > test_c)) then
             es(1:3,i) = vec_b/sqrt(test_b)
          else
             es(1:3,i) = vec_c/sqrt(test_c)
          endif
       enddo

       ! calculate the distances for the anti-podal pairs
       sum_temp = 0.
       do m = 1,n_cell
          i = child_list(m)
          r_part = dot_product(xyzh_merge(1:3,i),xyzh_merge(1:3,i))
          sum_temp = sum_temp + 0.5/12. * r_part
       enddo

       s_min = 3.* abs(minval(lamb)) / (12.*pmassi)

       S = max(sum_temp, s_min*1.05)

       ! now calculate the distance vectors along these eigenvectors
       do i = 1,3
          dist(i) = sqrt((S/3.) + (lamb(i)/(12.*pmassi)))
       enddo

       ! merge the first six particles with the last six particles
       do m = 1,(n_cell/2)
          eldest = mergelist(inodeparts(inoderange(1,icell) + m - 1)) ! remember we're running off the mergelist
          tuther = mergelist(inodeparts(inoderange(1,icell) + m + 5)) ! + 5

          ! save the entropy - we need this saved for later
          rho_eldest = rhoh(xyzh(4,eldest),pmassi*0.5) ! I don't know why 0.5 is required here?!
          rho_tuther = rhoh(xyzh(4,tuther),pmassi*0.5)
          P_eldest = eos_vars(igasP,eldest)
          P_tuther = eos_vars(igasP,tuther)
          gammai = gamma
          ! check to see if this particle has already been merged and is on the list
          ientropy_tuther = 0.
          already_stored = -1
          !$omp atomic capture
          entropy_count = entropy_count + 1
          localtmp = entropy_count
          !$omp end atomic
          do k = 1, localtmp-1
             if (entropy_list(k) == iorig(eldest)) already_stored = k
             ! this is in case it's been merged before, it's about to be killed
             ! by setting it to -1, it shouldn't be identified in adjust_entropy routine
             if (entropy_list(k) == iorig(tuther)) then
                entropy_list(k) = -1
                ientropy_tuther = entropy_stored(k)
             end if
          enddo
          ! use stored ientropy when possible (instead of recomputing) to ensure entropy conservation
          if (ientropy_tuther == 0.) ientropy_tuther = ientropy_tuther + 0.5*pmassi*P_tuther*rho_tuther**(-gammai)
          if (already_stored < 0) then
             entropy_stored(localtmp) = 0.5*pmassi*P_eldest*rho_eldest**(-gammai) + ientropy_tuther
             entropy_list(localtmp) = iorig(eldest)
          else
             entropy_stored(already_stored) = entropy_stored(already_stored) + ientropy_tuther
             entropy_list(localtmp) = -1    ! date already stored in 'already_stored', so mark new space as ignored
          endif

          ! discard tuther ("the other")
          call combine_two_particles(eldest,tuther)
          parent_list(m) = eldest
          apr_level(eldest) = apr_level(eldest) - int(1,kind=1)
          xyzh(4,eldest) = (xyzh(4,eldest))*(2.0**(1./3.)) ! rescale for its new mass
          if (ind_timesteps) call put_in_smallest_bin(eldest)

          ! book-keeping
          localtmp = nrelax
          if (do_relax) then
             !$omp atomic capture
             nrelax = nrelax + 1
             localtmp = nrelax
             !$omp end atomic
             relaxlist(localtmp) = eldest
          endif

          ! If this particle was on the shuffle list previously, take it off
          do n = 1,localtmp
             if (relaxlist(n) == tuther) relaxlist(n) = 0
          enddo

       enddo

       nkilled = nkilled + 12 ! this refers to the number of children killed

       ! now adjust the particle positions accordingly - we adjust away from com later
       ! particle 1:
       testp = parent_list(1)
       xyzh(1:3,testp) = dist(1) * es(:,1)

       ! particle 2:
       testpp = parent_list(2)
       xyzh(1:3,testpp) = -dist(1) * es(:,1)

       ! particle 3:
       testp = parent_list(3)
       xyzh(1:3,testp) = dist(2) * es(:,2)

       ! particle 4:
       testpp = parent_list(4)
       xyzh(1:3,testpp) = -dist(2) * es(:,2)

       ! particle 5:
       testp = parent_list(5)
       xyzh(1:3,testp) = dist(3) * es(:,3)

       ! particle 6:
       testpp = parent_list(6)
       xyzh(1:3,testpp) = -dist(3) * es(:,3)

       ! now we set the velocities
       ! first calculate the inertia tensor of the 6 new particles
       iner(:,:) = 0.
       pmassi = 2.*pmassi
       do m = 1,6
          i = parent_list(m)
          iner(1,1) = iner(1,1) + pmassi * (xyzh(2,i)**2 + xyzh(3,i)**2)
          iner(2,2) = iner(2,2) + pmassi * (xyzh(1,i)**2 + xyzh(3,i)**2)
          iner(3,3) = iner(3,3) + pmassi * (xyzh(1,i)**2 + xyzh(2,i)**2)
          iner(1,2) = iner(1,2) - pmassi * (xyzh(1,i)*xyzh(2,i))
          iner(1,3) = iner(1,3) - pmassi * (xyzh(1,i)*xyzh(3,i))
          iner(2,3) = iner(2,3) - pmassi * (xyzh(2,i)*xyzh(3,i))
       enddo

       iner(2,1) = iner(1,2)
       iner(3,1) = iner(1,3)
       iner(3,2) = iner(2,3)

       ! now invert the matrix and set the individual velocities
       call matrixinvert3D(iner,inv_iner,ierr)
       omega(:) = matmul(inv_iner,am)
       do m = 1,6
          i = parent_list(m)
          call cross_product3D(omega,xyzh(1:3,i),term)
          vxyzu(1:3,i) = lm_ave(1:3) + term
       enddo


       ! and now account for kinetic energy:
       ! calculate the current kinetic energy
       ekin = 0.
       do m = 1,6
          i = parent_list(m)
          ekin = ekin + 0.5 * pmassi * dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       enddo
       ! and the difference between original and current is (what we need to match)
       delta_ekin = ekin - ogen

       ! we write this out as a quadratic and solve the quadratic formula
       B = 0.
       do m = 1,6
          i = parent_list(m)
          un = xyzh(1:3,i) / sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
          B = B + pmassi * dot_product(vxyzu(1:3,i),un(:))
       enddo
       C = delta_ekin
       A = 3.*pmassi

       ! now we solve the quadratic
       discriminant = B**2 - 4.*A*C
       if (discriminant > 0.) then
          alpha1 = (-B + sqrt(discriminant)) / (2.*A)
          alpha2 = (-B - sqrt(discriminant)) / (2.*A)

          ! and take the solution that has the same sign as delta_ekin
          if (alpha1 * delta_ekin > 0.) then
             alpha = alpha1
          else
             alpha = alpha2
          endif
       else
          alpha = 0.
       endif

       ! finally, adjust the velocities accordingly and move the particles back to the simulation frame
       do m = 1,6
          i = parent_list(m)
          un = xyzh(1:3,i) / sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
          vxyzu(1:3,i) = vxyzu(1:3,i) + alpha*un
          xyzh(1:3,i) = xyzh(1:3,i) + pos_com(:)
       enddo

    endif

 enddo over_cells
 !$omp end parallel do

 deallocate(cells_com,apri_at_cells_com)

end subroutine merge_with_special_tree

end module apr
