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

 ! initialise the regions
 if (.not.allocated(apr_regions)) allocate(apr_regions(apr_max),npart_regions(apr_max))
 call set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)
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
                    shuffle_part,iphase,iactive,maxp,npartoftype,&
 use part,       only:igasP,rho,eos_vars,iorig
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
 integer :: n_ref,nrelax,nmerge,nkilled,nmerge_total,mm,n_to_split,iclosest
 real, allocatable :: xyzh_ref(:,:),force_ref(:,:),pmass_ref(:),rho_ref(:)
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
    allocate(xyzh_ref(4,maxp),force_ref(3,maxp),pmass_ref(maxp),rho_ref(maxp),relaxlist(maxp))
    relaxlist = -1

    n_ref = 0
    xyzh_ref = 0.
    force_ref = 0.
    pmass_ref = 0.
    rho_ref = 0.

    do ii = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,ii))) then ! ignore dead particles
          n_ref = n_ref + 1
          xyzh_ref(1:4,n_ref) = xyzh(1:4,ii)
          pmass_ref(n_ref) = aprmassoftype(igas,apr_level(ii))
          rho_ref(n_ref) = rho(ii)
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
          pmassi = aprmassoftype(igas,apr_level(mm))
          P_i = eos_vars(igasP,mm)
          rhoi = rho(mm)
          ientropy = pmassi*(P_i*rhoi**(-gamma))
          if (adjusted_split) then
             call splitpart(mm,kk,rneigh=rneighs(ii))
          else
             call splitpart(mm,kk)
          endif
          if (relax_in_loop) then
             relaxlist(nrelax + ii) = mm
             relaxlist(nrelax + n_to_split + ii) = kk
          endif
         !  entropy_count = entropy_count + 2
         !  entropy_stored(entropy_count - 1:entropy_count) = 0.5*ientropy
         !  ! because we share it across both evenly
         !  entropy_list(entropy_count - 1) = iorig(mm)
         !  entropy_list(entropy_count) = iorig(kk)
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
             if (maxvxyzu > 3) vxyzu_merge(4,mm) = vxyzu(4,ii)
             npart_regions(kk) = npart_regions(kk) + 1
          endif
       enddo
       !$omp end parallel do

       if (apr_verbose) print*,nmerge,'particles selected for merge'
       ! Now send them to be merged
       if (nmerge >= 4) call merge_with_special_tree(nmerge,idx_merge,xyzh_merge(:,1:nmerge),&
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
 if (nrelax > 0 .and. do_relax) call relax_particles(npart,n_ref,xyzh_ref,force_ref,rho_ref,nrelax,relaxlist)
 ! Turn it off now because we only want to do this on first splits
 if (.not. gr) do_relax = .false.

 ! As we may have killed particles, time to do an array shuffle
 call shuffle_part(npart)

 ! Tidy up
 if (do_relax) then
    deallocate(xyzh_ref,force_ref,pmass_ref,rho_ref)
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
 use part,          only:combine_two_particles
 use dim,           only:ind_timesteps
 use io,            only:fatal,warning
 use get_apr_level, only:get_apr,put_in_smallest_bin
 use sortutils,    only:indexx
 use vectorutils,   only:cross_product3D
 integer,         intent(inout) :: nmerge,nkilled,nrelax,relaxlist(:),npartnew,entropy_count
 integer(kind=8), intent(inout) :: entropy_list(:)
 integer(kind=1), intent(inout) :: apr_level(:)
 integer,         intent(in)    :: current_apr,mergelist(:)
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),entropy_stored(:)
 real,            intent(inout) :: xyzh_merge(:,:),vxyzu_merge(:,:)
 integer :: remainder,icell,n_cell,apri,m,i,j,k,n,localtmp,ia,ib,ic,id
 integer :: keep1,keep2,kill1,kill2,child_list(12),closest(4)
 integer,         allocatable :: apri_at_cells_com(:)
 real    :: am2,ekin,rl2,r_sep,vperp,vpar,sv,slen
 real    :: pos_com(3),vel_com(3),am(3),am_hat(3),am_term(3)
 real    :: vec(3),tvec(3),v_rel(3),r_rel(3)
 real    :: ppos(4,3),svec(3),svec_best(3),d2(4)
 type(cellforce)        :: cell

 ! First ensure that we're only sending in groups of 4 to the tree
 remainder = modulo(nmerge,4)
 nmerge = nmerge - remainder

 call build_tree(nmerge,nmerge,xyzh_merge(:,1:nmerge),vxyzu_merge(:,1:nmerge),&
                      for_apr=.true.)

 allocate(apri_at_cells_com(ncells))
 apri_at_cells_com = 0

 ! Get the apr level at the centre of mass of each leaf
 ! The group is merged once its com has crossed the boundary
 over_cells_part1: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells_part1 !--skip empty cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1
    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    pos_com = cell%xpos
    call get_apr(pos_com,icentre,apri)
    apri_at_cells_com(icell) = apri
    if (apri >= current_apr) cycle over_cells_part1

    ! make sure all particles-to-merge have crossed the boundary
    ! to make sure new children particles do NOT spawn in the finer side and immediately got split again
    do m = 1,n_cell
       i = inodeparts(inoderange(1,icell) + m - 1)
       call get_apr(xyzh_merge(1:3,i),icentre,apri)
       if (apri_at_cells_com(icell) < apri) apri_at_cells_com(icell) = apri
    enddo
 enddo over_cells_part1

 ! Now use the centre of mass of each cell to check whether it should
 ! be merged or not
 !$omp parallel do default(none) schedule(dynamic) &
 !$omp shared(xyzh,vxyzu,ncells,leaf_is_active,inoderange,inodeparts) &
 !$omp shared(apri_at_cells_com,do_relax,nrelax,relaxlist) &
 !$omp shared(current_apr,mergelist,apr_level) &
 !$omp shared(xyzh_merge,vxyzu_merge,get_apr,icentre) &
 !$omp private(icell,n_cell,i,m,j,k,n,apri,localtmp,ia,ib,ic,id) &
 !$omp private(keep1,keep2,kill1,kill2,child_list,closest) &
 !$omp private(pos_com,vel_com,am,am_hat,am_term,am2,ekin) &
 !$omp private(vec,tvec,v_rel,r_rel) &
 !$omp private(ppos,svec,svec_best,d2,slen,rl2,r_sep,vperp,vpar,sv) &
 !$omp reduction(+:nkilled)
 over_cells: do icell=1,int(ncells)
    if (leaf_is_active(icell) == 0) cycle over_cells !--skip empty cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1

    apri = apri_at_cells_com(icell)

    ! If the apr level based on the com is lower than the current level,
    ! we merge!
    if (apri < current_apr) then
       ! here we take the 4 particles in each leaf of the tree and combine these into 2 new particles,
       ! the new particles are constructed to conserve the average properties of the children

       ! the kdtree is built so that every leaf holds exactly 4 particles;
       ! skip anything else (should not happen) rather than crash
       if (n_cell /= 4) call fatal('merge_with_special_tree','Unexpected n_cell: should be 4, received ',var='n_cell',ival=n_cell)

       !  pmassi = aprmassoftype(igas,apr_level(mergelist(inodeparts(inoderange(1,icell))))) ! this *current* mass is correct
       !  ! because only particles to merge are sent in

       ! start by calculating the average properties of the 4 children
       pos_com = 0.
       vel_com(:) = 0.
       do m = 1,n_cell
          i = inodeparts(inoderange(1,icell) + m - 1)
          child_list(m) = i ! save these for later
          vel_com(:) = vel_com(:) + vxyzu_merge(1:3,i)
          pos_com(:) = pos_com(:) + xyzh_merge(1:3,i)
       enddo
       vel_com(:) = vel_com(:)/real(n_cell)
       pos_com(:) = pos_com(:)/real(n_cell)

       ! everything from here on is in the centre of mass reference frame of the group:
       ! am = angular momentum about the com (per unit mass)
       ! ekin = kinetic energy in the com frame (per unit mass)
       am(:) = 0.
       ekin = 0.
       do m = 1,n_cell
          i = child_list(m)
          r_rel(1:3) = xyzh_merge(1:3,i) - pos_com(1:3)
          v_rel(1:3) = vxyzu_merge(1:3,i) - vel_com(1:3)
          call cross_product3D(r_rel,v_rel,am_term)
          am(:) = am(:) + am_term(:)
          ekin = ekin + 0.5*dot_product(v_rel,v_rel)
       enddo

       ! the 2 new particles must lie on a plane perpendicular to the total angular momentum of the group,
       ! so the pair axis is restricted to that plane;
       ! choose the axis and separation that
       ! minimise the sum of the squared distances between each new particle and the 2 parents it inherits from:
       !   D^2 = |r_1-p_a|^2 + |r_1-p_b|^2 + |r_2-p_c|^2 + |r_2-p_d|^2
       ! with r_2 = -r_1 in the com frame.
       ! For a fixed pairing {a,b},{c,d}:
       !   D^2 = 4*r^2 + sum_i|p_i|^2 - 4*r*|P(p_a+p_b)|
       ! where P projects onto the plane perpendicular to L,
       ! which is minimised at r = |P(p_a+p_b)|/2,
       ! i.e. each child sits at the projection of its own pair's centre of mass onto that plane.
       ! The best pairing is the one with the longest projected pair-sum.
       am2 = dot_product(am,am)
       if (am2 > 0.) then
          am_hat(:) = am(:)/sqrt(am2)
       endif

       ! parent positions in the com frame
       do m = 1,n_cell
          i = child_list(m)
          ppos(m,1:3) = xyzh_merge(1:3,i) - pos_com(1:3)
       enddo

       ! try each of the 3 pairings of the 4 parents, keeping the one with
       ! the longest projected pair-sum
       slen = 0.
       svec_best(:) = 0.
       do n = 1,3
          select case(n)
          case(1)
             ia = 1; ib = 2; ic = 3; id = 4
          case(2)
             ia = 1; ib = 3; ic = 2; id = 4
          case(3)
             ia = 1; ib = 4; ic = 2; id = 3
          end select
          svec(1:3) = ppos(ia,1:3) + ppos(ib,1:3)
          if (am2 > 0.) svec(:) = svec(:) - dot_product(svec,am_hat)*am_hat(:)
          if (dot_product(svec,svec) > slen) then
             slen = dot_product(svec,svec)
             svec_best(1:3) = svec(1:3)
          endif
       enddo

       ! separation of the pair: half the projected pair-sum
       ! (children at the projection of their parents' centre of mass),
       ! but not smaller than the value required
       ! for the tangential velocities (implied by angular momentum conservation) to be real
       rl2 = 0.
       if (ekin > tiny(ekin)) rl2 = am2/(8.*ekin)
       if (slen > tiny(slen)) then
          vec(1:3) = svec_best(1:3)/sqrt(slen)
          r_sep = max(0.5*sqrt(slen),sqrt(rl2))
       else
          call warning('merge_with_special_tree','all 4 particles-to-merge seem to be lumped together')
          vec(1:3) = 0.
          r_sep = 0.
       endif


       ! do not merge if the resulting pair would straddle the apr boundary,
       ! as the outside member would simply be re-split on the next step (merge -> re-split churn);
       ! the group is left until it has fully crossed over
       call get_apr(pos_com(1:3) + r_sep*vec(1:3),icentre,apri)
       if (apri >= current_apr) cycle over_cells
       call get_apr(pos_com(1:3) - r_sep*vec(1:3),icentre,apri)
       if (apri >= current_apr) cycle over_cells

       ! relative velocity of the pair:
       ! the tangential component is fixed by angular momentum conservation (|am| = 4*r*v_tangential per unit mass),
       ! and the radial component by kinetic energy conservation,
       ! with its sign chosen to match the mean radial motion of the old group
       vperp = 0.
       if (r_sep > tiny(r_sep)) vperp = sqrt(am2)/(4.*r_sep)
       vpar = sqrt(max(0.,0.5*ekin - vperp*vperp))
       sv = 0.
       do m = 1,n_cell
          i = child_list(m)
          sv = sv + dot_product(vxyzu_merge(1:3,i) - vel_com(1:3),vec(1:3))
       enddo
       if (sv < 0.) vpar = -vpar
       v_rel(:) = vpar*vec(:)
       if (am2 > 0.) then
          call cross_product3D(am_hat,vec,tvec) ! tangential direction in the plane
          v_rel(:) = v_rel(:) + vperp*tvec(:)
       endif

       ! the pair configuration must be finite (catches NaN in the input)
       if (isnan(r_sep) .or. isnan(vperp) .or. isnan(vpar) .or. &
           isnan(vec(1)) .or. isnan(vec(2)) .or. isnan(vec(3))) &
          call fatal('merge_with_special_tree','non-finite pair configuration in merge')

       ! each child inherits the properties of the 2 parents closest to it
       do m = 1,n_cell
          i = child_list(m)
          d2(m) = dot_product(ppos(m,1:3) - r_sep*vec(1:3),ppos(m,1:3) - r_sep*vec(1:3))
       enddo
       call indexx(n_cell,d2(1:n_cell),closest(1:n_cell))

       keep1 = mergelist(inodeparts(inoderange(1,icell) + closest(1) - 1))
       kill1 = mergelist(inodeparts(inoderange(1,icell) + closest(2) - 1))
       keep2 = mergelist(inodeparts(inoderange(1,icell) + closest(3) - 1))
       kill2 = mergelist(inodeparts(inoderange(1,icell) + closest(4) - 1))
       call combine_two_particles(keep1,kill1)
       call combine_two_particles(keep2,kill2)

       ! now set the new positions and velocities of the pair,
       ! which conserve the total mass, linear momentum, angular momentum and kinetic energy
       ! of the 4 children by construction
       xyzh(1:3,keep1) = pos_com(1:3) + r_sep*vec(1:3)
       xyzh(1:3,keep2) = pos_com(1:3) - r_sep*vec(1:3)
       vxyzu(1:3,keep1) = vel_com(1:3) + v_rel(1:3)
       vxyzu(1:3,keep2) = vel_com(1:3) - v_rel(1:3)

       ! rescale smoothing length for the new particle mass
       xyzh(4,keep1) = xyzh(4,keep1)*(2.0**(1./3.))
       xyzh(4,keep2) = xyzh(4,keep2)*(2.0**(1./3.))

       apr_level(keep1) = apr_level(keep1) - int(1,kind=1)
       apr_level(keep2) = apr_level(keep2) - int(1,kind=1)

       if (ind_timesteps) then
          call put_in_smallest_bin(keep1)
          call put_in_smallest_bin(keep2)
       endif

       ! book-keeping
       localtmp = nrelax
       if (do_relax) then
          !$omp critical
          ! use critical instead of atomic capture here to ensure relaxlist is fully written before the check loop next
          nrelax = nrelax + 2
          localtmp = nrelax
          relaxlist(localtmp-1) = keep1
          relaxlist(localtmp) = keep2
          !$omp end critical
       endif

       ! If these particles were on the shuffle list previously, take them off
       do n = 1,localtmp
          if (relaxlist(n) == kill1 .or. relaxlist(n) == kill2) relaxlist(n) = 0
       enddo

       nkilled = nkilled + 4 ! this refers to the number of children killed

    endif

 enddo over_cells
 !$omp end parallel do

 deallocate(apri_at_cells_com)

end subroutine merge_with_special_tree

end module apr
