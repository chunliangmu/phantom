!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module implicit
!
! Utility routines for implicit radiative diffusion
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: io, neighcompact
!
 implicit none
 integer              :: ncompact,ncompactlocal,icompactmax,nneigh_average
 real,    allocatable :: vari(:,:),varij(:,:),varij2(:,:),varinew(:,:)
 integer, allocatable :: ivar(:,:),ijvar(:)
 logical, allocatable :: mask(:)
 logical              :: done_allocation = .false.
 integer              :: npart_alloc = 0
 integer              :: icompact_alloc = 0

 public :: allocate_memory_implicit,grow_memory_implicit

 private :: deallocate_implicit

contains

!---------------------------------------------------------
!+
!  allocate arrays for compacted neighbour lists.
!
!  NB: this is called every radiation step and MUST grow with
!  npart: with APR the number of particles increases during a
!  calculation, and writing beyond the initial allocation would
!  silently corrupt memory.
!+
!---------------------------------------------------------
subroutine allocate_memory_implicit(npart,radkern,hfact,ierr)
 use io,           only:fatal
 use neighcompact, only:compact_neigh_counts
 integer, intent(in)  :: npart
 real,    intent(in)  :: radkern,hfact
 integer, intent(out) :: ierr
 integer :: npart_alloc_new,npart_alloc_old

 ierr = 0
 if (done_allocation .and. npart <= npart_alloc) return

 npart_alloc_old = npart_alloc
 if (done_allocation) call deallocate_implicit
 if (npart_alloc_old > 0) then
    ! grow with a margin: APR adds particles one at a time, and reallocating
    ! the whole neighbour list for every extra particle would be wasteful
    npart_alloc_new = max(npart,npart_alloc_old + npart_alloc_old/10 + 1000)
 else
    npart_alloc_new = npart
 endif

 call compact_neigh_counts(npart_alloc_new,radkern,hfact,nneigh_average,icompactmax)
 if (icompactmax <= 0) then
    call fatal('radiation_implicit','invalid compacted neighbour-list size',&
               var='icompactmax',ival=icompactmax)
 endif
 allocate(ivar(3,npart_alloc_new),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ivar')
 allocate(ijvar(icompactmax),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for ijvar')
 allocate(vari(2,npart_alloc_new),varij(2,icompactmax),varij2(4,icompactmax),&
          varinew(3,npart_alloc_new),mask(npart_alloc_new),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot allocate memory for vari, varij, varij2, varinew, mask')
 npart_alloc = npart_alloc_new
 icompact_alloc = icompactmax
 done_allocation = .true.

end subroutine allocate_memory_implicit

!---------------------------------------------------------
!+
!  grow the edge-indexed arrays if the (symmetric) neighbour
!  list needs more space than the initial estimate. This can
!  happen with APR, where the neighbour count of a coarse
!  particle surrounded by refined particles is larger than the
!  uniform-density estimate used in compact_neigh_counts.
!+
!---------------------------------------------------------
subroutine grow_memory_implicit(nedges,ierr)
 use io, only:fatal
 integer, intent(in)  :: nedges
 integer, intent(out) :: ierr

 ierr = 0
 if (nedges <= icompact_alloc) return
 if (.not.done_allocation) call fatal('radiation_implicit','grow_memory_implicit called before allocation')

 ! grow with a margin so that a slowly growing neighbour list is not
 ! reallocated on every step; integer arithmetic also avoids rounding
 ! the allocation back below nedges for large lists
 if (nedges > huge(icompactmax) - nedges/10 - 1) then
    call fatal('radiation_implicit','compacted neighbour-list size exceeds default integer range',&
               var='nedges',ival=nedges)
 endif
 icompactmax = nedges + nedges/10 + 1
 deallocate(ijvar,varij,varij2)
 allocate(ijvar(icompactmax),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot grow memory for ijvar')
 allocate(varij(2,icompactmax),varij2(4,icompactmax),stat=ierr)
 if (ierr/=0) call fatal('radiation_implicit','cannot grow memory for varij, varij2')
 icompact_alloc = icompactmax

end subroutine grow_memory_implicit

!---------------------------------------------------------
!+
!  free all arrays (only used internally, but kept separate
!  for clarity)
!+
!---------------------------------------------------------
subroutine deallocate_implicit()
 if (done_allocation) deallocate(ivar,ijvar,vari,varij,varij2,varinew,mask)
 done_allocation = .false.
 npart_alloc = 0
 icompact_alloc = 0
 icompactmax = 0

end subroutine deallocate_implicit

end module implicit
