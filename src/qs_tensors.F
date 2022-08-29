!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************
MODULE qs_tensors
   USE OMP_LIB,                         ONLY: omp_get_num_threads,&
                                              omp_get_thread_num
   USE ai_contraction,                  ONLY: block_add
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              libxsmm_abc_contract
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE,&
                                              cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_split
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_filter, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_matrix_type, &
        dbcsr_has_symmetry, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE dbt_api,                         ONLY: &
        dbt_blk_sizes, dbt_clear, dbt_copy, dbt_create, dbt_destroy, dbt_distribution_destroy, &
        dbt_distribution_new, dbt_distribution_type, dbt_filter, dbt_get_block, dbt_get_info, &
        dbt_get_num_blocks, dbt_get_nze_total, dbt_get_stored_coordinates, &
        dbt_iterator_next_block, dbt_iterator_num_blocks, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_ndims, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, &
        dbt_put_block, dbt_reserve_blocks, dbt_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                              distribute_molecules_2d
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                              hfx_add_single_cache_element,&
                                              hfx_decompress_first_cache,&
                                              hfx_flush_last_cache,&
                                              hfx_get_mult_cache_elements,&
                                              hfx_get_single_cache_element,&
                                              hfx_reset_cache_and_container
   USE hfx_types,                       ONLY: alloc_containers,&
                                              dealloc_containers,&
                                              hfx_cache_type,&
                                              hfx_compression_type,&
                                              hfx_container_type,&
                                              hfx_init_container
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor,&
                                              eri_2center,&
                                              eri_2center_derivs,&
                                              eri_3center,&
                                              eri_3center_derivs,&
                                              libint_potential_type
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_2eri, cp_libint_cleanup_2eri1, cp_libint_cleanup_3eri, &
        cp_libint_cleanup_3eri1, cp_libint_init_2eri, cp_libint_init_2eri1, cp_libint_init_3eri, &
        cp_libint_init_3eri1, cp_libint_set_contrdepth, cp_libint_t
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE mp2_eri_gpw,                     ONLY: cleanup_gpw,&
                                              mp2_eri_3c_integrate_gpw,&
                                              prepare_gpw
   USE mp2_types,                       ONLY: mp2_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: &
        get_iterator_info, get_neighbor_list_set_p, neighbor_list_iterate, &
        neighbor_list_iterator_create, neighbor_list_iterator_p_type, &
        neighbor_list_iterator_release, neighbor_list_set_p_type, nl_sub_iterate, &
        release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_tensors_types,                ONLY: &
        distribution_3d_destroy, distribution_3d_type, neighbor_list_3c_iterator_type, &
        neighbor_list_3c_type, symmetric_ij, symmetric_ijk, symmetric_jk, symmetric_none, &
        symmetrik_ik
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: get_limit

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tensors'

   PUBLIC :: build_3c_neighbor_lists, &
             neighbor_list_3c_destroy, neighbor_list_3c_iterate, neighbor_list_3c_iterator_create, &
             neighbor_list_3c_iterator_destroy, get_3c_iterator_info, build_3c_integrals, &
             build_2c_neighbor_lists, build_2c_integrals, cutoff_screen_factor, &
             get_tensor_occupancy, compress_tensor, decompress_tensor, &
             build_3c_derivatives, build_2c_derivatives, build_3c_integrals_gpw

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE    :: array
   END TYPE

   ! cache size for integral compression
   INTEGER, PARAMETER, PRIVATE :: cache_size = 1024

CONTAINS

! **************************************************************************************************
!> \brief Build 2-center neighborlists adapted to different operators
!>        This mainly wraps build_neighbor_lists for consistency with build_3c_neighbor_lists
!> \param ij_list 2c neighbor list for atom pairs i, j
!> \param basis_i basis object for atoms i
!> \param basis_j basis object for atoms j
!> \param potential_parameter ...
!> \param name name of 2c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .TRUE.)
!> \param molecular ...
!> \param dist_2d optionally a custom 2d distribution
!> \param pot_to_rad which radius (1 for i, 2 for j) should be adapted to incorporate potential
! **************************************************************************************************
   SUBROUTINE build_2c_neighbor_lists(ij_list, basis_i, basis_j, potential_parameter, name, qs_env, &
                                      sym_ij, molecular, dist_2d, pot_to_rad)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ij_list
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, molecular
      TYPE(distribution_2d_type), OPTIONAL, POINTER      :: dist_2d
      INTEGER, INTENT(IN), OPTIONAL                      :: pot_to_rad

      INTEGER                                            :: ikind, nkind, pot_to_rad_prv
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: i_present, j_present
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: i_radius, j_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(distribution_2d_type), POINTER                :: dist_2d_prv
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, cell, local_particles, molecule_set, &
               particle_set, dist_2d_prv)

      IF (PRESENT(pot_to_rad)) THEN
         pot_to_rad_prv = pot_to_rad
      ELSE
         pot_to_rad_prv = 1
      END IF

      CALL get_qs_env(qs_env, &
                      nkind=nkind, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      local_particles=local_particles, &
                      distribution_2d=dist_2d_prv, &
                      molecule_set=molecule_set)

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ALLOCATE (i_present(nkind), j_present(nkind))
      ALLOCATE (i_radius(nkind), j_radius(nkind))

      i_present = .FALSE.
      j_present = .FALSE.
      i_radius = 0.0_dp
      j_radius = 0.0_dp

      IF (PRESENT(dist_2d)) dist_2d_prv => dist_2d

      !  Set up the radii, depending on the operator type
      IF (potential_parameter%potential_type == do_potential_id) THEN

         !overlap => use the kind radius for both i and j
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=i_radius(ikind))
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_j(ikind)%gto_basis_set, kind_radius=j_radius(ikind))
            END IF
         END DO

      ELSE IF (potential_parameter%potential_type == do_potential_coulomb) THEN

         !Coulomb operator, virtually infinite range => set j_radius to arbitrarily large number
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               IF (pot_to_rad_prv == 1) i_radius(ikind) = 1000000.0_dp
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               IF (pot_to_rad_prv == 2) j_radius(ikind) = 1000000.0_dp
            END IF
         END DO !ikind

      ELSE IF (potential_parameter%potential_type == do_potential_truncated .OR. &
               potential_parameter%potential_type == do_potential_short) THEN

         !Truncated coulomb/short range: set j_radius to r_cutoff + the kind_radii
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=i_radius(ikind))
               IF (pot_to_rad_prv == 1) i_radius(ikind) = i_radius(ikind) + cutoff_screen_factor*potential_parameter%cutoff_radius
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_j(ikind)%gto_basis_set, kind_radius=j_radius(ikind))
               IF (pot_to_rad_prv == 2) j_radius(ikind) = j_radius(ikind) + cutoff_screen_factor*potential_parameter%cutoff_radius
            END IF
         END DO

      ELSE
         CPABORT("Operator not implemented.")
      END IF

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(i_present, j_present, i_radius, j_radius, pair_radius)

      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, local_particles, dist_2d_prv, atomic_kind_set, &
                        molecule_set, molecule_only=.FALSE., particle_set=particle_set)
      CALL build_neighbor_lists(ij_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                symmetric=sym_ij, molecular=molecular, nlname=TRIM(name))

      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Build a 3-center neighbor list
!> \param ijk_list 3c neighbor list for atom triples i, j, k
!> \param basis_i basis object for atoms i
!> \param basis_j basis object for atoms j
!> \param basis_k basis object for atoms k
!> \param dist_3d 3d distribution object
!> \param potential_parameter ...
!> \param name name of 3c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .FALSE.)
!> \param sym_jk Symmetry in j, k (default .FALSE.)
!> \param sym_ik Symmetry in i, k (default .FALSE.)
!> \param molecular ??? not tested
!> \param op_pos ...
!> \param own_dist ...
! **************************************************************************************************
   SUBROUTINE build_3c_neighbor_lists(ijk_list, basis_i, basis_j, basis_k, &
                                      dist_3d, potential_parameter, name, qs_env, &
                                      sym_ij, sym_jk, sym_ik, molecular, op_pos, own_dist)
      TYPE(neighbor_list_3c_type), INTENT(OUT)           :: ijk_list
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(distribution_3d_type), INTENT(IN)             :: dist_3d
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, sym_jk, sym_ik, molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist

      CHARACTER(len=*), PARAMETER :: routineN = 'build_3c_neighbor_lists'

      INTEGER                                            :: handle, op_pos_prv, sym_level
      TYPE(libint_potential_type)                        :: pot_par_1, pot_par_2

      CALL timeset(routineN, handle)

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      END IF

      SELECT CASE (op_pos_prv)
      CASE (1)
         pot_par_1 = potential_parameter
         pot_par_2%potential_type = do_potential_id
      CASE (2)
         pot_par_2 = potential_parameter
         pot_par_1%potential_type = do_potential_id
      END SELECT

      CALL build_2c_neighbor_lists(ijk_list%ij_list, basis_i, basis_j, pot_par_1, TRIM(name)//"_sub_1", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, &
                                   dist_2d=dist_3d%dist_2d_1, pot_to_rad=1)

      CALL build_2c_neighbor_lists(ijk_list%jk_list, basis_j, basis_k, pot_par_2, TRIM(name)//"_sub_2", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, &
                                   dist_2d=dist_3d%dist_2d_2, pot_to_rad=2)

      ijk_list%sym = symmetric_none

      sym_level = 0
      IF (PRESENT(sym_ij)) THEN
         IF (sym_ij) THEN
            ijk_list%sym = symmetric_ij
            sym_level = sym_level + 1
         END IF
      END IF

      IF (PRESENT(sym_jk)) THEN
         IF (sym_jk) THEN
            ijk_list%sym = symmetric_jk
            sym_level = sym_level + 1
         END IF
      END IF

      IF (PRESENT(sym_ik)) THEN
         IF (sym_ik) THEN
            ijk_list%sym = symmetrik_ik
            sym_level = sym_level + 1
         END IF
      END IF

      IF (sym_level >= 2) THEN
         ijk_list%sym = symmetric_ijk
      END IF

      ijk_list%dist_3d = dist_3d
      IF (PRESENT(own_dist)) THEN
         ijk_list%owns_dist = own_dist
      ELSE
         ijk_list%owns_dist = .FALSE.
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Symmetry criterion
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION include_symmetric(a, b)
      INTEGER, INTENT(IN)                                :: a, b
      LOGICAL                                            :: include_symmetric

      IF (a > b) THEN
         include_symmetric = (MODULO(a + b, 2) /= 0)
      ELSE
         include_symmetric = (MODULO(a + b, 2) == 0)
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief Destroy 3c neighborlist
!> \param ijk_list ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_destroy(ijk_list)
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: ijk_list

      CALL release_neighbor_list_sets(ijk_list%ij_list)
      CALL release_neighbor_list_sets(ijk_list%jk_list)

      IF (ijk_list%owns_dist) THEN
         CALL distribution_3d_destroy(ijk_list%dist_3d)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 3-center neighborlist iterator
!> \param iterator ...
!> \param ijk_nl ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_iterator_create(iterator, ijk_nl)
      TYPE(neighbor_list_3c_iterator_type), INTENT(OUT)  :: iterator
      TYPE(neighbor_list_3c_type), INTENT(IN)            :: ijk_nl

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_iterator_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL neighbor_list_iterator_create(iterator%iter_ij, ijk_nl%ij_list)
      CALL neighbor_list_iterator_create(iterator%iter_jk, ijk_nl%jk_list, search=.TRUE.)
      iterator%iter_level = 0
      iterator%ijk_nl = ijk_nl

      iterator%bounds_i = 0
      iterator%bounds_j = 0
      iterator%bounds_k = 0

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief impose atomic upper and lower bounds
!> \param iterator ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE nl_3c_iter_set_bounds(iterator, bounds_i, bounds_j, bounds_k)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      IF (PRESENT(bounds_i)) iterator%bounds_i = bounds_i
      IF (PRESENT(bounds_j)) iterator%bounds_j = bounds_j
      IF (PRESENT(bounds_k)) iterator%bounds_k = bounds_k

   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy 3c-nl iterator
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_iterator_destroy(iterator)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_iterator_destroy'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL neighbor_list_iterator_release(iterator%iter_ij)
      CALL neighbor_list_iterator_release(iterator%iter_jk)
      NULLIFY (iterator%iter_ij)
      NULLIFY (iterator%iter_jk)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Iterate 3c-nl iterator
!> \param iterator ...
!> \return 0 if successful; 1 if end was reached
! **************************************************************************************************
   RECURSIVE FUNCTION neighbor_list_3c_iterate(iterator) RESULT(iter_stat)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER                                            :: iter_stat

      INTEGER                                            :: iatom, iter_level, jatom, jatom_1, &
                                                            jatom_2, katom
      LOGICAL                                            :: skip_this

      iter_level = iterator%iter_level

      IF (iter_level == 0) THEN
         iter_stat = neighbor_list_iterate(iterator%iter_ij)

         IF (iter_stat /= 0) THEN
            RETURN
         END IF

         CALL get_iterator_info(iterator%iter_ij, iatom=iatom, jatom=jatom)
         skip_this = .FALSE.
         IF ((iterator%bounds_i(1) > 0 .AND. iatom < iterator%bounds_i(1)) &
             .OR. (iterator%bounds_i(2) > 0 .AND. iatom > iterator%bounds_i(2))) skip_this = .TRUE.
         IF ((iterator%bounds_j(1) > 0 .AND. jatom < iterator%bounds_j(1)) &
             .OR. (iterator%bounds_j(2) > 0 .AND. jatom > iterator%bounds_j(2))) skip_this = .TRUE.

         IF (skip_this) THEN
            iter_stat = neighbor_list_3c_iterate(iterator)
            RETURN
         END IF

      END IF
      iter_stat = nl_sub_iterate(iterator%iter_jk, iterator%iter_ij)
      IF (iter_stat /= 0) THEN
         iterator%iter_level = 0
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      ELSE
         iterator%iter_level = 1
      END IF

      CPASSERT(iter_stat == 0)
      CPASSERT(iterator%iter_level == 1)
      CALL get_iterator_info(iterator%iter_ij, iatom=iatom, jatom=jatom_1)
      CALL get_iterator_info(iterator%iter_jk, iatom=jatom_2, jatom=katom)

      CPASSERT(jatom_1 == jatom_2)
      jatom = jatom_1

      skip_this = .FALSE.
      IF ((iterator%bounds_k(1) > 0 .AND. katom < iterator%bounds_k(1)) &
          .OR. (iterator%bounds_k(2) > 0 .AND. katom > iterator%bounds_k(2))) skip_this = .TRUE.

      IF (skip_this) THEN
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      END IF

      SELECT CASE (iterator%ijk_nl%sym)
      CASE (symmetric_none)
         skip_this = .FALSE.
      CASE (symmetric_ij)
         skip_this = .NOT. include_symmetric(iatom, jatom)
      CASE (symmetric_jk)
         skip_this = .NOT. include_symmetric(jatom, katom)
      CASE (symmetrik_ik)
         skip_this = .NOT. include_symmetric(iatom, katom)
      CASE (symmetric_ijk)
         skip_this = .NOT. include_symmetric(iatom, jatom) .OR. .NOT. include_symmetric(jatom, katom)
      CASE DEFAULT
         CPABORT("should not happen")
      END SELECT

      IF (skip_this) THEN
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief Get info of current iteration
!> \param iterator ...
!> \param ikind ...
!> \param jkind ...
!> \param kkind ...
!> \param nkind ...
!> \param iatom ...
!> \param jatom ...
!> \param katom ...
!> \param rij ...
!> \param rjk ...
!> \param rik ...
!> \param cell_j ...
!> \param cell_k ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE get_3c_iterator_info(iterator, ikind, jkind, kkind, nkind, iatom, jatom, katom, &
                                   rij, rjk, rik, cell_j, cell_k)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER, INTENT(OUT), OPTIONAL                     :: ikind, jkind, kkind, nkind, iatom, &
                                                            jatom, katom
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: rij, rjk, rik
      INTEGER, DIMENSION(3), INTENT(OUT), OPTIONAL       :: cell_j, cell_k

      INTEGER, DIMENSION(2)                              :: atoms_1, atoms_2, kinds_1, kinds_2
      INTEGER, DIMENSION(3)                              :: cell_1, cell_2
      REAL(KIND=dp), DIMENSION(3)                        :: r_1, r_2

      CPASSERT(iterator%iter_level == 1)

      CALL get_iterator_info(iterator%iter_ij, &
                             ikind=kinds_1(1), jkind=kinds_1(2), nkind=nkind, &
                             iatom=atoms_1(1), jatom=atoms_1(2), r=r_1, &
                             cell=cell_1)

      CALL get_iterator_info(iterator%iter_jk, &
                             ikind=kinds_2(1), jkind=kinds_2(2), &
                             iatom=atoms_2(1), jatom=atoms_2(2), r=r_2, &
                             cell=cell_2)

      IF (PRESENT(ikind)) ikind = kinds_1(1)
      IF (PRESENT(jkind)) jkind = kinds_1(2)
      IF (PRESENT(kkind)) kkind = kinds_2(2)
      IF (PRESENT(iatom)) iatom = atoms_1(1)
      IF (PRESENT(jatom)) jatom = atoms_1(2)
      IF (PRESENT(katom)) katom = atoms_2(2)

      IF (PRESENT(rij)) rij = r_1
      IF (PRESENT(rjk)) rjk = r_2
      IF (PRESENT(rik)) rik = r_1 + r_2

      IF (PRESENT(cell_j)) cell_j = cell_1
      IF (PRESENT(cell_k)) cell_k = cell_1 + cell_2

   END SUBROUTINE

! **************************************************************************************************
!> \brief Allocate blocks of a 3-center tensor based on neighborlist
!> \param t3c empty DBCSR tensor
!>            Should be of shape (1,1) if no kpoints are used and of shape (nimages, nimages)
!>            if k-points are used
!> \param nl_3c 3-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param qs_env ...
!> \param potential_parameter ...
!> \param op_pos ...
!> \param do_kpoints ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, potential_parameter, op_pos, &
                             do_kpoints, bounds_i, bounds_j, bounds_k)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t3c
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'alloc_block_3c'

      INTEGER                                            :: handle, i_img, iatom, ikind, j_img, &
                                                            jatom, jcell, jkind, katom, kcell, &
                                                            kkind, natom, nimg, op_ij, op_jk, &
                                                            op_pos_prv
      INTEGER(int_8)                                     :: a, b, nblk_per_thread
      INTEGER(int_8), ALLOCATABLE, DIMENSION(:, :)       :: nblk
      INTEGER, DIMENSION(3)                              :: cell_j, cell_k, kp_index_lbounds, &
                                                            kp_index_ubounds
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(one_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: alloc_i, alloc_j, alloc_k
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      NULLIFY (qs_kind_set, atomic_kind_set)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      END IF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control, kpoints=kpoints, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      END IF

      !Do a first loop over the nl and count the blocks present
      ALLOCATE (nblk(nimg, nimg))
      nblk(:, :) = 0

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bounds_k)

      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE
         ELSE
            jcell = 1; kcell = 1
         END IF

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=kind_radius_i)
         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, kind_radius=kind_radius_j)
         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, kind_radius=kind_radius_k)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         nblk(jcell, kcell) = nblk(jcell, kcell) + 1
      END DO
      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

      !Do a second loop over the nl to give block indices
      ALLOCATE (alloc_i(nimg, nimg))
      ALLOCATE (alloc_j(nimg, nimg))
      ALLOCATE (alloc_k(nimg, nimg))
      DO j_img = 1, nimg
         DO i_img = 1, nimg
            ALLOCATE (alloc_i(i_img, j_img)%array(nblk(i_img, j_img)))
            ALLOCATE (alloc_j(i_img, j_img)%array(nblk(i_img, j_img)))
            ALLOCATE (alloc_k(i_img, j_img)%array(nblk(i_img, j_img)))
         END DO
      END DO
      nblk(:, :) = 0

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bounds_k)

      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE
         ELSE
            jcell = 1; kcell = 1
         END IF

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=kind_radius_i)
         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, kind_radius=kind_radius_j)
         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, kind_radius=kind_radius_k)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         nblk(jcell, kcell) = nblk(jcell, kcell) + 1

         !Note: there may be repeated indices due to periodic images => dbt_reserve_blocks takes care of it
         alloc_i(jcell, kcell)%array(nblk(jcell, kcell)) = iatom
         alloc_j(jcell, kcell)%array(nblk(jcell, kcell)) = jatom
         alloc_k(jcell, kcell)%array(nblk(jcell, kcell)) = katom

      END DO
      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

!TODO: Parallelize creation of block list.
!$OMP PARALLEL DEFAULT(NONE) SHARED(t3c,nimg,nblk,alloc_i,alloc_j,alloc_k) &
!$OMP PRIVATE(i_img,j_img,nblk_per_thread,A,b)
      DO j_img = 1, nimg
         DO i_img = 1, nimg
            IF (ALLOCATED(alloc_i(i_img, j_img)%array)) THEN
               nblk_per_thread = nblk(i_img, j_img)/omp_get_num_threads() + 1
               a = omp_get_thread_num()*nblk_per_thread + 1
               b = MIN(a + nblk_per_thread, nblk(i_img, j_img))
               CALL dbt_reserve_blocks(t3c(i_img, j_img), &
                                       alloc_i(i_img, j_img)%array(a:b), &
                                       alloc_j(i_img, j_img)%array(a:b), &
                                       alloc_k(i_img, j_img)%array(a:b))
            END IF
         END DO
      END DO
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t3c ...
!> \param nl_3c ...
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param qs_env ...
!> \param potential_parameter ...
!> \param op_pos ...
!> \param do_kpoints ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE alloc_block_3c_old(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, potential_parameter, op_pos, &
                                 do_kpoints, bounds_i, bounds_j, bounds_k)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t3c
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_block_3c_old'

      INTEGER :: blk_cnt, handle, i, i_img, iatom, iblk, ikind, iproc, j_img, jatom, jcell, jkind, &
         katom, kcell, kkind, natom, nimg, op_ij, op_jk, op_pos_prv
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp
      INTEGER, DIMENSION(3)                              :: cell_j, cell_k, kp_index_lbounds, &
                                                            kp_index_ubounds
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, new_block
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(one_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: alloc_i, alloc_j, alloc_k
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      NULLIFY (qs_kind_set, atomic_kind_set)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      END IF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control, kpoints=kpoints, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      ALLOCATE (alloc_i(nimg, nimg))
      ALLOCATE (alloc_j(nimg, nimg))
      ALLOCATE (alloc_k(nimg, nimg))

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      END IF

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bounds_k)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE
         ELSE
            jcell = 1; kcell = 1
         END IF

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=kind_radius_i)
         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, kind_radius=kind_radius_j)
         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, kind_radius=kind_radius_k)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         ! tensor is not symmetric therefore need to allocate rows and columns in
         ! correspondence with neighborlist. Note that this only allocates half
         ! of the blocks (since neighborlist is symmetric). After filling the blocks,
         ! tensor will be added to its transposed

         ASSOCIATE (ai => alloc_i(jcell, kcell))
            ASSOCIATE (aj => alloc_j(jcell, kcell))
               ASSOCIATE (ak => alloc_k(jcell, kcell))

                  new_block = .TRUE.
                  IF (ALLOCATED(aj%array)) THEN
                     DO iblk = 1, SIZE(aj%array)
                        IF (ai%array(iblk) == iatom .AND. &
                            aj%array(iblk) == jatom .AND. &
                            ak%array(iblk) == katom) THEN
                           new_block = .FALSE.
                           EXIT
                        END IF
                     END DO
                  END IF
                  IF (.NOT. new_block) CYCLE

                  IF (ALLOCATED(ai%array)) THEN
                     blk_cnt = SIZE(ai%array)
                     ALLOCATE (tmp(blk_cnt))
                     tmp(:) = ai%array(:)
                     DEALLOCATE (ai%array)
                     ALLOCATE (ai%array(blk_cnt + 1))
                     ai%array(1:blk_cnt) = tmp(:)
                     ai%array(blk_cnt + 1) = iatom
                  ELSE
                     ALLOCATE (ai%array(1))
                     ai%array(1) = iatom
                  END IF

                  IF (ALLOCATED(aj%array)) THEN
                     tmp(:) = aj%array(:)
                     DEALLOCATE (aj%array)
                     ALLOCATE (aj%array(blk_cnt + 1))
                     aj%array(1:blk_cnt) = tmp(:)
                     aj%array(blk_cnt + 1) = jatom
                  ELSE
                     ALLOCATE (aj%array(1))
                     aj%array(1) = jatom
                  END IF

                  IF (ALLOCATED(ak%array)) THEN
                     tmp(:) = ak%array(:)
                     DEALLOCATE (ak%array)
                     ALLOCATE (ak%array(blk_cnt + 1))
                     ak%array(1:blk_cnt) = tmp(:)
                     ak%array(blk_cnt + 1) = katom
                  ELSE
                     ALLOCATE (ak%array(1))
                     ak%array(1) = katom
                  END IF

                  IF (ALLOCATED(tmp)) DEALLOCATE (tmp)
               END ASSOCIATE
            END ASSOCIATE
         END ASSOCIATE
      END DO

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

      DO i_img = 1, nimg
         DO j_img = 1, nimg
            IF (ALLOCATED(alloc_i(i_img, j_img)%array)) THEN
               DO i = 1, SIZE(alloc_i(i_img, j_img)%array)
                  CALL dbt_get_stored_coordinates(t3c(i_img, j_img), &
                                                  [alloc_i(i_img, j_img)%array(i), alloc_j(i_img, j_img)%array(i), &
                                                   alloc_k(i_img, j_img)%array(i)], &
                                                  iproc)
                  CPASSERT(iproc .EQ. para_env%mepos)
               END DO

               CALL dbt_reserve_blocks(t3c(i_img, j_img), &
                                       alloc_i(i_img, j_img)%array, &
                                       alloc_j(i_img, j_img)%array, &
                                       alloc_k(i_img, j_img)%array)
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Build 3-center derivative tensors
!> \param t3c_der_i empty DBCSR tensor which will contain the 1st center derivatives
!> \param t3c_der_k empty DBCSR tensor which will contain the 3rd center derivatives
!> \param filter_eps Filter threshold for tensor blocks
!> \param qs_env ...
!> \param nl_3c 3-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param potential_parameter ...
!> \param der_eps neglect integrals smaller than der_eps
!> \param op_pos operator position.
!>        1: calculate (i|jk) integrals,
!>        2: calculate (ij|k) integrals
!> \param do_kpoints ...
!> this routine requires that libint has been static initialised somewhere else
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE build_3c_derivatives(t3c_der_i, t3c_der_k, filter_eps, qs_env, &
                                   nl_3c, basis_i, basis_j, basis_k, &
                                   potential_parameter, &
                                   der_eps, &
                                   op_pos, do_kpoints, &
                                   bounds_i, bounds_j, bounds_k)

      TYPE(dbt_type), DIMENSION(:, :, :), INTENT(INOUT)  :: t3c_der_i, t3c_der_k
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: der_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_derivatives'

      INTEGER :: block_end_i, block_end_j, block_end_k, block_start_i, block_start_j, &
         block_start_k, egfi, handle, handle2, i, i_img, i_xyz, iatom, ibasis, ikind, ilist, imax, &
         iset, j_img, jatom, jcell, jkind, jset, katom, kcell, kkind, kset, m_max, max_ncoi, &
         max_ncoj, max_ncok, max_nset, max_nsgfi, max_nsgfj, max_nsgfk, maxli, maxlj, maxlk, &
         mepos, natom, nbasis, ncoi, ncoj, ncok, nimg, nseti, nsetj, nsetk, nthread, op_ij, op_jk, &
         op_pos_prv, sgfi, sgfj, sgfk, unit_id
      INTEGER, DIMENSION(2)                              :: bo
      INTEGER, DIMENSION(3)                              :: blk_size, cell_j, cell_k, &
                                                            kp_index_lbounds, kp_index_ubounds, sp
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, found, skip
      LOGICAL, DIMENSION(3)                              :: block_j_not_zero, block_k_not_zero, &
                                                            der_j_zero, der_k_zero
      REAL(dp), DIMENSION(3)                             :: der_ext_i, der_ext_j, der_ext_k
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k, prefac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ccp_buffer, cpp_buffer, &
                                                            max_contraction_i, max_contraction_j, &
                                                            max_contraction_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dijk_contr, dummy_block_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: block_t_i, block_t_j, block_t_k, dijk_j, &
                                                            dijk_k
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, sphi_i, sphi_j, &
                                                            sphi_k, zeti, zetj, zetk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spi, spk, tspj
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbt_type)                                     :: t3c_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t3c_template
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :, :)    :: t3c_der_j
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      END IF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      NULLIFY (qs_kind_set, atomic_kind_set)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t3c_der_i) == [nimg, nimg, 3]))
      CPASSERT(ALL(SHAPE(t3c_der_k) == [nimg, nimg, 3]))

      ALLOCATE (t3c_template(nimg, nimg))
      DO j_img = 1, nimg
         DO i_img = 1, nimg
            CALL dbt_create(t3c_der_i(i_img, j_img, 1), t3c_template(i_img, j_img))
         END DO
      END DO

      CALL alloc_block_3c(t3c_template, nl_3c, basis_i, basis_j, basis_k, qs_env, &
                          potential_parameter, op_pos=op_pos_prv, do_kpoints=do_kpoints, &
                          bounds_i=bounds_i, bounds_j=bounds_j, bounds_k=bounds_k)
      DO i_xyz = 1, 3
         DO j_img = 1, nimg
            DO i_img = 1, nimg
               CALL dbt_copy(t3c_template(i_img, j_img), t3c_der_i(i_img, j_img, i_xyz))
               CALL dbt_copy(t3c_template(i_img, j_img), t3c_der_k(i_img, j_img, i_xyz))
            END DO
         END DO
      END DO

      DO j_img = 1, nimg
         DO i_img = 1, nimg
            CALL dbt_destroy(t3c_template(i_img, j_img))
         END DO
      END DO
      DEALLOCATE (t3c_template)

      IF (nl_3c%sym == symmetric_jk) THEN
         ALLOCATE (t3c_der_j(nimg, nimg, 3))
         DO i_xyz = 1, 3
            DO j_img = 1, nimg
               DO i_img = 1, nimg
                  CALL dbt_create(t3c_der_k(i_img, j_img, i_xyz), t3c_der_j(i_img, j_img, i_xyz))
                  CALL dbt_copy(t3c_der_k(i_img, j_img, i_xyz), t3c_der_j(i_img, j_img, i_xyz))
               END DO
            END DO
         END DO
      END IF

      !Need the max l for each basis for libint and max nset, nco and nsgf for LIBXSMM contraction
      nbasis = SIZE(basis_i)
      max_nsgfi = 0
      max_ncoi = 0
      max_nset = 0
      maxli = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax, &
                                nset=iset, nsgf_set=nsgfi, npgf=npgfi)
         maxli = MAX(maxli, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfi = MAX(max_nsgfi, MAXVAL(nsgfi))
         max_ncoi = MAX(max_ncoi, MAXVAL(npgfi)*ncoset(maxli))
      END DO
      max_nsgfj = 0
      max_ncoj = 0
      maxlj = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax, &
                                nset=jset, nsgf_set=nsgfj, npgf=npgfj)
         maxlj = MAX(maxlj, imax)
         max_nset = MAX(max_nset, jset)
         max_nsgfj = MAX(max_nsgfj, MAXVAL(nsgfj))
         max_ncoj = MAX(max_ncoj, MAXVAL(npgfj)*ncoset(maxlj))
      END DO
      max_nsgfk = 0
      max_ncok = 0
      maxlk = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_k(ibasis)%gto_basis_set, maxl=imax, &
                                nset=kset, nsgf_set=nsgfk, npgf=npgfk)
         maxlk = MAX(maxlk, imax)
         max_nset = MAX(max_nset, kset)
         max_nsgfk = MAX(max_nsgfk, MAXVAL(nsgfk))
         max_ncok = MAX(max_ncok, MAXVAL(npgfk)*ncoset(maxlk))
      END DO
      m_max = maxli + maxlj + maxlk + 1

      !To minimize expensive memory opsand generally optimize contraction, pre-allocate
      !contiguous sphi arrays (and transposed in the cas of sphi_i)

      NULLIFY (tspj, spi, spk)
      ALLOCATE (spi(max_nset, nbasis), tspj(max_nset, nbasis), spk(max_nset, nbasis))

      DO ibasis = 1, nbasis
         DO iset = 1, max_nset
            NULLIFY (spi(iset, ibasis)%array)
            NULLIFY (tspj(iset, ibasis)%array)

            NULLIFY (spk(iset, ibasis)%array)
         END DO
      END DO

      DO ilist = 1, 3
         DO ibasis = 1, nbasis
            IF (ilist == 1) basis_set => basis_i(ibasis)%gto_basis_set
            IF (ilist == 2) basis_set => basis_j(ibasis)%gto_basis_set
            IF (ilist == 3) basis_set => basis_k(ibasis)%gto_basis_set

            DO iset = 1, basis_set%nset

               ncoi = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
               sgfi = basis_set%first_sgf(1, iset)
               egfi = sgfi + basis_set%nsgf_set(iset) - 1

               IF (ilist == 1) THEN
                  ALLOCATE (spi(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spi(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)

               ELSE IF (ilist == 2) THEN
                  ALLOCATE (tspj(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoi))
                  tspj(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoi, sgfi:egfi))

               ELSE
                  ALLOCATE (spk(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spk(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               END IF

            END DO !iset
         END DO !ibasis
      END DO !ilist

      !Init the truncated Coulomb operator
      IF (op_ij == do_potential_truncated .OR. op_jk == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      END IF

      nthread = 1
!$    nthread = omp_get_max_threads()

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,do_kpoints_prv,kp_index_lbounds,kp_index_ubounds,maxli,maxlk,maxlj,bounds_i,&
!$OMP         bounds_j,bounds_k,nimg,basis_i,basis_j,basis_k,dr_ij,dr_jk,dr_ik,ncoset,&
!$OMP         potential_parameter,der_eps,tspj,spi,spk,cell_to_index,max_ncoi,max_nsgfk,&
!$OMP         max_nsgfj,max_ncok,natom,nl_3c,t3c_der_i,t3c_der_k,t3c_der_j) &
!$OMP PRIVATE (lib,nl_3c_iter,ikind,jkind,kkind,iatom,jatom,katom,rij,rjk,rik,cell_j,cell_k,&
!$OMP          prefac,jcell,kcell,first_sgf_i,lmax_i,lmin_i,npgfi,nseti,nsgfi,rpgf_i,set_radius_i,&
!$OMP          sphi_i,zeti,kind_radius_i,first_sgf_j,lmax_j,lmin_j,npgfj,nsetj,nsgfj,rpgf_j,&
!$OMP          set_radius_j,sphi_j,zetj,kind_radius_j,first_sgf_k,lmax_k,lmin_k,npgfk,nsetk,nsgfk,&
!$OMP          rpgf_k,set_radius_k,sphi_k,zetk,kind_radius_k,djk,dij,dik,ncoi,ncoj,ncok,sgfi,sgfj,&
!$OMP          sgfk,dijk_j,dijk_k,ri,rj,rk,max_contraction_i,max_contraction_j,&
!$OMP          max_contraction_k,iset,jset,kset,block_t_i,blk_size,dijk_contr,cpp_buffer,ccp_buffer,&
!$OMP          block_start_j,block_end_j,block_start_k,block_end_k,block_start_i,block_end_i,found,&
!$OMP          dummy_block_t,sp,handle2,mepos,bo,block_t_k,der_ext_i,der_ext_j,der_ext_k,&
!$OMP          block_k_not_zero,der_k_zero,skip,der_j_zero,block_t_j,block_j_not_zero)

      mepos = 0
!$    mepos = omp_get_thread_num()

      CALL cp_libint_init_3eri1(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      !pre-allocate contraction buffers
      ALLOCATE (cpp_buffer(max_nsgfj*max_ncok), ccp_buffer(max_nsgfj*max_nsgfk*max_ncoi))

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)

      !We split the provided bounds among the threads such that each threads works on a different set of atoms
      IF (PRESENT(bounds_i)) THEN
         bo = get_limit(bounds_i(2) - bounds_i(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_i(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bo, bounds_j, bounds_k)
      ELSE IF (PRESENT(bounds_j)) THEN
         bo = get_limit(bounds_j(2) - bounds_j(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_j(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bo, bounds_k)
      ELSE IF (PRESENT(bounds_k)) THEN
         bo = get_limit(bounds_k(2) - bounds_k(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_k(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bo)
      ELSE
         bo = get_limit(natom, nthread, mepos)
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bo, bounds_j, bounds_k)
      END IF

      skip = .FALSE.
      IF (bo(1) > bo(2)) skip = .TRUE.

      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)
         IF (skip) EXIT

         IF (do_kpoints_prv) THEN
            prefac = 0.5_dp
         ELSEIF (nl_3c%sym == symmetric_jk) THEN
            IF (jatom == katom) THEN
               prefac = 0.5_dp
            ELSE
               prefac = 1.0_dp
            END IF
         ELSE
            prefac = 1.0_dp
         END IF

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE

         ELSE
            jcell = 1; kcell = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti, kind_radius=kind_radius_i)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj, kind_radius=kind_radius_j)

         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, lmax=lmax_k, lmin=lmin_k, &
                                npgf=npgfk, nset=nsetk, nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                                sphi=sphi_k, zet=zetk, kind_radius=kind_radius_k)

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         ALLOCATE (max_contraction_i(nseti))
         max_contraction_i = 0.0_dp
         DO iset = 1, nseti
            sgfi = first_sgf_i(1, iset)
            max_contraction_i(iset) = MAXVAL((/(SUM(ABS(sphi_i(:, i))), i=sgfi, sgfi + nsgfi(iset) - 1)/))
         END DO

         ALLOCATE (max_contraction_j(nsetj))
         max_contraction_j = 0.0_dp
         DO jset = 1, nsetj
            sgfj = first_sgf_j(1, jset)
            max_contraction_j(jset) = MAXVAL((/(SUM(ABS(sphi_j(:, i))), i=sgfj, sgfj + nsgfj(jset) - 1)/))
         END DO

         ALLOCATE (max_contraction_k(nsetk))
         max_contraction_k = 0.0_dp
         DO kset = 1, nsetk
            sgfk = first_sgf_k(1, kset)
            max_contraction_k(kset) = MAXVAL((/(SUM(ABS(sphi_k(:, i))), i=sgfk, sgfk + nsgfk(kset) - 1)/))
         END DO

         CALL dbt_blk_sizes(t3c_der_i(jcell, kcell, 1), [iatom, jatom, katom], blk_size)

         ALLOCATE (block_t_i(blk_size(2), blk_size(3), blk_size(1), 3))
         ALLOCATE (block_t_j(blk_size(2), blk_size(3), blk_size(1), 3))
         ALLOCATE (block_t_k(blk_size(2), blk_size(3), blk_size(1), 3))

         block_t_i = 0.0_dp
         block_t_j = 0.0_dp
         block_t_k = 0.0_dp
         block_j_not_zero = .FALSE.
         block_k_not_zero = .FALSE.

         DO iset = 1, nseti

            DO jset = 1, nsetj

               IF (set_radius_j(jset) + set_radius_i(iset) + dr_ij < dij) CYCLE

               DO kset = 1, nsetk

                  IF (set_radius_j(jset) + set_radius_k(kset) + dr_jk < djk) CYCLE
                  IF (set_radius_k(kset) + set_radius_i(iset) + dr_ik < dik) CYCLE

                  ncoi = npgfi(iset)*ncoset(lmax_i(iset))
                  ncoj = npgfj(jset)*ncoset(lmax_j(jset))
                  ncok = npgfk(kset)*ncoset(lmax_k(kset))

                  sgfi = first_sgf_i(1, iset)
                  sgfj = first_sgf_j(1, jset)
                  sgfk = first_sgf_k(1, kset)

                  IF (ncoj*ncok*ncoi > 0) THEN
                     ALLOCATE (dijk_j(ncoj, ncok, ncoi, 3))
                     ALLOCATE (dijk_k(ncoj, ncok, ncoi, 3))
                     dijk_j(:, :, :, :) = 0.0_dp
                     dijk_k(:, :, :, :) = 0.0_dp

                     der_j_zero = .FALSE.
                     der_k_zero = .FALSE.

                     !need positions for libint. Only relative positions are needed => set ri to 0.0
                     ri = 0.0_dp
                     rj = rij ! ri + rij
                     rk = rik ! ri + rik

                     CALL eri_3center_derivs(dijk_j, dijk_k, &
                                             lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), rpgf_j(:, jset), rj, &
                                             lmin_k(kset), lmax_k(kset), npgfk(kset), zetk(:, kset), rpgf_k(:, kset), rk, &
                                             lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), rpgf_i(:, iset), ri, &
                                             djk, dij, dik, lib, potential_parameter, &
                                             der_abc_1_ext=der_ext_j, der_abc_2_ext=der_ext_k)

                     IF (PRESENT(der_eps)) THEN
                        DO i_xyz = 1, 3
                           IF (der_eps > der_ext_j(i_xyz)*(max_contraction_i(iset)* &
                                                           max_contraction_j(jset)* &
                                                           max_contraction_k(kset))) THEN
                              der_j_zero(i_xyz) = .TRUE.
                           END IF
                        END DO

                        DO i_xyz = 1, 3
                           IF (der_eps > der_ext_k(i_xyz)*(max_contraction_i(iset)* &
                                                           max_contraction_j(jset)* &
                                                           max_contraction_k(kset))) THEN
                              der_k_zero(i_xyz) = .TRUE.
                           END IF
                        END DO
                        IF (ALL(der_j_zero) .AND. ALL(der_k_zero)) THEN
                           DEALLOCATE (dijk_j, dijk_k)
                           CYCLE
                        END IF
                     END IF

                     ALLOCATE (dijk_contr(nsgfj(jset), nsgfk(kset), nsgfi(iset)))

                     block_start_j = sgfj
                     block_end_j = sgfj + nsgfj(jset) - 1
                     block_start_k = sgfk
                     block_end_k = sgfk + nsgfk(kset) - 1
                     block_start_i = sgfi
                     block_end_i = sgfi + nsgfi(iset) - 1

                     DO i_xyz = 1, 3
                        IF (der_j_zero(i_xyz)) CYCLE

                        block_j_not_zero(i_xyz) = .TRUE.
                        CALL libxsmm_abc_contract(dijk_contr, dijk_j(:, :, :, i_xyz), tspj(jset, jkind)%array, &
                                                  spk(kset, kkind)%array, spi(iset, ikind)%array, &
                                                  ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                                  nsgfi(iset), cpp_buffer, ccp_buffer)

                        block_t_j(block_start_j:block_end_j, &
                                  block_start_k:block_end_k, &
                                  block_start_i:block_end_i, i_xyz) = &
                           block_t_j(block_start_j:block_end_j, &
                                     block_start_k:block_end_k, &
                                     block_start_i:block_end_i, i_xyz) + &
                           prefac*dijk_contr(:, :, :)

                     END DO

                     DO i_xyz = 1, 3
                        IF (der_k_zero(i_xyz)) CYCLE

                        block_k_not_zero(i_xyz) = .TRUE.
                        CALL libxsmm_abc_contract(dijk_contr, dijk_k(:, :, :, i_xyz), tspj(jset, jkind)%array, &
                                                  spk(kset, kkind)%array, spi(iset, ikind)%array, &
                                                  ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                                  nsgfi(iset), cpp_buffer, ccp_buffer)

                        block_t_k(block_start_j:block_end_j, &
                                  block_start_k:block_end_k, &
                                  block_start_i:block_end_i, i_xyz) = &
                           block_t_k(block_start_j:block_end_j, &
                                     block_start_k:block_end_k, &
                                     block_start_i:block_end_i, i_xyz) + &
                           prefac*dijk_contr(:, :, :)

                     END DO

                     DEALLOCATE (dijk_j, dijk_k, dijk_contr)
                  END IF ! number of triples > 0
               END DO
            END DO
         END DO

         CALL timeset(routineN//"_put_dbcsr", handle2)
!$OMP CRITICAL
         sp = SHAPE(block_t_i(:, :, :, 1))
         sp([2, 3, 1]) = sp

         DO i_xyz = 1, 3
            !Derivatives wrt to center i can be obtained by translational invariance
            IF ((.NOT. block_j_not_zero(i_xyz)) .AND. (.NOT. block_k_not_zero(i_xyz))) CYCLE
            block_t_i(:, :, :, i_xyz) = -(block_t_j(:, :, :, i_xyz) + block_t_k(:, :, :, i_xyz))

            CALL dbt_put_block(t3c_der_i(jcell, kcell, i_xyz), [iatom, jatom, katom], sp, &
                               RESHAPE(block_t_i(:, :, :, i_xyz), SHAPE=sp, ORDER=[2, 3, 1]), &
                               summation=.TRUE.)
         END DO
!$OMP END CRITICAL
!$OMP CRITICAL
         IF (nl_3c%sym == symmetric_jk) THEN

            sp = SHAPE(block_t_j(:, :, :, 1))
            sp([2, 3, 1]) = sp

            DO i_xyz = 1, 3
               IF (.NOT. block_j_not_zero(i_xyz)) CYCLE
               CALL dbt_put_block(t3c_der_j(jcell, kcell, i_xyz), [iatom, jatom, katom], sp, &
                                  RESHAPE(block_t_j(:, :, :, i_xyz), SHAPE=sp, ORDER=[2, 3, 1]), &
                                  summation=.TRUE.)
            END DO
         END IF
!$OMP END CRITICAL
!$OMP CRITICAL
         sp = SHAPE(block_t_k(:, :, :, 1))
         sp([2, 3, 1]) = sp

         DO i_xyz = 1, 3
            IF (.NOT. block_k_not_zero(i_xyz)) CYCLE
            CALL dbt_put_block(t3c_der_k(jcell, kcell, i_xyz), [iatom, jatom, katom], sp, &
                               RESHAPE(block_t_k(:, :, :, i_xyz), SHAPE=sp, ORDER=[2, 3, 1]), &
                               summation=.TRUE.)
         END DO
!$OMP END CRITICAL

         CALL timestop(handle2)

         DEALLOCATE (block_t_i)
         DEALLOCATE (block_t_j)
         DEALLOCATE (block_t_k)

         DEALLOCATE (max_contraction_i, max_contraction_j, max_contraction_k)
      END DO

      CALL cp_libint_cleanup_3eri1(lib)

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)
!$OMP END PARALLEL

      IF (do_kpoints_prv) THEN
         DO i_xyz = 1, 3
            DO kcell = 1, nimg
               DO jcell = 1, nimg
                  ! need half of filter eps because afterwards we add transposed tensor
                  CALL dbt_filter(t3c_der_i(jcell, kcell, i_xyz), filter_eps/2)
                  CALL dbt_filter(t3c_der_k(jcell, kcell, i_xyz), filter_eps/2)
               END DO
            END DO
         END DO

      ELSEIF (nl_3c%sym == symmetric_jk) THEN
         !Add the transpose of t3c_der_j to t3c_der_k to get the fully populated tensor
         CALL dbt_create(t3c_der_k(1, 1, 1), t3c_tmp)
         DO i_xyz = 1, 3
            DO kcell = 1, nimg
               DO jcell = 1, nimg
                  CALL dbt_copy(t3c_der_j(jcell, kcell, i_xyz), t3c_der_k(jcell, kcell, i_xyz), &
                                order=[1, 3, 2], move_data=.TRUE., summation=.TRUE.)
                  CALL dbt_filter(t3c_der_k(jcell, kcell, i_xyz), filter_eps)

                  CALL dbt_copy(t3c_der_i(jcell, kcell, i_xyz), t3c_tmp)
                  CALL dbt_copy(t3c_tmp, t3c_der_i(jcell, kcell, i_xyz), &
                                order=[1, 3, 2], move_data=.TRUE., summation=.TRUE.)
                  CALL dbt_filter(t3c_der_i(jcell, kcell, i_xyz), filter_eps)
               END DO
            END DO
         END DO
         CALL dbt_destroy(t3c_tmp)

      ELSEIF (nl_3c%sym == symmetric_none) THEN
         DO i_xyz = 1, 3
            DO kcell = 1, nimg
               DO jcell = 1, nimg
                  CALL dbt_filter(t3c_der_i(jcell, kcell, i_xyz), filter_eps)
                  CALL dbt_filter(t3c_der_k(jcell, kcell, i_xyz), filter_eps)
               END DO
            END DO
         END DO
      ELSE
         CPABORT("requested symmetric case not implemented")
      END IF

      IF (nl_3c%sym == symmetric_jk) THEN
         DO i_xyz = 1, 3
            DO j_img = 1, nimg
               DO i_img = 1, nimg
                  CALL dbt_destroy(t3c_der_j(i_img, j_img, i_xyz))
               END DO
            END DO
         END DO
      END IF

      DO iset = 1, max_nset
         DO ibasis = 1, nbasis
            IF (ASSOCIATED(spi(iset, ibasis)%array)) DEALLOCATE (spi(iset, ibasis)%array)
            IF (ASSOCIATED(tspj(iset, ibasis)%array)) DEALLOCATE (tspj(iset, ibasis)%array)

            IF (ASSOCIATED(spk(iset, ibasis)%array)) DEALLOCATE (spk(iset, ibasis)%array)
         END DO
      END DO

      DEALLOCATE (spi, tspj, spk)

      CALL timestop(handle)
   END SUBROUTINE build_3c_derivatives

! **************************************************************************************************
!> \brief A test routine to check if GPW integrals can be viable in RI-HFX
!> \param t3c ...
!> \param filter_eps ...
!> \param qs_env ...
!> \param nl_3c ...
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param potential_parameter ...
! **************************************************************************************************
   SUBROUTINE build_3c_integrals_gpw(t3c, filter_eps, qs_env, &
                                     nl_3c, basis_i, basis_j, basis_k, &
                                     potential_parameter)

      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t3c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integrals_gpw'

      INTEGER :: blacs_grid_layout, blk, color, dims(3), group_size, handle, handle2, i, i_RI, &
         iatom, jatom, katom, lb_RI, map_1(1), map_2(2), natom, ngroup, pdims(3), ub_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsize_i, bsize_j, bsize_k, dist_i
      INTEGER, DIMENSION(:), POINTER                     :: dist_mu, dist_nu
      LOGICAL                                            :: afound, blacs_repeatable, pfound
      REAL(dp)                                           :: cutoff_old, fac, relative_cutoff_old
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: e_cutoff_old, wf_vector
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: ablock
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_ext
      TYPE(cp_fm_type), POINTER                          :: cp_fm_dummy
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_ext
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type)                                 :: mat_munu
      TYPE(dbcsr_type)                                   :: mat_munu_desymm
      TYPE(dbt_distribution_type)                        :: dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_type)                                     :: t3c_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_ext
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_ext
      TYPE(pw_p_type)                                    :: pot_g, psi_L, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(task_list_type), POINTER                      :: task_list_ext

      NULLIFY (para_env, para_env_ext, blacs_env_ext, dft_control, atomic_kind_set, qs_kind_set, &
               molecule_kind_set, molecule_set, particle_set, cell, sab_orb_ext, dist_mu, dist_mu, &
               pw_env_ext, poisson_env, auxbas_pw_pool, cp_fm_dummy, pblock)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, molecule_kind_set=molecule_kind_set, molecule_set=molecule_set, &
                      particle_set=particle_set, cell=cell, natom=natom)

      !Best route of action, for now, seems to be the creation of an intermediate tensor that follows
      !the distribution of GPW => the cummunicator is split on subcomms, each subcomm comes with
      !a matrix ditribution => create a tensor with the corresponding pdims and hope that works
      !TODO: if everything works and we go ahead with it, might want to the the memory cut on RI atom

      !For now, we assume the tensor is RI x AO x AO
      CALL dbt_get_info(t3c(1, 1), pdims=pdims, nblks_total=dims)
      ngroup = pdims(1)
      group_size = para_env%num_pe/ngroup
      color = para_env%mepos/group_size
      CALL cp_para_env_split(para_env_ext, para_env, color)

      blacs_grid_layout = BLACS_GRID_SQUARE
      blacs_repeatable = .TRUE.
      CALL cp_blacs_env_create(blacs_env_ext, para_env_ext, &
                               blacs_grid_layout, &
                               blacs_repeatable)

      CALL create_mat_munu(mat_munu, qs_env, qs_env%mp2_env, dft_control, atomic_kind_set, qs_kind_set, &
                           atom2d, molecule_kind_set, molecule_set, sab_orb_ext, particle_set, cell, &
                           blacs_env_ext)
      CALL dbcsr_create(matrix=mat_munu_desymm, template=mat_munu%matrix, matrix_type=dbcsr_type_no_symmetry)

      !Create a tmp tensor to get matching dists
      map_1 = [1]
      map_2 = [2, 3]
      !pdims = [ngroup, group_size, 1]
      CALL dbt_pgrid_create(para_env%group, pdims, pgrid)
      PRINT *, "3D pgrid dimensions: ", pdims

      ALLOCATE (dist_i(dims(1)))
      dist_i(:) = (/(MODULO(i, ngroup), i=0, dims(1) - 1)/) !compatible with sub-comm color
      CALL dbcsr_get_info(mat_munu%matrix, proc_row_dist=dist_mu, proc_col_dist=dist_nu)
      CALL dbt_distribution_new(dist, pgrid, dist_i, dist_mu, dist_nu)

      ALLOCATE (bsize_i(dims(1)))
      ALLOCATE (bsize_j(dims(2)))
      ALLOCATE (bsize_k(dims(3)))
      CALL dbt_get_info(t3c(1, 1), blk_size_1=bsize_i, blk_size_2=bsize_j, blk_size_3=bsize_k)

      CALL dbt_create(t3c_tmp, "TEST", dist, map_1, map_2, bsize_i, bsize_j, bsize_k)
      !Note: that all seems to work in terms of dsitributions, provided that we desymmetrize all of it

      !reserve blocks
      CALL alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, potential_parameter)
      CALL dbt_copy(t3c(1, 1), t3c_tmp, move_data=.TRUE.) !TODO: this seems to crash with certain config of cores

      !TODO: for now use a symmetric sab_orb and task_list because required by integrate_v_rspace
      !      this leads to using dbcsr_desymmetrize many times... need to check the cost of it
      !      and see if we can maybe create a non-symmetric task_list
      CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_ext, pw_env_ext, &
                       auxbas_pw_pool, poisson_env, task_list_ext, rho_r, rho_g, pot_g, psi_L, sab_orb_ext)

      !Loop over the RI basis elements and integrate
      ALLOCATE (wf_vector(SUM(bsize_i(:))))
      DO iatom = 1, natom

         !only act if on correct sub-communicator
         IF (dist_i(iatom) .NE. color) CYCLE

         lb_RI = SUM(bsize_i(1:iatom - 1))
         ub_RI = lb_RI + bsize_i(iatom)
         DO i_RI = lb_RI + 1, ub_RI

            wf_vector = 0.0_dp
            wf_vector(i_RI) = 1.0_dp

            CALL timeset("HFX_actual_integ", handle2)
            CALL mp2_eri_3c_integrate_gpw(cp_fm_dummy, psi_L, rho_g, atomic_kind_set, qs_kind_set, cell, &
                                          dft_control, particle_set, pw_env_ext, wf_vector, poisson_env, rho_r, &
                                          pot_g, potential_parameter, mat_munu, qs_env, task_list_ext)
            CALL timestop(handle2)
            CALL timeset("HFX_desymm", handle2)
            CALL dbcsr_desymmetrize(mat_munu%matrix, mat_munu_desymm)
            CALL dbcsr_filter(mat_munu_desymm, filter_eps)
            CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
            CALL timestop(handle2)

            !DBCSR iterator on mat_munu, only take upper diagonal elements and factor 0.5 on diag, put in t3c_tmp
            CALL timeset("HFX_blk_copy", handle2)
            CALL dbcsr_iterator_start(iter, mat_munu_desymm)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, row=jatom, column=katom, blk=blk)

               fac = 1.0_dp
               IF (jatom == katom) fac = 0.5_dp

               CALL dbcsr_get_block_p(mat_munu_desymm, jatom, katom, pblock, pfound)
               IF (pfound) THEN

                  CALL dbt_get_block(t3c_tmp, [iatom, jatom, katom], ablock, afound)

                  IF (afound) THEN

                     ablock(i_RI - lb_RI, :, :) = fac*pblock(:, :)
                     CALL dbt_put_block(t3c_tmp, [iatom, jatom, katom], &
                                        [SIZE(ablock, 1), SIZE(ablock, 2), SIZE(ablock, 3)], ablock)

                     DEALLOCATE (ablock)
                  END IF
               END IF !pfound

            END DO !iterator
            CALL dbcsr_iterator_stop(iter)
            CALL timestop(handle2)

         END DO !i_RI
      END DO !iatom

      CALL dbt_copy(t3c_tmp, t3c(1, 1), move_data=.TRUE.)

      CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_ext, pw_env_ext, &
                       task_list_ext, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)

      CALL dbt_destroy(t3c_tmp)
      CALL dbt_distribution_destroy(dist)
      CALL dbt_pgrid_destroy(pgrid)

      CALL release_neighbor_list_sets(sab_orb_ext)
      CALL cp_blacs_env_release(blacs_env_ext)
      CALL cp_para_env_release(para_env_ext)

      CALL dbcsr_release(mat_munu%matrix)
      CALL dbcsr_release(mat_munu_desymm)
      DEALLOCATE (mat_munu%matrix)

      CALL timestop(handle)

   END SUBROUTINE build_3c_integrals_gpw

! TODO: for now, and just for checks, this is just a copy from mp2_gpw.F
! **************************************************************************************************
!> \brief Encapsulate the building of dbcsr_matrix mat_munu
!> \param mat_munu ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param dft_control ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param atom2d ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param sab_orb_sub ...
!> \param particle_set ...
!> \param cell ...
!> \param blacs_env_sub ...
!> \param do_ri_aux_basis ...
!> \param do_mixed_basis ...
!> \param group_size_prim ...
!> \param do_im_time ...
!> \author Jan Wilhelm, code by Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE create_mat_munu(mat_munu, qs_env, mp2_env, dft_control, atomic_kind_set, qs_kind_set, &
                              atom2d, molecule_kind_set, &
                              molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub, &
                              do_ri_aux_basis, do_mixed_basis, group_size_prim, do_im_time)

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_ri_aux_basis, do_mixed_basis
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_prim
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_mat_munu'

      INTEGER                                            :: handle, ikind, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: do_kpoints, my_do_im_time, &
                                                            my_do_mixed_basis, my_do_ri_aux_basis
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_sub
      TYPE(distribution_1d_type), POINTER                :: local_molecules_sub, local_particles_sub
      TYPE(distribution_2d_type), POINTER                :: distribution_2d_sub
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      CALL timeset(routineN, handle)

      NULLIFY (basis_set_ri_aux)

      my_do_ri_aux_basis = .FALSE.
      IF (PRESENT(do_ri_aux_basis)) THEN
         my_do_ri_aux_basis = do_ri_aux_basis
      END IF

      my_do_mixed_basis = .FALSE.
      IF (PRESENT(do_mixed_basis)) THEN
         my_do_mixed_basis = do_mixed_basis
      END IF

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) THEN
         my_do_im_time = do_im_time
      END IF

!      do_kpoints = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints .OR. &
!                             SUM(qs_env%mp2_env%ri_rpa_im_time%kp_grid) > 0
      do_kpoints = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints) THEN
         ! please choose EPS_PGF_ORB in QS section smaller than EPS_GRID in WFC_GPW section
         IF (mp2_env%mp2_gpw%eps_grid < dft_control%qs_control%eps_pgf_orb) THEN
            mp2_env%mp2_gpw%eps_grid = dft_control%qs_control%eps_pgf_orb
            CPWARN("WFC_GPW%EPS_GRID has been set to QS%EPS_PGF_ORB")
         END IF
      END IF

      ! hack hack hack XXXXXXXXXXXXXXX ... to be fixed
      dft_control%qs_control%eps_pgf_orb = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_rho_rspace = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_gvg_rspace = mp2_env%mp2_gpw%eps_grid
      CALL init_interaction_radii(dft_control%qs_control, qs_kind_set)

      ! get a distribution_1d
      NULLIFY (local_particles_sub, local_molecules_sub)
      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles_sub, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules_sub, &
                                   force_env_section=qs_env%input)

      ! get a distribution_2d
      NULLIFY (distribution_2d_sub)
      CALL distribute_molecules_2d(cell=cell, &
                                   atomic_kind_set=atomic_kind_set, &
                                   qs_kind_set=qs_kind_set, &
                                   particle_set=particle_set, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   distribution_2d=distribution_2d_sub, &
                                   blacs_env=blacs_env_sub, &
                                   force_env_section=qs_env%input)

      ! Build the sub orbital-orbital overlap neighbor lists
      NULLIFY (sab_orb_sub)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, local_particles_sub, distribution_2d_sub, atomic_kind_set, &
                        molecule_set, molecule_only=.FALSE., particle_set=particle_set)

      ALLOCATE (orb_present(nkind))
      ALLOCATE (orb_radius(nkind))
      ALLOCATE (pair_radius(nkind, nkind))

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         END IF
      END DO

      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

      ! for cubic RPA/GW with kpoints, we need all neighbors and not only the symmetric ones
      !TODO: note that we make it non-symmetric. IT NEEDS TO BE SYMMETRIC FOR INTEGRATION
      CALL build_neighbor_lists(sab_orb_sub, particle_set, atom2d, cell, pair_radius, &
                                mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb_sub", &
                                symmetric=.TRUE.)

      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, orb_radius, pair_radius)

      ! a dbcsr_dist
      ALLOCATE (dbcsr_dist_sub)
      CALL cp_dbcsr_dist2d_to_dist(distribution_2d_sub, dbcsr_dist_sub)

      ! build a dbcsr matrix the hard way
      natom = SIZE(particle_set)
      ALLOCATE (row_blk_sizes(natom))
      IF (my_do_ri_aux_basis) THEN

         ALLOCATE (basis_set_ri_aux(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, basis=basis_set_ri_aux)
         DEALLOCATE (basis_set_ri_aux)

      ELSE IF (my_do_mixed_basis) THEN

         ALLOCATE (basis_set_ri_aux(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, basis=basis_set_ri_aux)
         DEALLOCATE (basis_set_ri_aux)

         ALLOCATE (col_blk_sizes(natom))

         CALL get_particle_set(particle_set, qs_kind_set, nsgf=col_blk_sizes)
         col_blk_sizes = col_blk_sizes*group_size_prim

      ELSE
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes)
      END IF

      NULLIFY (mat_munu%matrix)
      ALLOCATE (mat_munu%matrix)

      IF (my_do_ri_aux_basis) THEN

         CALL dbcsr_create(matrix=mat_munu%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)

      ELSE IF (my_do_mixed_basis) THEN

         CALL dbcsr_create(matrix=mat_munu%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                           nze=0)

      ELSE

         !TODO: note that we make it non-symmetric IT NEEDS TO BE SYMMETRIC FOR INTEGRATION
         CALL dbcsr_create(matrix=mat_munu%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_symmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)

         IF (.NOT. my_do_im_time) THEN

            CALL cp_dbcsr_alloc_block_from_nbl(mat_munu%matrix, sab_orb_sub)

         END IF

      END IF

      DEALLOCATE (row_blk_sizes)

      IF (my_do_mixed_basis) THEN
         DEALLOCATE (col_blk_sizes)
      END IF

      CALL dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (dbcsr_dist_sub)

      CALL distribution_2d_release(distribution_2d_sub)

      CALL distribution_1d_release(local_particles_sub)
      CALL distribution_1d_release(local_molecules_sub)

      CALL timestop(handle)

   END SUBROUTINE create_mat_munu

! **************************************************************************************************
!> \brief Build 3-center integral tensor
!> \param t3c empty DBCSR tensor
!>            Should be of shape (1,1) if no kpoints are used and of shape (nimages, nimages)
!>            if k-points are used
!> \param filter_eps Filter threshold for tensor blocks
!> \param qs_env ...
!> \param nl_3c 3-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param potential_parameter ...
!> \param int_eps neglect integrals smaller than int_eps
!> \param op_pos operator position.
!>        1: calculate (i|jk) integrals,
!>        2: calculate (ij|k) integrals
!> \param do_kpoints ...
!> this routine requires that libint has been static initialised somewhere else
!> \param desymmetrize ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE build_3c_integrals(t3c, filter_eps, qs_env, &
                                 nl_3c, basis_i, basis_j, basis_k, &
                                 potential_parameter, int_eps, &
                                 op_pos, do_kpoints, desymmetrize, &
                                 bounds_i, bounds_j, bounds_k)

      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t3c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: int_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints, desymmetrize
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integrals'

      INTEGER :: block_end_i, block_end_j, block_end_k, block_start_i, block_start_j, &
         block_start_k, egfi, handle, handle2, i, iatom, ibasis, ikind, ilist, imax, iset, jatom, &
         jcell, jkind, jset, katom, kcell, kkind, kset, m_max, max_ncoi, max_ncoj, max_ncok, &
         max_nset, max_nsgfi, max_nsgfj, max_nsgfk, maxli, maxlj, maxlk, mepos, natom, nbasis, &
         ncoi, ncoj, ncok, nimg, nseti, nsetj, nsetk, nthread, op_ij, op_jk, op_pos_prv, sgfi, &
         sgfj, sgfk, unit_id
      INTEGER, DIMENSION(2)                              :: bo
      INTEGER, DIMENSION(3)                              :: blk_size, cell_j, cell_k, &
                                                            kp_index_lbounds, kp_index_ubounds, sp
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: block_not_zero, debug, desymmetrize_prv, &
                                                            do_kpoints_prv, found, skip
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k, prefac, sijk_ext
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ccp_buffer, cpp_buffer, &
                                                            max_contraction_i, max_contraction_j, &
                                                            max_contraction_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_t, dummy_block_t, sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, sphi_i, sphi_j, &
                                                            sphi_k, zeti, zetj, zetk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spi, spk, tspj
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbt_type)                                     :: t_3c_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      debug = .FALSE.

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      IF (PRESENT(desymmetrize)) THEN
         desymmetrize_prv = desymmetrize
      ELSE
         desymmetrize_prv = .TRUE.
      END IF

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      END IF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      NULLIFY (qs_kind_set, atomic_kind_set)

      CALL alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, &
                          potential_parameter, op_pos=op_pos_prv, do_kpoints=do_kpoints, &
                          bounds_i=bounds_i, bounds_j=bounds_j, bounds_k=bounds_k)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t3c) == [nimg, nimg]))

      !Need the max l for each basis for libint and max nset, nco and nsgf for LIBXSMM contraction
      nbasis = SIZE(basis_i)
      max_nsgfi = 0
      max_ncoi = 0
      max_nset = 0
      maxli = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax, &
                                nset=iset, nsgf_set=nsgfi, npgf=npgfi)
         maxli = MAX(maxli, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfi = MAX(max_nsgfi, MAXVAL(nsgfi))
         max_ncoi = MAX(max_ncoi, MAXVAL(npgfi)*ncoset(maxli))
      END DO
      max_nsgfj = 0
      max_ncoj = 0
      maxlj = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax, &
                                nset=jset, nsgf_set=nsgfj, npgf=npgfj)
         maxlj = MAX(maxlj, imax)
         max_nset = MAX(max_nset, jset)
         max_nsgfj = MAX(max_nsgfj, MAXVAL(nsgfj))
         max_ncoj = MAX(max_ncoj, MAXVAL(npgfj)*ncoset(maxlj))
      END DO
      max_nsgfk = 0
      max_ncok = 0
      maxlk = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_k(ibasis)%gto_basis_set, maxl=imax, &
                                nset=kset, nsgf_set=nsgfk, npgf=npgfk)
         maxlk = MAX(maxlk, imax)
         max_nset = MAX(max_nset, kset)
         max_nsgfk = MAX(max_nsgfk, MAXVAL(nsgfk))
         max_ncok = MAX(max_ncok, MAXVAL(npgfk)*ncoset(maxlk))
      END DO
      m_max = maxli + maxlj + maxlk

      !To minimize expensive memory opsand generally optimize contraction, pre-allocate
      !contiguous sphi arrays (and transposed in the cas of sphi_i)

      NULLIFY (tspj, spi, spk)
      ALLOCATE (spi(max_nset, nbasis), tspj(max_nset, nbasis), spk(max_nset, nbasis))

      DO ibasis = 1, nbasis
         DO iset = 1, max_nset
            NULLIFY (spi(iset, ibasis)%array)
            NULLIFY (tspj(iset, ibasis)%array)

            NULLIFY (spk(iset, ibasis)%array)
         END DO
      END DO

      DO ilist = 1, 3
         DO ibasis = 1, nbasis
            IF (ilist == 1) basis_set => basis_i(ibasis)%gto_basis_set
            IF (ilist == 2) basis_set => basis_j(ibasis)%gto_basis_set
            IF (ilist == 3) basis_set => basis_k(ibasis)%gto_basis_set

            DO iset = 1, basis_set%nset

               ncoi = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
               sgfi = basis_set%first_sgf(1, iset)
               egfi = sgfi + basis_set%nsgf_set(iset) - 1

               IF (ilist == 1) THEN
                  ALLOCATE (spi(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spi(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)

               ELSE IF (ilist == 2) THEN
                  ALLOCATE (tspj(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoi))
                  tspj(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoi, sgfi:egfi))

               ELSE
                  ALLOCATE (spk(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spk(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               END IF

            END DO !iset
         END DO !ibasis
      END DO !ilist

      !Init the truncated Coulomb operator
      IF (op_ij == do_potential_truncated .OR. op_jk == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      END IF

      nthread = 1
!$    nthread = omp_get_max_threads()

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,do_kpoints_prv,kp_index_lbounds,kp_index_ubounds,maxli,maxlk,maxlj,bounds_i,&
!$OMP         bounds_j,bounds_k,nimg,basis_i,basis_j,basis_k,dr_ij,dr_jk,dr_ik,ncoset,&
!$OMP         potential_parameter,int_eps,t3c,tspj,spi,spk,debug,cell_to_index,max_ncoi,max_nsgfk,&
!$OMP         max_nsgfj,max_ncok,natom,nl_3c) &
!$OMP PRIVATE (lib,nl_3c_iter,ikind,jkind,kkind,iatom,jatom,katom,rij,rjk,rik,cell_j,cell_k,&
!$OMP          prefac,jcell,kcell,first_sgf_i,lmax_i,lmin_i,npgfi,nseti,nsgfi,rpgf_i,set_radius_i,&
!$OMP          sphi_i,zeti,kind_radius_i,first_sgf_j,lmax_j,lmin_j,npgfj,nsetj,nsgfj,rpgf_j,&
!$OMP          set_radius_j,sphi_j,zetj,kind_radius_j,first_sgf_k,lmax_k,lmin_k,npgfk,nsetk,nsgfk,&
!$OMP          rpgf_k,set_radius_k,sphi_k,zetk,kind_radius_k,djk,dij,dik,ncoi,ncoj,ncok,sgfi,sgfj,&
!$OMP          sgfk,sijk,ri,rj,rk,sijk_ext,block_not_zero,max_contraction_i,max_contraction_j,&
!$OMP          max_contraction_k,iset,jset,kset,block_t,blk_size,sijk_contr,cpp_buffer,ccp_buffer,&
!$OMP          block_start_j,block_end_j,block_start_k,block_end_k,block_start_i,block_end_i,found,&
!$OMP          dummy_block_t,sp,handle2,mepos,bo,skip)

      mepos = 0
!$    mepos = omp_get_thread_num()

      CALL cp_libint_init_3eri(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      !pre-allocate contraction buffers
      ALLOCATE (cpp_buffer(max_nsgfj*max_ncok), ccp_buffer(max_nsgfj*max_nsgfk*max_ncoi))

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)

      !We split the provided bounds among the threads such that each threads works on a different set of atoms
      IF (PRESENT(bounds_i)) THEN
         bo = get_limit(bounds_i(2) - bounds_i(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_i(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bo, bounds_j, bounds_k)
      ELSE IF (PRESENT(bounds_j)) THEN
         bo = get_limit(bounds_j(2) - bounds_j(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_j(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bo, bounds_k)
      ELSE IF (PRESENT(bounds_k)) THEN
         bo = get_limit(bounds_k(2) - bounds_k(1) + 1, nthread, mepos)
         bo(:) = bo(:) + bounds_k(1) - 1
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bo)
      ELSE
         bo = get_limit(natom, nthread, mepos)
         CALL nl_3c_iter_set_bounds(nl_3c_iter, bo, bounds_j, bounds_k)
      END IF

      skip = .FALSE.
      IF (bo(1) > bo(2)) skip = .TRUE.

      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)
         IF (skip) EXIT

         IF (do_kpoints_prv) THEN
            prefac = 0.5_dp
         ELSEIF (nl_3c%sym == symmetric_jk) THEN
            IF (jatom == katom) THEN
               ! factor 0.5 due to double-counting of diagonal blocks
               ! (we desymmetrize by adding transpose)
               prefac = 0.5_dp
            ELSE
               prefac = 1.0_dp
            END IF
         ELSE
            prefac = 1.0_dp
         END IF

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE

         ELSE
            jcell = 1; kcell = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti, kind_radius=kind_radius_i)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj, kind_radius=kind_radius_j)

         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, lmax=lmax_k, lmin=lmin_k, &
                                npgf=npgfk, nset=nsetk, nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                                sphi=sphi_k, zet=zetk, kind_radius=kind_radius_k)

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         ALLOCATE (max_contraction_i(nseti))
         max_contraction_i = 0.0_dp
         DO iset = 1, nseti
            sgfi = first_sgf_i(1, iset)
            max_contraction_i(iset) = MAXVAL((/(SUM(ABS(sphi_i(:, i))), i=sgfi, sgfi + nsgfi(iset) - 1)/))
         END DO

         ALLOCATE (max_contraction_j(nsetj))
         max_contraction_j = 0.0_dp
         DO jset = 1, nsetj
            sgfj = first_sgf_j(1, jset)
            max_contraction_j(jset) = MAXVAL((/(SUM(ABS(sphi_j(:, i))), i=sgfj, sgfj + nsgfj(jset) - 1)/))
         END DO

         ALLOCATE (max_contraction_k(nsetk))
         max_contraction_k = 0.0_dp
         DO kset = 1, nsetk
            sgfk = first_sgf_k(1, kset)
            max_contraction_k(kset) = MAXVAL((/(SUM(ABS(sphi_k(:, i))), i=sgfk, sgfk + nsgfk(kset) - 1)/))
         END DO

         CALL dbt_blk_sizes(t3c(jcell, kcell), [iatom, jatom, katom], blk_size)

         ALLOCATE (block_t(blk_size(2), blk_size(3), blk_size(1)))

         block_t = 0.0_dp
         block_not_zero = .FALSE.

         DO iset = 1, nseti

            DO jset = 1, nsetj

               IF (set_radius_j(jset) + set_radius_i(iset) + dr_ij < dij) CYCLE

               DO kset = 1, nsetk

                  IF (set_radius_j(jset) + set_radius_k(kset) + dr_jk < djk) CYCLE
                  IF (set_radius_k(kset) + set_radius_i(iset) + dr_ik < dik) CYCLE

                  ncoi = npgfi(iset)*ncoset(lmax_i(iset))
                  ncoj = npgfj(jset)*ncoset(lmax_j(jset))
                  ncok = npgfk(kset)*ncoset(lmax_k(kset))

                  sgfi = first_sgf_i(1, iset)
                  sgfj = first_sgf_j(1, jset)
                  sgfk = first_sgf_k(1, kset)

                  IF (ncoj*ncok*ncoi > 0) THEN
                     ALLOCATE (sijk(ncoj, ncok, ncoi))
                     sijk(:, :, :) = 0.0_dp

                     !need positions for libint. Only relative positions are needed => set ri to 0.0
                     ri = 0.0_dp
                     rj = rij ! ri + rij
                     rk = rik ! ri + rik

                     CALL eri_3center(sijk, &
                                      lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), rpgf_j(:, jset), rj, &
                                      lmin_k(kset), lmax_k(kset), npgfk(kset), zetk(:, kset), rpgf_k(:, kset), rk, &
                                      lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), rpgf_i(:, iset), ri, &
                                      djk, dij, dik, lib, potential_parameter, int_abc_ext=sijk_ext)

                     IF (PRESENT(int_eps)) THEN
                        IF (int_eps > sijk_ext*(max_contraction_i(iset)* &
                                                max_contraction_j(jset)* &
                                                max_contraction_k(kset))) THEN
                           DEALLOCATE (sijk)
                           CYCLE
                        END IF
                     END IF

                     block_not_zero = .TRUE.

                     ALLOCATE (sijk_contr(nsgfj(jset), nsgfk(kset), nsgfi(iset)))
                     CALL libxsmm_abc_contract(sijk_contr, sijk, tspj(jset, jkind)%array, &
                                               spk(kset, kkind)%array, spi(iset, ikind)%array, &
                                               ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                               nsgfi(iset), cpp_buffer, ccp_buffer)
                     DEALLOCATE (sijk)

                     block_start_j = sgfj
                     block_end_j = sgfj + nsgfj(jset) - 1
                     block_start_k = sgfk
                     block_end_k = sgfk + nsgfk(kset) - 1
                     block_start_i = sgfi
                     block_end_i = sgfi + nsgfi(iset) - 1

                     block_t(block_start_j:block_end_j, &
                             block_start_k:block_end_k, &
                             block_start_i:block_end_i) = &
                        block_t(block_start_j:block_end_j, &
                                block_start_k:block_end_k, &
                                block_start_i:block_end_i) + &
                        prefac*sijk_contr(:, :, :)

                     DEALLOCATE (sijk_contr)
                  END IF ! number of triples > 0

               END DO

            END DO

         END DO

         IF (block_not_zero) THEN
!$OMP CRITICAL
            CALL timeset(routineN//"_put_dbcsr", handle2)
            IF (debug) THEN
               CALL dbt_get_block(t3c(jcell, kcell), &
                                  [iatom, jatom, katom], dummy_block_t, found=found)
               CPASSERT(found)
            END IF
            sp = SHAPE(block_t)

            sp([2, 3, 1]) = sp

            CALL dbt_put_block(t3c(jcell, kcell), &
                               [iatom, jatom, katom], sp, RESHAPE(block_t, SHAPE=sp, ORDER=[2, 3, 1]), summation=.TRUE.)

            CALL timestop(handle2)
!$OMP END CRITICAL
         END IF

         DEALLOCATE (block_t)
         DEALLOCATE (max_contraction_i, max_contraction_j, max_contraction_k)
      END DO

      CALL cp_libint_cleanup_3eri(lib)

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)
!$OMP END PARALLEL

      IF (nl_3c%sym == symmetric_jk .OR. do_kpoints_prv) THEN
         DO kcell = 1, nimg
            DO jcell = 1, nimg
               ! need half of filter eps because afterwards we add transposed tensor
               CALL dbt_filter(t3c(jcell, kcell), filter_eps/2)
            END DO
         END DO

         IF (desymmetrize_prv) THEN
            ! add transposed of overlap integrals
            CALL dbt_create(t3c(1, 1), t_3c_tmp)
            DO kcell = 1, jcell
               DO jcell = 1, nimg
                  CALL dbt_copy(t3c(jcell, kcell), t_3c_tmp)
                  CALL dbt_copy(t_3c_tmp, t3c(kcell, jcell), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)
                  CALL dbt_filter(t3c(kcell, jcell), filter_eps)
               END DO
            END DO
            DO kcell = jcell + 1, nimg
               DO jcell = 1, nimg
                  CALL dbt_copy(t3c(jcell, kcell), t_3c_tmp)
                  CALL dbt_copy(t_3c_tmp, t3c(kcell, jcell), order=[1, 3, 2], summation=.FALSE., move_data=.TRUE.)
                  CALL dbt_filter(t3c(kcell, jcell), filter_eps)
               END DO
            END DO
            CALL dbt_destroy(t_3c_tmp)
         END IF
      ELSEIF (nl_3c%sym == symmetric_none) THEN
         DO kcell = 1, nimg
            DO jcell = 1, nimg
               CALL dbt_filter(t3c(jcell, kcell), filter_eps)
            END DO
         END DO
      ELSE
         CPABORT("requested symmetric case not implemented")
      END IF

      DO iset = 1, max_nset
         DO ibasis = 1, nbasis
            IF (ASSOCIATED(spi(iset, ibasis)%array)) DEALLOCATE (spi(iset, ibasis)%array)
            IF (ASSOCIATED(tspj(iset, ibasis)%array)) DEALLOCATE (tspj(iset, ibasis)%array)

            IF (ASSOCIATED(spk(iset, ibasis)%array)) DEALLOCATE (spk(iset, ibasis)%array)
         END DO
      END DO

      DEALLOCATE (spi, tspj, spk)

      CALL timestop(handle)
   END SUBROUTINE build_3c_integrals

! **************************************************************************************************
!> \brief Calculates the derivatives of 2-center integrals, wrt to the first center
!> \param t2c_der ...
!> this routine requires that libint has been static initialised somewhere else
!> \param filter_eps Filter threshold for matrix blocks
!> \param qs_env ...
!> \param nl_2c 2-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param potential_parameter ...
!> \param do_kpoints ...
! **************************************************************************************************
   SUBROUTINE build_2c_derivatives(t2c_der, filter_eps, qs_env, &
                                   nl_2c, basis_i, basis_j, &
                                   potential_parameter, do_kpoints)

      TYPE(dbcsr_type), DIMENSION(:, :), INTENT(INOUT)   :: t2c_der
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'build_2c_derivatives'

      INTEGER :: handle, i_xyz, iatom, ibasis, icol, ikind, imax, img, irow, iset, jatom, jkind, &
         jset, m_max, maxli, maxlj, natom, ncoi, ncoj, nimg, nseti, nsetj, op_prv, sgfi, sgfj, &
         unit_id
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmin_i, lmin_j, npgfi, &
                                                            npgfj, nsgfi, nsgfj
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, do_symmetric, found, &
                                                            trans
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dij_contr, dij_rs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dij
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, sphi_i, sphi_j, zeti, &
                                                            zetj
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type), DIMENSION(3)                   :: block_t
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      op_prv = potential_parameter%potential_type

      NULLIFY (qs_kind_set, atomic_kind_set, block_t(1)%block, block_t(2)%block, block_t(3)%block, cell_to_index)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t2c_der) == [nimg, 3]))

      ! check for symmetry
      CPASSERT(SIZE(nl_2c) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=nl_2c, symmetric=do_symmetric)

      IF (do_symmetric) THEN
         DO img = 1, nimg
            !Derivtive matrix is assymetric
            DO i_xyz = 1, 3
               CPASSERT(dbcsr_get_matrix_type(t2c_der(img, i_xyz)) == dbcsr_type_antisymmetric)
            END DO
         END DO
      ELSE
         DO img = 1, nimg
            DO i_xyz = 1, 3
               CPASSERT(dbcsr_get_matrix_type(t2c_der(img, i_xyz)) == dbcsr_type_no_symmetry)
            END DO
         END DO
      END IF

      DO img = 1, nimg
         DO i_xyz = 1, 3
            CALL cp_dbcsr_alloc_block_from_nbl(t2c_der(img, i_xyz), nl_2c)
         END DO
      END DO

      maxli = 0
      DO ibasis = 1, SIZE(basis_i)
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax)
         maxli = MAX(maxli, imax)
      END DO
      maxlj = 0
      DO ibasis = 1, SIZE(basis_j)
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax)
         maxlj = MAX(maxlj, imax)
      END DO

      m_max = maxli + maxlj + 1

      !Init the truncated Coulomb operator
      IF (op_prv == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      CALL cp_libint_init_2eri1(lib, MAX(maxli, maxlj))
      CALL cp_libint_set_contrdepth(lib, 1)

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         IF (do_kpoints_prv) THEN
            img = cell_to_index(cell(1), cell(2), cell(3))
            IF (img > nimg) CYCLE
         ELSE
            img = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj)

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         dab = NORM2(rij)
         trans = do_symmetric .AND. (iatom > jatom)

         DO i_xyz = 1, 3
            CALL dbcsr_get_block_p(matrix=t2c_der(img, i_xyz), &
                                   row=irow, col=icol, BLOCK=block_t(i_xyz)%block, found=found)
            CPASSERT(found)
         END DO

         DO iset = 1, nseti

            ncoi = npgfi(iset)*ncoset(lmax_i(iset))
            sgfi = first_sgf_i(1, iset)

            DO jset = 1, nsetj

               ncoj = npgfj(jset)*ncoset(lmax_j(jset))
               sgfj = first_sgf_j(1, jset)

               IF (ncoi*ncoj > 0) THEN
                  ALLOCATE (dij_contr(nsgfi(iset), nsgfj(jset)))
                  ALLOCATE (dij(ncoi, ncoj, 3))
                  dij(:, :, :) = 0.0_dp

                  ri = 0.0_dp
                  rj = rij

                  CALL eri_2center_derivs(dij, lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), &
                                          rpgf_i(:, iset), ri, lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), &
                                          rpgf_j(:, jset), rj, dab, lib, potential_parameter)

                  DO i_xyz = 1, 3

                     dij_contr(:, :) = 0.0_dp
                     CALL ab_contract(dij_contr, dij(:, :, i_xyz), &
                                      sphi_i(:, sgfi:), sphi_j(:, sgfj:), &
                                      ncoi, ncoj, nsgfi(iset), nsgfj(jset))

                     IF (trans) THEN
                        !if transpose, then -1 factor for antisymmetry
                        ALLOCATE (dij_rs(nsgfj(jset), nsgfi(iset)))
                        dij_rs(:, :) = -1.0_dp*TRANSPOSE(dij_contr)
                     ELSE
                        ALLOCATE (dij_rs(nsgfi(iset), nsgfj(jset)))
                        dij_rs(:, :) = dij_contr
                     END IF

                     CALL block_add("IN", dij_rs, &
                                    nsgfi(iset), nsgfj(jset), block_t(i_xyz)%block, &
                                    sgfi, sgfj, trans=trans)
                     DEALLOCATE (dij_rs)
                  END DO

                  DEALLOCATE (dij, dij_contr)
               END IF
            END DO
         END DO
      END DO

      CALL cp_libint_cleanup_2eri1(lib)

      CALL neighbor_list_iterator_release(nl_iterator)
      DO img = 1, nimg
         DO i_xyz = 1, 3
            CALL dbcsr_finalize(t2c_der(img, i_xyz))
            CALL dbcsr_filter(t2c_der(img, i_xyz), filter_eps)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_2c_derivatives

! **************************************************************************************************
!> \brief ...
!> \param t2c empty DBCSR matrix
!> \param filter_eps Filter threshold for matrix blocks
!> \param qs_env ...
!> \param nl_2c 2-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param potential_parameter ...
!> \param do_kpoints ...
!> this routine requires that libint has been static initialised somewhere else
! **************************************************************************************************
   SUBROUTINE build_2c_integrals(t2c, filter_eps, qs_env, &
                                 nl_2c, basis_i, basis_j, &
                                 potential_parameter, do_kpoints)
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: t2c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'build_2c_integrals'

      INTEGER :: handle, i_diag, iatom, ibasis, icol, ikind, imax, img, irow, iset, jatom, jkind, &
         jset, m_max, maxli, maxlj, natom, ncoi, ncoj, nimg, nseti, nsetj, op_prv, sgfi, sgfj, &
         unit_id
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmin_i, lmin_j, npgfi, &
                                                            npgfj, nsgfi, nsgfj
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, do_symmetric, found, &
                                                            trans
      REAL(KIND=dp)                                      :: dab, regularization_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sij, sij_contr, sij_rs
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, sphi_i, sphi_j, zeti, &
                                                            zetj
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type)                                 :: block_t
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      END IF

      op_prv = potential_parameter%potential_type

      NULLIFY (qs_kind_set, atomic_kind_set, block_t%block, cell_to_index)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
         regularization_RI = qs_env%mp2_env%ri_rpa_im_time%regularization_RI
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t2c) == [nimg]))

      ! check for symmetry
      CPASSERT(SIZE(nl_2c) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=nl_2c, symmetric=do_symmetric)

      IF (do_symmetric) THEN
         DO img = 1, nimg
            CPASSERT(dbcsr_has_symmetry(t2c(img)))
         END DO
      ELSE
         DO img = 1, nimg
            CPASSERT(.NOT. dbcsr_has_symmetry(t2c(img)))
         END DO
      END IF

      DO img = 1, nimg
         CALL cp_dbcsr_alloc_block_from_nbl(t2c(img), nl_2c)
      END DO

      maxli = 0
      DO ibasis = 1, SIZE(basis_i)
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax)
         maxli = MAX(maxli, imax)
      END DO
      maxlj = 0
      DO ibasis = 1, SIZE(basis_j)
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax)
         maxlj = MAX(maxlj, imax)
      END DO

      m_max = maxli + maxlj

      !Init the truncated Coulomb operator
      IF (op_prv == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      CALL cp_libint_init_2eri(lib, MAX(maxli, maxlj))
      CALL cp_libint_set_contrdepth(lib, 1)

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         IF (do_kpoints_prv) THEN
            img = cell_to_index(cell(1), cell(2), cell(3))
            IF (img > nimg) CYCLE
         ELSE
            img = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj)

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         dab = NORM2(rij)

         CALL dbcsr_get_block_p(matrix=t2c(img), &
                                row=irow, col=icol, BLOCK=block_t%block, found=found)
         CPASSERT(found)
         trans = do_symmetric .AND. (iatom > jatom)

         DO iset = 1, nseti

            ncoi = npgfi(iset)*ncoset(lmax_i(iset))
            sgfi = first_sgf_i(1, iset)

            DO jset = 1, nsetj

               ncoj = npgfj(jset)*ncoset(lmax_j(jset))
               sgfj = first_sgf_j(1, jset)

               IF (ncoi*ncoj > 0) THEN
                  ALLOCATE (sij_contr(nsgfi(iset), nsgfj(jset)))
                  sij_contr(:, :) = 0.0_dp

                  ALLOCATE (sij(ncoi, ncoj))
                  sij(:, :) = 0.0_dp

                  ri = 0.0_dp
                  rj = rij

                  CALL eri_2center(sij, lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), &
                                   rpgf_i(:, iset), ri, lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), &
                                   rpgf_j(:, jset), rj, dab, lib, potential_parameter)

                  CALL ab_contract(sij_contr, sij, &
                                   sphi_i(:, sgfi:), sphi_j(:, sgfj:), &
                                   ncoi, ncoj, nsgfi(iset), nsgfj(jset))

                  DEALLOCATE (sij)
                  IF (trans) THEN
                     ALLOCATE (sij_rs(nsgfj(jset), nsgfi(iset)))
                     sij_rs(:, :) = TRANSPOSE(sij_contr)
                  ELSE
                     ALLOCATE (sij_rs(nsgfi(iset), nsgfj(jset)))
                     sij_rs(:, :) = sij_contr
                  END IF

                  DEALLOCATE (sij_contr)

                  ! RI regularization
                  IF (do_kpoints_prv .AND. cell(1) == 0 .AND. cell(2) == 0 .AND. cell(3) == 0 .AND. &
                      iatom == jatom .AND. iset == jset) THEN
                     DO i_diag = 1, nsgfi(iset)
                        sij_rs(i_diag, i_diag) = sij_rs(i_diag, i_diag) + &
                                                 regularization_RI*MAX(1.0_dp, 1.0_dp/MINVAL(zeti(:, iset)))
                     END DO
                  END IF

                  CALL block_add("IN", sij_rs, &
                                 nsgfi(iset), nsgfj(jset), block_t%block, &
                                 sgfi, sgfj, trans=trans)
                  DEALLOCATE (sij_rs)

               END IF
            END DO
         END DO
      END DO

      CALL cp_libint_cleanup_2eri(lib)

      CALL neighbor_list_iterator_release(nl_iterator)
      DO img = 1, nimg
         CALL dbcsr_finalize(t2c(img))
         CALL dbcsr_filter(t2c(img), filter_eps)
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_2c_integrals

! **************************************************************************************************
!> \brief ...
!> \param tensor tensor with data. Data is cleared after compression.
!> \param blk_indices ...
!> \param compressed compressed tensor data
!> \param eps all entries < eps are discarded
!> \param memory ...
! **************************************************************************************************
   SUBROUTINE compress_tensor(tensor, blk_indices, compressed, eps, memory)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: blk_indices
      TYPE(hfx_compression_type), INTENT(INOUT)          :: compressed
      REAL(dp), INTENT(IN)                               :: eps
      REAL(dp), INTENT(INOUT)                            :: memory

      INTEGER                                            :: buffer_left, buffer_size, buffer_start, &
                                                            i, iblk, memory_usage, nbits, nblk, &
                                                            nints, offset, shared_offset
      INTEGER(int_8)                                     :: estimate_to_store_int, &
                                                            storage_counter_integrals
      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: found
      REAL(dp)                                           :: spherical_estimate
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: blk_data
      REAL(dp), DIMENSION(:), POINTER                    :: blk_data_1d
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(hfx_cache_type), DIMENSION(:), POINTER        :: integral_caches
      TYPE(hfx_cache_type), POINTER                      :: maxval_cache
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container

      CALL dealloc_containers(compressed, memory_usage)
      CALL alloc_containers(compressed, 1)

      maxval_container => compressed%maxval_container(1)
      integral_containers => compressed%integral_containers(:, 1)

      CALL hfx_init_container(maxval_container, memory_usage, .FALSE.)
      DO i = 1, 64
         CALL hfx_init_container(integral_containers(i), memory_usage, .FALSE.)
      END DO

      maxval_cache => compressed%maxval_cache(1)
      integral_caches => compressed%integral_caches(:, 1)

      IF (ALLOCATED(blk_indices)) DEALLOCATE (blk_indices)
      ALLOCATE (blk_indices(dbt_get_num_blocks(tensor), 3))
      shared_offset = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor,blk_indices,shared_offset) &
!$OMP PRIVATE(iter,ind,offset,nblk,iblk)
      CALL dbt_iterator_start(iter, tensor)
      nblk = dbt_iterator_num_blocks(iter)
!$OMP CRITICAL
      offset = shared_offset
      shared_offset = shared_offset + nblk
!$OMP END CRITICAL
      DO iblk = 1, nblk
         CALL dbt_iterator_next_block(iter, ind)
         blk_indices(offset + iblk, :) = ind(:)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      ! Can not use the tensor iterator here because the order of the blocks is not guaranteed.
      DO i = 1, SIZE(blk_indices, 1)
         ind = blk_indices(i, :)
         CALL dbt_get_block(tensor, ind, blk_data, found)
         CPASSERT(found)
         nints = SIZE(blk_data)
         blk_data_1d(1:nints) => blk_data
         spherical_estimate = MAXVAL(ABS(blk_data_1d))
         IF (spherical_estimate == 0.0_dp) spherical_estimate = TINY(spherical_estimate)
         estimate_to_store_int = EXPONENT(spherical_estimate)
         estimate_to_store_int = MAX(estimate_to_store_int, -15_int_8)

         CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                           maxval_cache, maxval_container, memory_usage, &
                                           .FALSE.)

         spherical_estimate = SET_EXPONENT(1.0_dp, estimate_to_store_int + 1)

         nbits = EXPONENT(ANINT(spherical_estimate/eps)) + 1
         IF (nbits > 64) THEN
            CALL cp_abort(__LOCATION__, &
                          "Overflow during tensor compression. Please use a larger EPS_FILTER or EPS_STORAGE_SCALING")
         END IF

         buffer_left = nints
         buffer_start = 1

         DO WHILE (buffer_left > 0)
            buffer_size = MIN(buffer_left, cache_size)
            CALL hfx_add_mult_cache_elements(blk_data_1d(buffer_start:), &
                                             buffer_size, nbits, &
                                             integral_caches(nbits), &
                                             integral_containers(nbits), &
                                             eps, 1.0_dp, &
                                             memory_usage, &
                                             .FALSE.)
            buffer_left = buffer_left - buffer_size
            buffer_start = buffer_start + buffer_size
         END DO

         NULLIFY (blk_data_1d); DEALLOCATE (blk_data)
      END DO

      CALL dbt_clear(tensor)

      storage_counter_integrals = memory_usage*cache_size
      memory = memory + REAL(storage_counter_integrals, dp)/1024/128
      !WRITE (UNIT=iw, FMT="((T3,A,T60,I21))") &
      !   "HFX_MEM_INFO| Total memory consumption ERI's RAM [MiB]:            ", memory

      CALL hfx_flush_last_cache(6, maxval_cache, maxval_container, memory_usage, &
                                .FALSE.)
      DO i = 1, 64
         CALL hfx_flush_last_cache(i, integral_caches(i), integral_containers(i), &
                                   memory_usage, .FALSE.)
      END DO

      CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, memory_usage, .FALSE.)
      DO i = 1, 64
         CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                            memory_usage, .FALSE.)
      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor empty tensor which is filled by decompressed data
!> \param blk_indices indices of blocks to be reserved
!> \param compressed compressed data
!> \param eps all entries < eps are discarded
! **************************************************************************************************
   SUBROUTINE decompress_tensor(tensor, blk_indices, compressed, eps)

      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, DIMENSION(:, :)                           :: blk_indices
      TYPE(hfx_compression_type), INTENT(INOUT)          :: compressed
      REAL(dp), INTENT(IN)                               :: eps

      INTEGER                                            :: A, b, buffer_left, buffer_size, &
                                                            buffer_start, i, memory_usage, nbits, &
                                                            nblk_per_thread, nints
      INTEGER(int_8)                                     :: estimate_to_store_int
      INTEGER, DIMENSION(3)                              :: blk_size, ind
      REAL(dp)                                           :: spherical_estimate
      REAL(dp), ALLOCATABLE, DIMENSION(:), TARGET        :: blk_data
      REAL(dp), DIMENSION(:, :, :), POINTER              :: blk_data_3d
      TYPE(hfx_cache_type), DIMENSION(:), POINTER        :: integral_caches
      TYPE(hfx_cache_type), POINTER                      :: maxval_cache
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container

      maxval_cache => compressed%maxval_cache(1)
      maxval_container => compressed%maxval_container(1)
      integral_caches => compressed%integral_caches(:, 1)
      integral_containers => compressed%integral_containers(:, 1)

      memory_usage = 0

      CALL hfx_decompress_first_cache(6, maxval_cache, maxval_container, memory_usage, .FALSE.)

      DO i = 1, 64
         CALL hfx_decompress_first_cache(i, integral_caches(i), integral_containers(i), &
                                         memory_usage, .FALSE.)
      END DO

!TODO: Parallelize creation of block list.
!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor,blk_indices) PRIVATE(nblk_per_thread,A,b)
      nblk_per_thread = SIZE(blk_indices, 1)/omp_get_num_threads() + 1
      a = omp_get_thread_num()*nblk_per_thread + 1
      b = MIN(a + nblk_per_thread, SIZE(blk_indices, 1))
      CALL dbt_reserve_blocks(tensor, blk_indices(a:b, :))
!$OMP END PARALLEL

      ! Can not use the tensor iterator here because the order of the blocks is not guaranteed.
      DO i = 1, SIZE(blk_indices, 1)
         ind = blk_indices(i, :)
         CALL dbt_blk_sizes(tensor, ind, blk_size)
         nints = PRODUCT(blk_size)
         CALL hfx_get_single_cache_element( &
            estimate_to_store_int, 6, &
            maxval_cache, maxval_container, memory_usage, &
            .FALSE.)

         spherical_estimate = SET_EXPONENT(1.0_dp, estimate_to_store_int + 1)

         nbits = EXPONENT(ANINT(spherical_estimate/eps)) + 1

         buffer_left = nints
         buffer_start = 1

         ALLOCATE (blk_data(nints))
         DO WHILE (buffer_left > 0)
            buffer_size = MIN(buffer_left, cache_size)
            CALL hfx_get_mult_cache_elements(blk_data(buffer_start), &
                                             buffer_size, nbits, &
                                             integral_caches(nbits), &
                                             integral_containers(nbits), &
                                             eps, 1.0_dp, &
                                             memory_usage, &
                                             .FALSE.)
            buffer_left = buffer_left - buffer_size
            buffer_start = buffer_start + buffer_size
         END DO

         blk_data_3d(1:blk_size(1), 1:blk_size(2), 1:blk_size(3)) => blk_data
         CALL dbt_put_block(tensor, ind, blk_size, blk_data_3d)
         NULLIFY (blk_data_3d); DEALLOCATE (blk_data)
      END DO

      CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, memory_usage, .FALSE.)
      DO i = 1, 64
         CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                            memory_usage, .FALSE.)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param nze ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE get_tensor_occupancy(tensor, nze, occ)
      TYPE(dbt_type), INTENT(IN)                         :: tensor
      INTEGER(int_8), INTENT(OUT)                        :: nze
      REAL(dp), INTENT(OUT)                              :: occ

      INTEGER, DIMENSION(dbt_ndims(tensor))              :: dims

      nze = dbt_get_nze_total(tensor)
      CALL dbt_get_info(tensor, nfull_total=dims)
      occ = REAL(nze, dp)/PRODUCT(REAL(dims, dp))

   END SUBROUTINE

END MODULE
