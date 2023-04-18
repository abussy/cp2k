!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX and K-points.
!> \auhtor Augustin Bussy (01.2023)
! **************************************************************************************************

MODULE hfx_ri_kp
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_logical_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_clear, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_dot, dbcsr_filter, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_scale, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_clear, dbt_contract, &
        dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, dbt_destroy, &
        dbt_distribution_destroy, dbt_distribution_new, dbt_distribution_type, dbt_filter, &
        dbt_finalize, dbt_get_block, dbt_get_info, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, &
        dbt_pgrid_type, dbt_put_block, dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE hfx_ri,                          ONLY: get_idx_to_atom,&
                                              hfx_ri_pre_scf_calc_tensors
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: hfx_ri_do_2c_cholesky,&
                                              hfx_ri_do_2c_diag,&
                                              hfx_ri_do_2c_iter
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE iterate_matrix,                  ONLY: invert_hotelling
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_density_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tensors,                      ONLY: build_2c_neighbor_lists,&
                                              build_3c_neighbor_lists,&
                                              get_3c_iterator_info,&
                                              get_tensor_occupancy,&
                                              neighbor_list_3c_destroy,&
                                              neighbor_list_3c_iterate,&
                                              neighbor_list_3c_iterator_create,&
                                              neighbor_list_3c_iterator_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              distribution_2d_create,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_iterator_type,&
                                              neighbor_list_3c_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_ri_update_ks_kp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_kp'
CONTAINS

! NOTES: for a start, we do not seek performance, but accuracy. So in this first implementation,
!        we give little consideration to batching, load balance and such.
!        We also put everything here, even if there is some code replication with the original RI_HFX
!        We will only work in the RHO flavor
!        For now, we will also always assume that there is a single para_env, and that there is no
!        K-point subgroup. This might change in the future

! **************************************************************************************************
!> \brief I_1nitialize the ri_data for K-point. For now, we take the normal, usual existing ri_data
!>        and we adapt it to our needs
!> \param t_3c_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE adapt_ri_data_to_kp(t_3c_template, ri_data, qs_env)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i_img, i_RI, iatom, natom, nblks_RI, &
                                                            nimg, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, dist1, dist2, dist3, &
                                                            dist_RI_ext
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pdims
      TYPE(dbt_distribution_type)                        :: dist
      TYPE(dbt_pgrid_type)                               :: pgrid, pgrid_2d
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (dft_control, para_env)

      !The main thing that we need to do is to allocate more space for the integrals, such that there
      !is room for each periodic image. Note that we only go in 1D, i.e. we store (mu^0 sigma^a|P^0),
      !and (P^0|Q^a) => the RI basis is always in the main cell.

      !Get kpoint info
      CALL get_qs_env(qs_env, dft_control=dft_control, natom=natom, para_env=para_env, nkind=nkind)
      nimg = ri_data%nimg

      CALL dbt_destroy(ri_data%t_3c_int_ctr_1(1, 1))
      DEALLOCATE (ri_data%t_3c_int_ctr_1)
      ALLOCATE (ri_data%t_3c_int_ctr_1(1, nimg))

      !Along the RI direction we have basis elements spread accross ncell_RI images. Take that into
      !account, while making sure that an atom tuple iatom, jatom, katom is always on the same proc
      CALL create_3c_tensor(t_3c_template, dist1, dist2, dist3, &
                            ri_data%pgrid_1, ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, &
                            ri_data%bsizes_AO_split, [1, 2], [3], name="(AO RI | AO)")

      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(nblks_RI*ri_data%ncell_RI))
      ALLOCATE (dist_RI_ext(nblks_RI*ri_data%ncell_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI_split(:)
         dist_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = dist2(:)
      END DO

      CALL dbt_distribution_new(dist, ri_data%pgrid_1, dist1, dist_RI_ext, dist3)
      CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), "KP 3c ints", dist, [1, 2], [3], &
                      ri_data%bsizes_AO_split, bsizes_RI_ext, ri_data%bsizes_AO_split)

      DO i_img = 2, nimg
         CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_1(1, i_img))
      END DO
      CALL dbt_distribution_destroy(dist)
      DEALLOCATE (dist1, dist2, dist3, dist_RI_ext)

      !Because 3c tensor are atom specific, there is only one block along one of the three direction
      !we therefore create new pgrids with this in mind
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env, pdims_2d, pgrid_2d, &
                            tensor_dims=[SIZE(bsizes_RI_ext), SIZE(ri_data%bsizes_AO_split)])
      CALL dbt_pgrid_destroy(pgrid_2d)

      pdims(1) = 1 !dimension along which there is mu^0, defining the atom specific RI basis
      pdims(2) = pdims_2d(1)
      pdims(3) = pdims_2d(2)
      CALL dbt_pgrid_create(para_env, pdims, pgrid)

      CALL dbt_destroy(ri_data%t_3c_int_ctr_2(1, 1))
      DEALLOCATE (ri_data%t_3c_int_ctr_2)
      ALLOCATE (ri_data%t_3c_int_ctr_2(1, 1))
      CALL create_3c_tensor(ri_data%t_3c_int_ctr_2(1, 1), dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO_split, bsizes_RI_ext, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbt_pgrid_destroy(pgrid)

      CALL dbt_destroy(ri_data%t_2c_inv(1, 1))
      DEALLOCATE (ri_data%t_2c_inv)
      ALLOCATE (ri_data%t_2c_inv(1, natom))

      !We use full block sizes for the 2c quantities
      DEALLOCATE (bsizes_RI_ext)
      nblks_RI = SIZE(ri_data%bsizes_RI)
      ALLOCATE (bsizes_RI_ext(nblks_RI*ri_data%ncell_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI(:)
      END DO

      CALL create_2c_tensor(ri_data%t_2c_inv(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)
      DO iatom = 2, natom
         CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_inv(1, iatom))
      END DO

      ALLOCATE (ri_data%kp_cost(natom, natom, nimg))
      ri_data%kp_cost = 0.0_dp

   END SUBROUTINE adapt_ri_data_to_kp

! **************************************************************************************************
!> \brief The pre-scf steps for RI-HFX k-points calculation. Namely the calculation of the integrals
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_kp(ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ri_pre_scf_kp'

      INTEGER                                            :: handle, i_img, i_RI, iatom, ind(3), &
                                                            j_img, natom, nblks_RI, ncell_RI, &
                                                            nimg, nkind
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: t_2c_op_pot, t_2c_op_RI
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(dbt_type)                                     :: t_3c_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_int
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, natom=natom, nkind=nkind)

      !We do all the checks on what we allow in this initial implementation
      IF (ri_data%flavor .NE. ri_pmat) CPABORT("K-points RI-HFX only with RHO flavor")
      IF (ri_data%same_op) ri_data%same_op = .FALSE. !force the full calculation with RI metric
      IF (ABS(ri_data%eps_pgf_orb - dft_control%qs_control%eps_pgf_orb) > 1.0E-16_dp) &
         CPABORT("RI%EPS_PGF_ORB and QS%EPS_PGF_ORB must be identical for RI-HFX k-points")

      CALL get_kp_and_ri_images(ri_data, qs_env)
      nimg = ri_data%nimg
      ncell_RI = ri_data%ncell_RI

      !Calculate the integrals
      ALLOCATE (t_2c_op_pot(nimg), t_2c_op_RI(nimg))
      ALLOCATE (t_3c_int(ncell_RI, nimg))
      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int, do_kpoints=.TRUE.)

      !Make sure the internals have the k-point format
      CALL adapt_ri_data_to_kp(t_3c_tmp, ri_data, qs_env)

      !We copy the 3c integrals, while mapping the image of P^i to the correct RI block
      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      DO j_img = 1, nimg
         DO i_RI = 1, ncell_RI
            CALL dbt_copy(t_3c_int(i_RI, j_img), t_3c_tmp, order=[2, 1, 3], move_data=.TRUE.) !now (mu^0 | P^i sigma^j)
            CALL dbt_destroy(t_3c_int(i_RI, j_img))

            !loop over t_3c_tmp and copy blocks over
!$OMP PARALLEL DEFAULT(NONE) SHARED(t_3c_tmp,ri_data,nblks_RI,j_img,i_RI) PRIVATE(iter,ind,blk,found)
            CALL dbt_iterator_start(iter, t_3c_tmp)
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, ind)
               CALL dbt_get_block(t_3c_tmp, ind, blk, found)
               IF (.NOT. found) CYCLE

               !distributions guaranteed to match
               CALL dbt_put_block(ri_data%t_3c_int_ctr_1(1, j_img), &
                                  [ind(1), (i_RI - 1)*nblks_RI + ind(2), ind(3)], SHAPE(blk), blk)
               DEALLOCATE (blk)
            END DO
            CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
         END DO
      END DO
      CALL dbt_destroy(t_3c_tmp)

      !For each atom i, we calculate the inverse RI metric (P^0 | Q^0)^-1 according to the MIC
      DO iatom = 1, natom
         CALL get_ext_2c_int(ri_data%t_2c_inv(1, iatom), t_2c_op_RI, iatom, iatom, 1, ri_data, qs_env, .TRUE.)
      END DO

      DO i_img = 1, nimg
         CALL dbcsr_release(t_2c_op_RI(i_img))
      END DO

      !Store the HFX pot for each periodic image
      ALLOCATE (ri_data%mat_2c_pot(1, nimg))
      DO i_img = 1, nimg
         CALL dbcsr_create(ri_data%mat_2c_pot(1, i_img), template=t_2c_op_pot(i_img))
         CALL dbcsr_copy(ri_data%mat_2c_pot(1, i_img), t_2c_op_pot(i_img))
         CALL dbcsr_release(t_2c_op_pot(i_img))
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_pre_scf_kp

! **************************************************************************************************
!> \brief Update the KS matrices for each real-space image
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_kp(qs_env, ri_data, ks_matrix, ehfx, &
                                  geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_kp'

      INTEGER                                            :: apc_img, b_img, group_size, handle, &
                                                            handle2, i_img, i_spin, iatom, iblk, &
                                                            igroup, jatom, mb_img, natom, ngroups, &
                                                            nimg
      INTEGER(int_8)                                     :: nflop, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_ranges, dist1, dist2, idx_to_at_AO
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ac_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: sparsity_pattern
      REAL(dp)                                           :: etmp, fac, occ, pref, t1, t2, t3, t4
      TYPE(cp_1d_logical_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: iatom_to_subgroup
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(dbcsr_type)                                   :: ks_desymm, rho_desymm
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_2c_pot
      TYPE(dbt_type)                                     :: t_3c_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: ints_stack, ks_t_split, rho_stack, &
                                                            t_2c_ao_tmp, t_2c_inv, t_2c_work, &
                                                            t_3c_int, t_3c_work_2, t_3c_work_3
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: ks_t, ks_t_sub, rho_ao_t, t_3c_apc, &
                                                            t_3c_apc_sub
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(section_vals_type), POINTER                   :: hfx_section

      NULLIFY (kpoints, dft_control, para_env, para_env_sub, blacs_env_sub, hfx_section)

      CALL timeset(routineN, handle)

      !TODO: deal with anti-symmetric matrices as in standard RI-HFX (for TDDFT) ?

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env, natom=natom, kpoints=kpoints)

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_kp(ri_data, qs_env)
      END IF
      nimg = ri_data%nimg

      !We need to calculate the KS matrix for each periodic cell with index b: F_mu^0,nu^b
      !F_mu^0,nu^b = -0.5 sum_a,c P_sigma^0,lambda^c (mu^0, sigma^a| P^0) V_P^0,Q^b (Q^b| nu^b lambda^a+c)

      !We use a local RI basis set for each atom in the system, which inlcudes RI basis elements for
      !each neighboring atom standing within the KIND radius (decay of Gaussian with smallest exponent)

      !We also limit the number of periodic images we consider accorrding to the HFX potentail in the
      !RI basis, because if V_P^0,Q^b is zero everywhere, then image b can be ignored (RI basis less diffuse)

      !We manage to calculate each KS matrix doing a double loop on iamges, and a double loop on atoms
      !First, we pre-contract and store P_sigma^0,lambda^c (mu^0, sigma^a| P^0) into T_mu^0,lambda^a+c,P^0
      !Then, we loop over b_img, iatom, jatom to calculate V_P^0,Q^b
      !Finally, we do an additional loop over a+c images where we do V_P^0,Q^b (Q^b| nu^b lambda^a+c)
      !and the final contraction with T_mu^0,lambda^a+c,P^0

      ALLOCATE (rho_ao_t(nspins, nimg), ks_t(nspins, nimg))
      CALL create_2c_tensor(rho_ao_t(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      IF (nspins == 2) CALL dbt_create(rho_ao_t(1, 1), rho_ao_t(2, 1))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ks_matrix(i_spin, i_img)%matrix, ks_t(i_spin, i_img))
            IF (i_img > 1) CALL dbt_create(rho_ao_t(1, 1), rho_ao_t(i_spin, i_img))
         END DO
      END DO
      CALL get_pmat_images(rho_ao_t, ri_data, qs_env)

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      !First we calculate and store T^1_mu^0,lambda^a+c,P = P_mu^0,lambda^c * (mu_0 sigma^a | P)
      !To avoid doing nimg**2 tiny contractions that do not scale well with a large number of CPUs,
      !we instead do a single loop over the a+c image index. For each a+c, we get a list of allowed
      !combination of a,c indices. Then we build TAS tensors P_mu^0,lambda^c with all concerned c's
      !and (mu^0 sigma^a | P) with all a's. Then we perform a single contraction with larger tensors,
      !were the sum over a,c is automatically taken care of

      !Create the stack tensors in the approriate distribution
      ALLOCATE (rho_stack(2), ints_stack(2))
      CALL get_stack_tensors(rho_stack, ints_stack, rho_ao_t(1, 1), ri_data%t_3c_int_ctr_1(1, 1), ri_data, qs_env)

      ALLOCATE (t_3c_apc(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_apc(i_spin, i_img))
         END DO
      END DO
      CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), t_3c_tmp)
      t1 = m_walltime()
      CALL timeset(routineN//"_3c_1", handle2)

      !Loop over a+c images
      !TODO: n_mem batching? probably unnecessary
      ALLOCATE (ac_pairs(2, nimg))
      DO apc_img = 1, nimg
         IF (ri_data%present_images(apc_img) == 0) CYCLE
         CALL get_ac_pairs(ac_pairs, apc_img, qs_env)
         CALL fill_ints_stack(ints_stack(1), ri_data%t_3c_int_ctr_1(1, :), ac_pairs, ri_data)
         CALL dbt_copy(ints_stack(1), ints_stack(2), move_data=.TRUE.)

         DO i_spin = 1, nspins
            CALL fill_rho_stack(rho_stack(1), rho_ao_t(i_spin, :), ac_pairs, ri_data)
            CALL dbt_copy(rho_stack(1), rho_stack(2), move_data=.TRUE.)

            CALL dbt_batched_contract_init(t_3c_tmp)
            CALL dbt_contract(1.0_dp, ints_stack(2), rho_stack(2), &
                              0.0_dp, t_3c_tmp, map_1=[1, 2], map_2=[3], &
                              contract_1=[3], notcontract_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], &
                              filter_eps=ri_data%filter_eps, flop=nflop, &
                              move_data=.TRUE.)
            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
            CALL dbt_batched_contract_finalize(t_3c_tmp)

            CALL dbt_copy(t_3c_tmp, t_3c_apc(i_spin, apc_img), move_data=.TRUE.)
         END DO
      END DO

      CALL dbt_destroy(rho_stack(1))
      CALL dbt_destroy(rho_stack(2))
      CALL dbt_destroy(ints_stack(1))
      CALL dbt_destroy(ints_stack(2))
      CALL dbt_destroy(t_3c_tmp)
      CALL timestop(handle2)

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF%RI")
      CALL section_vals_val_get(hfx_section, "KP_NGROUPS", i_val=ngroups)

      IF (MOD(para_env%num_pe, ngroups) .NE. 0) THEN
         CPWARN("KP_NGROUPS must be an integer divisor of the total number of MPI ranks. It was set to 1.")
         ngroups = 1
         CALL section_vals_val_set(hfx_section, "KP_NGROUPS", i_val=ngroups)
      END IF
      group_size = para_env%num_pe/ngroups
      igroup = para_env%mepos/group_size

      ALLOCATE (para_env_sub)
      CALL para_env_sub%from_split(para_env, igroup)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! The sparsity pattern of each iatom, jatom pair, on each b_img, and on which subgroup
      ALLOCATE (sparsity_pattern(natom, natom, nimg))
      CALL get_sparsity_pattern(sparsity_pattern, ri_data, qs_env)
      ALLOCATE (iatom_to_subgroup(natom))
      CALL get_sub_dist(iatom_to_subgroup, sparsity_pattern, ngroups, ri_data)

      !Get all the required tensors in the subgroups
      ALLOCATE (t_2c_inv(natom), mat_2c_pot(nimg), ks_t_sub(nspins, nimg), t_2c_ao_tmp(1), &
                ks_t_split(2), t_2c_work(3))
      CALL get_subgroup_2c_tensors(t_2c_inv, mat_2c_pot, t_2c_work, t_2c_ao_tmp, ks_t_split, ks_t_sub, &
                                   group_size, ngroups, para_env, para_env_sub, ri_data)

      ALLOCATE (t_3c_int(nimg), t_3c_apc_sub(nspins, nimg), t_3c_work_2(3), t_3c_work_3(3))
      CALL get_subgroup_3c_tensors(t_3c_int, t_3c_work_2, t_3c_work_3, t_3c_apc, t_3c_apc_sub, &
                                   iatom_to_subgroup, group_size, ngroups, para_env, para_env_sub, ri_data)

      ALLOCATE (batch_ranges(natom + 1))
      batch_ranges(natom + 1) = SIZE(ri_data%bsizes_AO_split) + 1
      iatom = 0
      DO iblk = 1, SIZE(ri_data%bsizes_AO_split)
         IF (idx_to_at_AO(iblk) == iatom + 1) THEN
            iatom = iatom + 1
            batch_ranges(iatom) = iblk
         END IF
      END DO
      CALL dbt_batched_contract_init(t_3c_work_3(1), batch_range_2=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_work_3(2), batch_range_2=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_work_2(1), batch_range_1=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_work_2(2), batch_range_1=batch_ranges)
      DO b_img = 1, nimg
         CALL dbt_batched_contract_init(ks_t_split(1))
         CALL dbt_batched_contract_init(ks_t_split(2))
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (.NOT. sparsity_pattern(iatom, jatom, b_img) == igroup) CYCLE
               pref = 1.0_dp
               IF (iatom == jatom .AND. b_img == 1) pref = 0.5_dp

               !measure the cost of the given i, j, b configuration
               t3 = m_walltime()

               !Get the proper HFX potential 2c integrals (P_i^0|Q_j^b), and multply by RI metric
               CALL timeset(routineN//"_2c", handle2)
               CALL get_ext_2c_int(t_2c_work(1), mat_2c_pot, iatom, jatom, b_img, ri_data, qs_env, &
                                   blacs_env_ext=blacs_env_sub, para_env_ext=para_env_sub)
               CALL dbt_contract(1.0_dp, t_2c_work(1), t_2c_inv(jatom), &
                                 0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[1], notcontract_2=[2], &
                                 filter_eps=ri_data%filter_eps, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_contract(1.0_dp, t_2c_inv(iatom), t_2c_work(2), &
                                 0.0_dp, t_2c_work(1), map_1=[1], map_2=[2], &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[1], notcontract_2=[2], &
                                 filter_eps=ri_data%filter_eps, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_copy(t_2c_work(1), t_2c_work(3), move_data=.TRUE.) !move to split blocks
               CALL dbt_filter(t_2c_work(3), ri_data%filter_eps)
               CALL timestop(handle2)

               CALL dbt_batched_contract_init(t_2c_work(3))
               !Loop over a+c images
               CALL timeset(routineN//"_3c_2", handle2)
               DO i_img = 1, nimg
                  apc_img = get_apc_index_from_ib(i_img, b_img, qs_env)
                  IF (apc_img < 1 .OR. apc_img > nimg) CYCLE

                  !Contract the integrals with the 2c quantity: V_PQ * (Q^b| nu^b lambda^a+c)
                  CALL get_atom_3c_ints(t_3c_work_3(3), t_3c_int(i_img), jatom, 2, idx_to_at_AO)
                  CALL get_tensor_occupancy(t_3c_work_3(3), nze, occ)
                  IF (nze == 0) CYCLE
                  CALL dbt_copy(t_3c_work_3(3), t_3c_work_3(1), move_data=.TRUE.)
                  CALL dbt_contract(1.0_dp, t_2c_work(3), t_3c_work_3(1), &
                                    0.0_dp, t_3c_work_3(2), map_1=[1], map_2=[2, 3], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2, 3], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  CALL dbt_copy(t_3c_work_3(2), t_3c_work_2(2), order=[2, 1, 3], move_data=.TRUE.)
                  CALL dbt_copy(t_3c_work_3(3), t_3c_work_3(1))

                  !Finally get the iatom, jatom block of the Fock matrix
                  DO i_spin = 1, nspins
                     CALL get_atom_3c_ints(t_3c_work_2(3), t_3c_apc_sub(i_spin, apc_img), iatom, 1, idx_to_at_AO)
                     CALL get_tensor_occupancy(t_3c_work_2(3), nze, occ)
                     IF (nze == 0) CYCLE
                     CALL dbt_copy(t_3c_work_2(3), t_3c_work_2(1), move_data=.TRUE.)
                     CALL dbt_contract(-pref*fac, t_3c_work_2(1), t_3c_work_2(2), &
                                       1.0_dp, ks_t_split(i_spin), map_1=[1], map_2=[2], &
                                       contract_1=[2, 3], notcontract_1=[1], &
                                       contract_2=[2, 3], notcontract_2=[1], &
                                       filter_eps=ri_data%filter_eps, &
                                       move_data=i_spin == nspins, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  END DO

               END DO !i_img
               CALL timestop(handle2)
               CALL dbt_batched_contract_finalize(t_2c_work(3))

               t4 = m_walltime()
               ri_data%kp_cost(iatom, jatom, b_img) = t4 - t3
            END DO !iatom
         END DO !jatom
         CALL dbt_batched_contract_finalize(ks_t_split(1))
         CALL dbt_batched_contract_finalize(ks_t_split(2))

         DO i_spin = 1, nspins
            CALL dbt_copy(ks_t_split(i_spin), t_2c_ao_tmp(1), move_data=.TRUE.)
            CALL dbt_copy(t_2c_ao_tmp(1), ks_t_sub(i_spin, b_img), summation=.TRUE.)
         END DO
      END DO !b_img
      CALL dbt_batched_contract_finalize(t_3c_work_3(1))
      CALL dbt_batched_contract_finalize(t_3c_work_3(2))
      CALL dbt_batched_contract_finalize(t_3c_work_2(1))
      CALL dbt_batched_contract_finalize(t_3c_work_2(2))
      CALL para_env%sync()
      CALL para_env%sum(ri_data%dbcsr_nflop)
      CALL para_env%sum(ri_data%kp_cost)
      t2 = m_walltime()
      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      !transfer KS tensor from subgroup to main group
      CALL gather_ks_matrix(ks_t, ks_t_sub, group_size, sparsity_pattern, para_env, ri_data)

      !clean-up subgroup tensors
      CALL dbt_destroy(t_2c_ao_tmp(1))
      CALL dbt_destroy(ks_t_split(1))
      CALL dbt_destroy(ks_t_split(2))
      CALL dbt_destroy(t_2c_work(1))
      CALL dbt_destroy(t_2c_work(2))
      CALL dbt_destroy(t_2c_work(3))
      CALL dbt_destroy(t_3c_work_2(1))
      CALL dbt_destroy(t_3c_work_2(2))
      CALL dbt_destroy(t_3c_work_2(3))
      CALL dbt_destroy(t_3c_work_3(1))
      CALL dbt_destroy(t_3c_work_3(2))
      CALL dbt_destroy(t_3c_work_3(3))
      DO i_img = 1, nimg
         CALL dbt_destroy(t_3c_int(i_img))
         CALL dbcsr_release(mat_2c_pot(i_img))
         DO i_spin = 1, nspins
            CALL dbt_destroy(t_3c_apc(i_spin, i_img))
            CALL dbt_destroy(t_3c_apc_sub(i_spin, i_img))
            CALL dbt_destroy(ks_t_sub(i_spin, i_img))
         END DO
      END DO
      DO iatom = 1, natom
         CALL dbt_destroy(t_2c_inv(iatom))
         DEALLOCATE (iatom_to_subgroup(iatom)%array)
      END DO
      DEALLOCATE (t_2c_ao_tmp)

      !End of subgroup parallelization
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL para_env_sub%free()
      DEALLOCATE (para_env_sub)

      !symmetrize the KS matrices
      ALLOCATE (t_2c_ao_tmp(nimg))
      DO i_spin = 1, nspins
         DO b_img = 1, nimg
            CALL dbt_create(ks_t(1, 1), t_2c_ao_tmp(b_img))
            CALL dbt_copy(ks_t(i_spin, b_img), t_2c_ao_tmp(b_img))
         END DO
         DO b_img = 1, nimg
            mb_img = get_opp_index(b_img, qs_env)
            IF (mb_img > 0 .AND. mb_img .LE. nimg) THEN
               CALL dbt_copy(t_2c_ao_tmp(b_img), ks_t(i_spin, mb_img), order=[2, 1], summation=.TRUE.)
            END IF
            CALL dbt_destroy(t_2c_ao_tmp(b_img))
         END DO
      END DO

      !calculate the energy
      CALL dbt_create(ks_t(1, 1), t_2c_ao_tmp(1))
      CALL dbcsr_create(ks_desymm, template=ks_matrix(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(rho_desymm, template=ks_matrix(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      ehfx = 0.0_dp
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_copy(ks_t(i_spin, i_img), t_2c_ao_tmp(1))
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), ks_desymm)
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), ks_matrix(i_spin, i_img)%matrix)
            CALL dbcsr_filter(ks_matrix(i_spin, i_img)%matrix, ri_data%filter_eps)

            CALL dbt_copy(rho_ao_t(i_spin, i_img), t_2c_ao_tmp(1))
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), rho_desymm)

            CALL dbcsr_dot(ks_desymm, rho_desymm, etmp)
            ehfx = ehfx + 0.5_dp*etmp
         END DO
      END DO
      CALL dbcsr_release(rho_desymm)
      CALL dbcsr_release(ks_desymm)

      CALL dbt_destroy(t_2c_ao_tmp(1))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_destroy(rho_ao_t(i_spin, i_img))
            CALL dbt_destroy(ks_t(i_spin, i_img))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_ks_kp

! **************************************************************************************************
!> \brief return the cell index a+c corresponding to given cell index i and b, with i = a+c-b
!> \param i_index ...
!> \param b_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_apc_index_from_ib(i_index, b_index, qs_env) RESULT(apc_index)
      INTEGER, INTENT(IN)                                :: i_index, b_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: apc_index

      INTEGER, DIMENSION(3)                              :: cell_apc
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      !i = a+c-b => a+c = i+b
      cell_apc(:) = index_to_cell(:, i_index) + index_to_cell(:, b_index)

      IF (ANY([cell_apc(1), cell_apc(2), cell_apc(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([cell_apc(1), cell_apc(2), cell_apc(3)] > UBOUND(cell_to_index))) THEN

         apc_index = 0
      ELSE
         apc_index = cell_to_index(cell_apc(1), cell_apc(2), cell_apc(3))
      END IF

   END FUNCTION get_apc_index_from_ib

! **************************************************************************************************
!> \brief return the cell index i corresponding to the summ of cell_a and cell_c
!> \param a_index ...
!> \param c_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_apc_index(a_index, c_index, qs_env) RESULT(i_index)
      INTEGER, INTENT(IN)                                :: a_index, c_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: i_index

      INTEGER, DIMENSION(3)                              :: cell_i
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      cell_i(:) = index_to_cell(:, a_index) + index_to_cell(:, c_index)

      IF (ANY([cell_i(1), cell_i(2), cell_i(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([cell_i(1), cell_i(2), cell_i(3)] > UBOUND(cell_to_index))) THEN

         i_index = 0
      ELSE
         i_index = cell_to_index(cell_i(1), cell_i(2), cell_i(3))
      END IF

   END FUNCTION get_apc_index

! **************************************************************************************************
!> \brief A routine that returns all allowed a,c pairs such that a+c images corresponds to input
!> \param ac_pairs ...
!> \param apc_index ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_ac_pairs(ac_pairs, apc_index, qs_env)
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: ac_pairs
      INTEGER, INTENT(IN)                                :: apc_index
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: a_index, c_index, i_pair, nimg

      nimg = SIZE(ac_pairs, 2)

      ac_pairs(:, :) = 0
      i_pair = 0
      DO a_index = 1, nimg
         DO c_index = 1, nimg
            IF (get_apc_index(a_index, c_index, qs_env) == apc_index) THEN
               i_pair = i_pair + 1
               ac_pairs(1, i_pair) = a_index
               ac_pairs(2, i_pair) = c_index
            END IF
         END DO
      END DO

   END SUBROUTINE get_ac_pairs

! **************************************************************************************************
!> \brief A function that, given a cell index a, returun the index corresponding to -a, and zero if
!>        if out of bounds
!> \param a_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_opp_index(a_index, qs_env) RESULT(opp_index)
      INTEGER, INTENT(IN)                                :: a_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: opp_index

      INTEGER, DIMENSION(3)                              :: opp_cell
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (kpoints, cell_to_index, index_to_cell)

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      opp_cell(:) = -index_to_cell(:, a_index)

      IF (ANY([opp_cell(1), opp_cell(2), opp_cell(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([opp_cell(1), opp_cell(2), opp_cell(3)] > UBOUND(cell_to_index))) THEN

         opp_index = 0
      ELSE
         opp_index = cell_to_index(opp_cell(1), opp_cell(2), opp_cell(3))
      END IF

   END FUNCTION get_opp_index

! **************************************************************************************************
!> \brief A routine that returns the actual non-symemtric density matrix for each image, by Fourier
!>        transforming the kpoint density matrix
!> \param rho_ao_t ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_pmat_images(rho_ao_t, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: rho_ao_t
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: cell_j(3), i_img, i_spin, iatom, icol, &
                                                            irow, j_img, jatom, mi_img, mj_img, &
                                                            nimg, nspins
      INTEGER(int_8)                                     :: max_nze, nze
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(dp)                                           :: fac, occ
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock, pblock_desymm
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, rho_ao, rho_desymm
      TYPE(dbt_type)                                     :: tmp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl, sab_nl_nosym
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      NULLIFY (rho_desymm, kpoints, sab_nl_nosym, scf_env, fmwork, rho_ao, matrix_ks, rho, sab_nl, &
               nl_iterator, cell_to_index, pblock, pblock_desymm)

      CALL get_qs_env(qs_env, kpoints=kpoints, scf_env=scf_env, matrix_ks_kp=matrix_ks)
      fmwork => scf_env%scf_work1
      CALL get_kpoint_info(kpoints, sab_nl_nosym=sab_nl_nosym, cell_to_index=cell_to_index)

      nspins = SIZE(matrix_ks, 1)
      nimg = ri_data%nimg

      ALLOCATE (rho_desymm(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            ALLOCATE (rho_desymm(i_spin, i_img)%matrix)
            CALL dbcsr_create(rho_desymm(i_spin, i_img)%matrix, template=matrix_ks(i_spin, i_img)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(rho_desymm(i_spin, i_img)%matrix, sab_nl_nosym)
         END DO
      END DO

      !TODO: within the next routine, symtrans is called, and it is not sure it works in all KP symm cases
      CALL kpoint_density_transform(kpoints, rho_desymm, .FALSE., matrix_ks(1, 1)%matrix, &
                                    sab_nl_nosym, fmwork)

      max_nze = 0
      CALL dbt_create(rho_desymm(1, 1)%matrix, tmp)
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_copy_matrix_to_tensor(rho_desymm(i_spin, i_img)%matrix, tmp)
            CALL dbcsr_scale(rho_desymm(i_spin, i_img)%matrix, 0.0_dp)
            CALL dbt_copy(tmp, rho_ao_t(i_spin, i_img), move_data=.TRUE.)
            CALL dbt_filter(rho_ao_t(i_spin, i_img), ri_data%filter_eps)

            CALL get_tensor_occupancy(rho_ao_t(i_spin, i_img), nze, occ)
            IF (nze > max_nze) max_nze = nze
         END DO
      END DO

      !If all densities are empty, then we are at the first SCF step and there is no Kpoint internal
      !density yet. In this case, we construct the non-symmetric pmats from the symmetric ones
      IF (max_nze == 0) THEN
         CALL get_qs_env(qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
         CALL get_kpoint_info(kpoints, sab_nl=sab_nl)

         CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)
            j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (j_img > nimg .OR. j_img < 1) CYCLE

            fac = 1.0_dp
            IF (iatom == jatom) fac = 0.5_dp
            mj_img = get_opp_index(j_img, qs_env)
            !if no opposite image, then no sum of P^j + P^-j => need full diag
            IF (mj_img == 0) fac = 1.0_dp

            irow = iatom
            icol = jatom
            IF (iatom > jatom) THEN
               !because symmetric nl. Value for atom pair i,j is actually stored in j,i if i > j
               irow = jatom
               icol = iatom
            END IF

            DO i_spin = 1, nspins
               CALL dbcsr_get_block_p(rho_ao(i_spin, j_img)%matrix, irow, icol, pblock, found)
               IF (.NOT. found) CYCLE

               !distribution of symm and non-symm matrix match in that way
               CALL dbcsr_get_block_p(rho_desymm(i_spin, j_img)%matrix, iatom, jatom, pblock_desymm, found)
               IF (.NOT. found) CYCLE

               IF (iatom > jatom) THEN
                  pblock_desymm(:, :) = fac*TRANSPOSE(pblock(:, :))
               ELSE
                  pblock_desymm(:, :) = fac*pblock(:, :)
               END IF
            END DO
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         DO i_img = 1, nimg
            DO i_spin = 1, nspins
               CALL dbt_copy_matrix_to_tensor(rho_desymm(i_spin, i_img)%matrix, tmp)
               CALL dbt_copy(tmp, rho_ao_t(i_spin, i_img), summation=.TRUE., move_data=.TRUE.)

               !symmetrize by addin transpose of opp img
               mi_img = get_opp_index(i_img, qs_env)
               IF (mi_img > 0 .AND. mi_img .LE. nimg) THEN
                  CALL dbt_copy_matrix_to_tensor(rho_desymm(i_spin, mi_img)%matrix, tmp)
                  CALL dbt_copy(tmp, rho_ao_t(i_spin, i_img), order=[2, 1], summation=.TRUE., move_data=.TRUE.)
               END IF
               CALL dbt_filter(rho_ao_t(i_spin, i_img), ri_data%filter_eps)
            END DO
         END DO
      END IF

      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbcsr_release(rho_desymm(i_spin, i_img)%matrix)
            DEALLOCATE (rho_desymm(i_spin, i_img)%matrix)
         END DO
      END DO

      CALL dbt_destroy(tmp)
      DEALLOCATE (rho_desymm)

   END SUBROUTINE get_pmat_images

! **************************************************************************************************
!> \brief A routine that, given a cell index b and atom indices ij, returns a 2c tensor with the HFX
!>        potential (P_i^0|Q_j^b), within the extended RI basis
!> \param t_2c_pot ...
!> \param mat_orig ...
!> \param atom_i ...
!> \param atom_j ...
!> \param img_b ...
!> \param ri_data ...
!> \param qs_env ...
!> \param do_inverse ...
!> \param para_env_ext ...
!> \param blacs_env_ext ...
! **************************************************************************************************
   SUBROUTINE get_ext_2c_int(t_2c_pot, mat_orig, atom_i, atom_j, img_b, ri_data, qs_env, do_inverse, &
                             para_env_ext, blacs_env_ext)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_pot
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: mat_orig
      INTEGER, INTENT(IN)                                :: atom_i, atom_j, img_b
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_inverse
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_ext
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_ext_2c_int'

      INTEGER :: blk, group, handle, handle2, i_img, i_RI, iatom, iblk, ikind, img_tot, j_img, &
         j_RI, jatom, jblk, jkind, n_dependent, natom, nblks_RI, nimg, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: present_atoms_i, present_atoms_j
      INTEGER, DIMENSION(3)                              :: cell_b, cell_i, cell_j, cell_tot
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_ext, ri_blk_size_ext, &
                                                            row_dist, row_dist_ext
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell, pgrid
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_inverse_prv, found
      REAL(dp)                                           :: dij, threshold
      REAL(dp), DIMENSION(3)                             :: rij
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist, dbcsr_dist_ext
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter
      TYPE(dbcsr_type)                                   :: work, work_tight, work_tight_inv
      TYPE(dbt_type)                                     :: t_2c_tmp
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, nl_2c, nl_iterator, cell, kpoints, cell_to_index, index_to_cell, dist_2d, &
               para_env, pblock, blacs_env, particle_set, col_dist, row_dist, pgrid, &
               col_dist_ext, row_dist_ext)

      CALL timeset(routineN, handle)

      !Idea: run over the neighbor list once for i and once for j, and record in which cell the MIC
      !      atoms are. Then loop over the atoms and only take the pairs the we need

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, cell=cell, &
                      kpoints=kpoints, para_env=para_env, blacs_env=blacs_env, particle_set=particle_set)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      do_inverse_prv = .FALSE.
      IF (PRESENT(do_inverse)) do_inverse_prv = do_inverse
      IF (do_inverse_prv) THEN
         CPASSERT(atom_i == atom_j)
      END IF

      IF (PRESENT(para_env_ext)) para_env => para_env_ext
      IF (PRESENT(blacs_env_ext)) blacs_env => blacs_env_ext

      nimg = SIZE(mat_orig)

      CALL timeset(routineN//"_nl_iter", handle2)

      !create our own dist_2d in the subgroup
      ALLOCATE (dist1(natom), dist2(natom))
      DO iatom = 1, natom
         dist1(iatom) = MOD(iatom, blacs_env%num_pe(1))
         dist2(iatom) = MOD(iatom, blacs_env%num_pe(2))
      END DO
      CALL distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, blacs_env_ext=blacs_env)

      ALLOCATE (basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%ri_metric, &
                                   "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)

      ALLOCATE (present_atoms_i(natom, nimg), present_atoms_j(natom, nimg))
      present_atoms_i = 0
      present_atoms_j = 0

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=rij, cell=cell_j, &
                                ikind=ikind, jkind=jkind)

         dij = NORM2(rij)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         !TODO: this we could do once and for all for each atom, and carry it around
         !      but only if it is costly. which it probably isn't
         IF (iatom == atom_i .AND. dij .LE. ri_data%kp_RI_range) present_atoms_i(jatom, j_img) = 1
         IF (iatom == atom_j .AND. dij .LE. ri_data%kp_RI_range) present_atoms_j(jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(nl_2c)
      CALL distribution_2d_release(dist_2d)
      CALL timestop(handle2)

      CALL para_env%sum(present_atoms_i)
      CALL para_env%sum(present_atoms_j)

      !Need to build a work matrix with matching distribution to mat_orig
      CALL dbcsr_get_info(mat_orig(1), distribution=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist, row_dist=row_dist, col_dist=col_dist, group=group, pgrid=pgrid)
      ALLOCATE (row_dist_ext(ri_data%ncell_RI*natom), col_dist_ext(ri_data%ncell_RI*natom))
      ALLOCATE (ri_blk_size_ext(ri_data%ncell_RI*natom))
      DO i_RI = 1, ri_data%ncell_RI
         row_dist_ext((i_RI - 1)*natom + 1:i_RI*natom) = row_dist(:)
         col_dist_ext((i_RI - 1)*natom + 1:i_RI*natom) = col_dist(:)
         RI_blk_size_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
      END DO

      CALL dbcsr_distribution_new(dbcsr_dist_ext, group=group, pgrid=pgrid, &
                                  row_dist=row_dist_ext, col_dist=col_dist_ext)
      CALL dbcsr_create(work, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
      CALL dbcsr_distribution_release(dbcsr_dist_ext)
      DEALLOCATE (col_dist_ext, row_dist_ext, RI_blk_size_ext)

      cell_b(:) = index_to_cell(:, img_b)
      DO i_img = 1, nimg
         i_RI = ri_data%img_to_RI_cell(i_img)
         IF (i_RI == 0) CYCLE
         cell_i(:) = index_to_cell(:, i_img)
         DO j_img = 1, nimg
            j_RI = ri_data%img_to_RI_cell(j_img)
            IF (j_RI == 0) CYCLE
            cell_j(:) = index_to_cell(:, j_img)
            cell_tot = cell_j - cell_i + cell_b

            IF (ANY([cell_tot(1), cell_tot(2), cell_tot(3)] < LBOUND(cell_to_index)) .OR. &
                ANY([cell_tot(1), cell_tot(2), cell_tot(3)] > UBOUND(cell_to_index))) CYCLE
            img_tot = cell_to_index(cell_tot(1), cell_tot(2), cell_tot(3))
            IF (img_tot > nimg .OR. img_tot < 1) CYCLE

            CALL dbcsr_iterator_start(dbcsr_iter, mat_orig(img_tot))
            DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
               CALL dbcsr_iterator_next_block(dbcsr_iter, row=iatom, column=jatom, blk=blk)
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               IF (present_atoms_j(jatom, j_img) == 0) CYCLE

               CALL dbcsr_get_block_p(mat_orig(img_tot), iatom, jatom, pblock, found)
               IF (.NOT. found) CYCLE

               CALL dbcsr_put_block(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock)

            END DO
            CALL dbcsr_iterator_stop(dbcsr_iter)

         END DO !j_img
      END DO !i_img
      CALL dbcsr_finalize(work)

      IF (do_inverse_prv) THEN

         !Because there are a lot of empty rows/cols in work, we need to get rid of them for inversion
         nblks_RI = SUM(present_atoms_i)
         ALLOCATE (col_dist_ext(nblks_RI), row_dist_ext(nblks_RI), RI_blk_size_ext(nblks_RI))
         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1
               col_dist_ext(iblk) = col_dist(iatom)
               row_dist_ext(iblk) = row_dist(iatom)
               RI_blk_size_ext(iblk) = ri_data%bsizes_RI(iatom)
            END DO
         END DO

         CALL dbcsr_distribution_new(dbcsr_dist_ext, group=group, pgrid=pgrid, &
                                     row_dist=row_dist_ext, col_dist=col_dist_ext)
         CALL dbcsr_create(work_tight, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
         CALL dbcsr_create(work_tight_inv, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
         CALL dbcsr_distribution_release(dbcsr_dist_ext)
         DEALLOCATE (col_dist_ext, row_dist_ext, RI_blk_size_ext)

         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1

               jblk = 0
               DO j_img = 1, nimg
                  j_RI = ri_data%img_to_RI_cell(j_img)
                  IF (j_RI == 0) CYCLE
                  DO jatom = 1, natom
                     IF (present_atoms_j(jatom, j_img) == 0) CYCLE
                     jblk = jblk + 1

                     CALL dbcsr_get_block_p(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock, found)
                     IF (.NOT. found) CYCLE
                     CALL dbcsr_put_block(work_tight, iblk, jblk, pblock)
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_finalize(work_tight)
         CALL dbcsr_clear(work)

         SELECT CASE (ri_data%t2c_method)
         CASE (hfx_ri_do_2c_iter)
            threshold = MAX(ri_data%filter_eps, 1.0e-12_dp)
            CALL invert_hotelling(work_tight_inv, work_tight, threshold=threshold, silent=.FALSE.)
         CASE (hfx_ri_do_2c_cholesky)
            CALL dbcsr_copy(work_tight_inv, work_tight)
            CALL cp_dbcsr_cholesky_decompose(work_tight_inv, para_env=para_env, blacs_env=blacs_env)
            CALL cp_dbcsr_cholesky_invert(work_tight_inv, para_env=para_env, blacs_env=blacs_env, &
                                          upper_to_full=.TRUE.)
         CASE (hfx_ri_do_2c_diag)
            CALL dbcsr_copy(work_tight_inv, work_tight)
            CALL cp_dbcsr_power(work_tight_inv, -1.0_dp, ri_data%eps_eigval, n_dependent, &
                                para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
         END SELECT

         !move back data to standard extended RI pattern
         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1

               jblk = 0
               DO j_img = 1, nimg
                  j_RI = ri_data%img_to_RI_cell(j_img)
                  IF (j_RI == 0) CYCLE
                  DO jatom = 1, natom
                     IF (present_atoms_j(jatom, j_img) == 0) CYCLE
                     jblk = jblk + 1

                     CALL dbcsr_get_block_p(work_tight_inv, iblk, jblk, pblock, found)
                     IF (.NOT. found) CYCLE
                     CALL dbcsr_put_block(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock)
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_finalize(work)

         CALL dbcsr_release(work_tight)
         CALL dbcsr_release(work_tight_inv)
      END IF

      CALL dbt_create(work, t_2c_tmp)
      CALL dbt_copy_matrix_to_tensor(work, t_2c_tmp)
      CALL dbt_copy(t_2c_tmp, t_2c_pot, move_data=.TRUE.)
      CALL dbt_filter(t_2c_pot, ri_data%filter_eps)

      CALL dbt_destroy(t_2c_tmp)
      CALL dbcsr_release(work)

      CALL timestop(handle)

   END SUBROUTINE get_ext_2c_int

! **************************************************************************************************
!> \brief copies the 3c integrals correspoinding to a single atom mu from the general (P^0| mu^0 sigam^a)
!> \param t_3c_at ...
!> \param t_3c_ints ...
!> \param iatom ...
!> \param dim_at ...
!> \param idx_to_at ...
! **************************************************************************************************
   SUBROUTINE get_atom_3c_ints(t_3c_at, t_3c_ints, iatom, dim_at, idx_to_at)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_at, t_3c_ints
      INTEGER, INTENT(IN)                                :: iatom, dim_at
      INTEGER, DIMENSION(:), INTENT(IN)                  :: idx_to_at

      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbt_iterator_type)                            :: iter

!$OMP PARALLEL DEFAULT(NONE) SHARED(t_3c_ints,t_3c_at,iatom,idx_to_at,dim_at) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t_3c_ints)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         IF (.NOT. idx_to_at(ind(dim_at)) == iatom) CYCLE

         CALL dbt_get_block(t_3c_ints, ind, blk, found)
         IF (.NOT. found) CYCLE

         CALL dbt_put_block(t_3c_at, ind, SHAPE(blk), blk)
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      CALL dbt_finalize(t_3c_at)

   END SUBROUTINE get_atom_3c_ints

! **************************************************************************************************
!> \brief Copy the data of a 2D tensor living in the main MPI group to a sub-group, given the proc
!>        mapping from one to the other (e.g. for a proc idx in the subgroup, we get the idx in the main)
!> \param t2c_sub ...
!> \param t2c_main ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE copy_2c_to_subgroup(t2c_sub, t2c_main, group_size, ngroups, para_env)
      TYPE(dbt_type), INTENT(INOUT)                      :: t2c_sub, t2c_main
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, iblk, igroup, iproc, ir, is, jblk, &
                                                            nocc, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes1, bsizes2, idx1, idx2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: block_dest, block_source
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: current_dest
      INTEGER, DIMENSION(2)                              :: ind, nblks
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      !Stategy: we loop over the main tensor, and send all the data. Then we loop over the sub tensor
      !         and receive it. We do all of it with async MPI communication. The sub tensor needs
      !         to have blocks pre-reserved though

      CALL dbt_get_info(t2c_main, nblks_total=nblks)

      !Loop over the main tensor, count how many blocks are there, which ones, and on which proc
      ALLOCATE (block_source(nblks(1), nblks(2)))
      block_source = -1
      nocc = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(t2c_main,para_env,nocc,block_source) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t2c_main)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(t2c_main, ind, blk, found)
         IF (.NOT. found) CYCLE

         block_source(ind(1), ind(2)) = para_env%mepos
!$OMP ATOMIC
         nocc = nocc + 1
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL para_env%sum(nocc)
      CALL para_env%sum(block_source)
      block_source = block_source + para_env%num_pe - 1

      !Loop over the sub tensor, get the block destination + get blk indices for reservation
      igroup = para_env%mepos/group_size
      ALLOCATE (idx1(nocc), idx2(nocc))
      ALLOCATE (block_dest(nblks(1), nblks(2)))
      block_dest = -1
      nocc = 0
      DO jblk = 1, nblks(2)
         DO iblk = 1, nblks(1)
            IF (block_source(iblk, jblk) == -1) CYCLE
            nocc = nocc + 1
            idx1(nocc) = iblk
            idx2(nocc) = jblk

            CALL dbt_get_stored_coordinates(t2c_sub, [iblk, jblk], iproc)
            block_dest(iblk, jblk) = igroup*group_size + iproc !mapping of iproc in subgroup to main group idx
         END DO
      END DO

      ALLOCATE (bsizes1(nblks(1)), bsizes2(nblks(2)))
      CALL dbt_get_info(t2c_main, blk_size_1=bsizes1, blk_size_2=bsizes2)

      ALLOCATE (current_dest(nblks(1), nblks(2), 0:ngroups - 1))
      DO igroup = 0, ngroups - 1
         !for a given subgroup, need to make the destination available to everyone in the main group
         current_dest(:, :, igroup) = block_dest(:, :)
         CALL para_env%bcast(current_dest(:, :, igroup), source=igroup*group_size) !bcast from first proc in sub-group
      END DO

      !Loop over groups, blocks and send/receive
      ALLOCATE (send_buff(ngroups*nocc), recv_buff(ngroups*nocc))
      ALLOCATE (send_req(ngroups*nocc), recv_req(ngroups*nocc))
      ir = 0
      is = 0
      tag = 0
      DO igroup = 0, ngroups - 1
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               IF (block_source(iblk, jblk) == -1) CYCLE

               !a unique tag per block
               tag = tag + 1

               found = .FALSE.
               IF (para_env%mepos == block_source(iblk, jblk)) THEN
                  CALL dbt_get_block(t2c_main, [iblk, jblk], blk, found)
               END IF

               !If blocks live on same proc, simply copy. Else MPI send/recv
               IF (block_source(iblk, jblk) == current_dest(iblk, jblk, igroup)) THEN
                  IF (found) CALL dbt_put_block(t2c_sub, [iblk, jblk], SHAPE(blk), blk)
               ELSE
                  IF (para_env%mepos == block_source(iblk, jblk) .AND. found) THEN
                     ALLOCATE (send_buff(tag)%array(bsizes1(iblk), bsizes2(jblk)))
                     send_buff(tag)%array(:, :) = blk(:, :)
                     is = is + 1
                     CALL para_env%isend(msgin=send_buff(tag)%array, dest=current_dest(iblk, jblk, igroup), &
                                         request=send_req(is), tag=tag)
                  END IF

                  IF (para_env%mepos == current_dest(iblk, jblk, igroup)) THEN
                     ALLOCATE (recv_buff(tag)%array(bsizes1(iblk), bsizes2(jblk)))
                     ir = ir + 1
                     CALL para_env%irecv(msgout=recv_buff(tag)%array, source=block_source(iblk, jblk), &
                                         request=recv_req(ir), tag=tag)
                  END IF
               END IF

               IF (found) DEALLOCATE (blk)
            END DO
         END DO
      END DO

      CALL mp_waitall(send_req(1:is))
      CALL mp_waitall(recv_req(1:ir))
      !clean-up
      DO i = 1, nocc*ngroups
         IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
      END DO

      !Finally copy the data from the buffer to the sub-tensor
      tag = 0
      DO igroup = 0, ngroups - 1
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               IF (block_source(iblk, jblk) == -1) CYCLE
               tag = tag + 1
               IF (para_env%mepos == current_dest(iblk, jblk, igroup) .AND. &
                   block_source(iblk, jblk) .NE. current_dest(iblk, jblk, igroup)) THEN

                  ALLOCATE (blk(bsizes1(iblk), bsizes2(jblk)))
                  blk(:, :) = recv_buff(tag)%array(:, :)
                  CALL dbt_put_block(t2c_sub, [iblk, jblk], SHAPE(blk), blk)
                  DEALLOCATE (blk)
               END IF
            END DO
         END DO
      END DO

      !clean-up
      DO i = 1, nocc*ngroups
         IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
      END DO
      DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      CALL dbt_finalize(t2c_sub)

   END SUBROUTINE copy_2c_to_subgroup

! **************************************************************************************************
!> \brief Copy the data of a 3D tensor living in the main MPI group to a sub-group, given the proc
!>        mapping from one to the other (e.g. for a proc idx in the subgroup, we get the idx in the main)
!> \param t3c_sub ...
!> \param t3c_main ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param iatom_to_subgroup ...
!> \param dim_at ...
! **************************************************************************************************
   SUBROUTINE copy_3c_to_subgroup(t3c_sub, t3c_main, group_size, ngroups, para_env, iatom_to_subgroup, dim_at)
      TYPE(dbt_type), INTENT(INOUT)                      :: t3c_sub, t3c_main
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_1d_logical_p_type), DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: iatom_to_subgroup
      INTEGER, INTENT(IN), OPTIONAL                      :: dim_at

      INTEGER                                            :: i, iatom, iblk, igroup, iproc, ir, is, &
                                                            jblk, kblk, nocc, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes1, bsizes2, bsizes3, idx1, idx2, &
                                                            idx3
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: block_dest, block_source
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: current_dest
      INTEGER, DIMENSION(3)                              :: ind, nblks
      LOGICAL                                            :: filter_at, found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      !Stategy: we loop over the main tensor, and send all the data. Then we loop over the sub tensor
      !         and receive it. We do all of it with async MPI communication. The sub tensor needs
      !         to have blocks pre-reserved though

      CALL dbt_get_info(t3c_main, nblks_total=nblks)

      !in some cases, only copy a fraction of the 3c tensor to a given subgroup (corresponding to some atoms)
      filter_at = .FALSE.
      IF (PRESENT(iatom_to_subgroup) .AND. PRESENT(dim_at)) THEN
         filter_at = .TRUE.
         CPASSERT(nblks(dim_at) == SIZE(iatom_to_subgroup))
      END IF

      !Loop over the main tensor, count how many blocks are there, which ones, and on which proc
      ALLOCATE (block_source(nblks(1), nblks(2), nblks(3)))
      block_source = -1
      nocc = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(t3c_main,para_env,nocc,block_source) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t3c_main)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(t3c_main, ind, blk, found)
         IF (.NOT. found) CYCLE

         block_source(ind(1), ind(2), ind(3)) = para_env%mepos
!$OMP ATOMIC
         nocc = nocc + 1
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL para_env%sum(nocc)
      CALL para_env%sum(block_source)
      block_source = block_source + para_env%num_pe - 1

      !Loop over the sub tensor, get the block destination + get blk indices for reservation
      igroup = para_env%mepos/group_size
      ALLOCATE (idx1(nocc), idx2(nocc), idx3(nocc))
      ALLOCATE (block_dest(nblks(1), nblks(2), nblks(3)))
      block_dest = -1
      nocc = 0
      DO kblk = 1, nblks(3)
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               IF (block_source(iblk, jblk, kblk) == -1) CYCLE
               nocc = nocc + 1
               idx1(nocc) = iblk
               idx2(nocc) = jblk
               idx3(nocc) = kblk

               CALL dbt_get_stored_coordinates(t3c_sub, [iblk, jblk, kblk], iproc)
               block_dest(iblk, jblk, kblk) = igroup*group_size + iproc !mapping of iproc in subgroup to main group idx
            END DO
         END DO
      END DO

      ALLOCATE (bsizes1(nblks(1)), bsizes2(nblks(2)), bsizes3(nblks(3)))
      CALL dbt_get_info(t3c_main, blk_size_1=bsizes1, blk_size_2=bsizes2, blk_size_3=bsizes3)

      ALLOCATE (current_dest(nblks(1), nblks(2), nblks(3), 0:ngroups - 1))
      DO igroup = 0, ngroups - 1
         !for a given subgroup, need to make the destination available to everyone in the main group
         current_dest(:, :, :, igroup) = block_dest(:, :, :)
         CALL para_env%bcast(current_dest(:, :, :, igroup), source=igroup*group_size) !bcast from first proc in subgroup
      END DO

      !Loop over groups, blocks and send/receive
      ALLOCATE (send_buff(ngroups*nocc), recv_buff(ngroups*nocc))
      ALLOCATE (send_req(ngroups*nocc), recv_req(ngroups*nocc))
      ir = 0
      is = 0
      tag = 0
      DO igroup = 0, ngroups - 1
         DO kblk = 1, nblks(3)
            DO jblk = 1, nblks(2)
               DO iblk = 1, nblks(1)
                  IF (block_source(iblk, jblk, kblk) == -1) CYCLE

                  !a unique tag per block
                  tag = tag + 1

                  IF (filter_at) THEN
                     ind(:) = [iblk, jblk, kblk]
                     iatom = ind(dim_at)
                     IF (.NOT. iatom_to_subgroup(iatom)%array(igroup + 1)) CYCLE
                  END IF

                  found = .FALSE.
                  IF (para_env%mepos == block_source(iblk, jblk, kblk)) THEN
                     CALL dbt_get_block(t3c_main, [iblk, jblk, kblk], blk, found)
                  END IF

                  !If blocks live on same proc, simply copy. Else MPI send/recv
                  IF (block_source(iblk, jblk, kblk) == current_dest(iblk, jblk, kblk, igroup)) THEN
                     IF (found) CALL dbt_put_block(t3c_sub, [iblk, jblk, kblk], SHAPE(blk), blk)
                  ELSE
                     IF (para_env%mepos == block_source(iblk, jblk, kblk) .AND. found) THEN
                        ALLOCATE (send_buff(tag)%array(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                        send_buff(tag)%array(:, :, :) = blk(:, :, :)
                        is = is + 1
                        CALL para_env%isend(msgin=send_buff(tag)%array, &
                                            dest=current_dest(iblk, jblk, kblk, igroup), &
                                            request=send_req(is), tag=tag)
                     END IF

                     IF (para_env%mepos == current_dest(iblk, jblk, kblk, igroup)) THEN
                        ALLOCATE (recv_buff(tag)%array(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                        ir = ir + 1
                        CALL para_env%irecv(msgout=recv_buff(tag)%array, source=block_source(iblk, jblk, kblk), &
                                            request=recv_req(ir), tag=tag)
                     END IF
                  END IF

                  IF (found) DEALLOCATE (blk)
               END DO
            END DO
         END DO
      END DO

      CALL mp_waitall(send_req(1:is))
      CALL mp_waitall(recv_req(1:ir))
      !clean-up
      DO i = 1, nocc*ngroups
         IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
      END DO

      !Finally copy the data from the buffer to the sub-tensor
      tag = 0
      DO igroup = 0, ngroups - 1
         DO kblk = 1, nblks(3)
            DO jblk = 1, nblks(2)
               DO iblk = 1, nblks(1)
                  IF (block_source(iblk, jblk, kblk) == -1) CYCLE
                  tag = tag + 1
                  IF (filter_at) THEN
                     ind(:) = [iblk, jblk, kblk]
                     iatom = ind(dim_at)
                     IF (.NOT. iatom_to_subgroup(iatom)%array(igroup + 1)) CYCLE
                  END IF
                  IF (para_env%mepos == current_dest(iblk, jblk, kblk, igroup) .AND. &
                      block_source(iblk, jblk, kblk) .NE. current_dest(iblk, jblk, kblk, igroup)) THEN

                     ALLOCATE (blk(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                     blk(:, :, :) = recv_buff(tag)%array(:, :, :)
                     CALL dbt_put_block(t3c_sub, [iblk, jblk, kblk], SHAPE(blk), blk)
                     DEALLOCATE (blk)
                  END IF
               END DO
            END DO
         END DO
      END DO

      !clean-up
      DO i = 1, nocc*ngroups
         IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
      END DO
      DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      CALL dbt_finalize(t3c_sub)

   END SUBROUTINE copy_3c_to_subgroup

! **************************************************************************************************
!> \brief A routine that gather the pieces of the KS matrix accross the subgroup and puts it in the
!>        main group. Each b_img, iatom, jatom tuple is one a single CPU
!> \param ks_t ...
!> \param ks_t_sub ...
!> \param group_size ...
!> \param sparsity_pattern ...
!> \param para_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE gather_ks_matrix(ks_t, ks_t_sub, group_size, sparsity_pattern, para_env, ri_data)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: ks_t, ks_t_sub
      INTEGER, INTENT(IN)                                :: group_size
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: sparsity_pattern
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gather_ks_matrix'

      INTEGER                                            :: b_img, dest, handle, i, i_spin, iatom, &
                                                            igroup, ir, is, jatom, n_mess, natom, &
                                                            nimg, nspins, source, tag
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      CALL timeset(routineN, handle)

      nimg = SIZE(sparsity_pattern, 3)
      natom = SIZE(sparsity_pattern, 2)
      nspins = SIZE(ks_t, 1)

      DO b_img = 1, nimg
         n_mess = 0
         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) > -1) n_mess = n_mess + 1
               END DO
            END DO
         END DO

         ALLOCATE (send_buff(n_mess), recv_buff(n_mess))
         ALLOCATE (send_req(n_mess), recv_req(n_mess))
         ir = 0
         is = 0
         n_mess = 0
         tag = 0

         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) < 0) CYCLE
                  n_mess = n_mess + 1
                  tag = tag + 1

                  !sending the message
                  CALL dbt_get_stored_coordinates(ks_t(i_spin, b_img), [iatom, jatom], dest)
                  CALL dbt_get_stored_coordinates(ks_t_sub(i_spin, b_img), [iatom, jatom], source) !source within sub
                  igroup = sparsity_pattern(iatom, jatom, b_img)
                  source = source + igroup*group_size
                  IF (para_env%mepos == source) THEN
                     CALL dbt_get_block(ks_t_sub(i_spin, b_img), [iatom, jatom], blk, found)
                     IF (source == dest) THEN
                        IF (found) CALL dbt_put_block(ks_t(i_spin, b_img), [iatom, jatom], SHAPE(blk), blk)
                     ELSE
                        ALLOCATE (send_buff(n_mess)%array(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                        send_buff(n_mess)%array(:, :) = 0.0_dp
                        IF (found) THEN
                           send_buff(n_mess)%array(:, :) = blk(:, :)
                        END IF
                        is = is + 1
                        CALL para_env%isend(msgin=send_buff(n_mess)%array, dest=dest, &
                                            request=send_req(is), tag=tag)
                     END IF
                     DEALLOCATE (blk)
                  END IF

                  !receiving the message
                  IF (para_env%mepos == dest .AND. source .NE. dest) THEN
                     ALLOCATE (recv_buff(n_mess)%array(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                     ir = ir + 1
                     CALL para_env%irecv(msgout=recv_buff(n_mess)%array, source=source, &
                                         request=recv_req(ir), tag=tag)
                  END IF
               END DO !iatom
            END DO !jatom
         END DO !ispin

         CALL mp_waitall(send_req(1:is))
         CALL mp_waitall(recv_req(1:ir))

         !Copy the messages received into the KS matrix
         n_mess = 0
         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) < 0) CYCLE
                  n_mess = n_mess + 1

                  CALL dbt_get_stored_coordinates(ks_t(i_spin, b_img), [iatom, jatom], dest)
                  IF (para_env%mepos == dest) THEN
                     IF (.NOT. ASSOCIATED(recv_buff(n_mess)%array)) CYCLE
                     ALLOCATE (blk(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                     blk(:, :) = recv_buff(n_mess)%array(:, :)
                     CALL dbt_put_block(ks_t(i_spin, b_img), [iatom, jatom], SHAPE(blk), blk)
                     DEALLOCATE (blk)
                  END IF
               END DO
            END DO
         END DO

         !clean-up
         DO i = 1, n_mess
            IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
            IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
         END DO
         DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      END DO !b_img

      CALL timestop(handle)

   END SUBROUTINE gather_ks_matrix

! **************************************************************************************************
!> \brief copy all required 2c tensors from the main MPI group to the subgroups
!> \param t_2c_inv ...
!> \param mat_2c_pot ...
!> \param t_2c_work ...
!> \param t_2c_ao_tmp ...
!> \param ks_t_split ...
!> \param ks_t_sub ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_subgroup_2c_tensors(t_2c_inv, mat_2c_pot, t_2c_work, t_2c_ao_tmp, ks_t_split, &
                                      ks_t_sub, group_size, ngroups, para_env, para_env_sub, ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_inv
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: mat_2c_pot
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_work, t_2c_ao_tmp, ks_t_split
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: ks_t_sub
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_2c_tensors'

      INTEGER                                            :: handle, i, i_img, i_RI, i_spin, iatom, &
                                                            iproc, j, natom, nblks, nimg, nspins
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            dist1, dist2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, RI_blk_size, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: dbcsr_pgrid
      REAL(dp)                                           :: occ
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist_sub
      TYPE(dbt_pgrid_type)                               :: pgrid_2d
      TYPE(dbt_type)                                     :: work, work_sub

      CALL timeset(routineN, handle)

      !Create the 2d pgrid
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env_sub, pdims_2d, pgrid_2d)

      natom = SIZE(ri_data%bsizes_RI)
      nblks = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*natom))
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks + 1:i_RI*nblks) = ri_data%bsizes_RI_split(:)
      END DO

      !nRI x nRI 2c tensors
      CALL create_2c_tensor(t_2c_inv(1), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      DO iatom = 2, natom
         CALL dbt_create(t_2c_inv(1), t_2c_inv(iatom))
      END DO
      CALL dbt_create(t_2c_inv(1), t_2c_work(1))
      CALL dbt_create(t_2c_inv(1), t_2c_work(2))

      CALL create_2c_tensor(t_2c_work(3), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext_split, bsizes_RI_ext_split, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      !cpoy the data from the main group.
      DO iatom = 1, natom
         CALL copy_2c_to_subgroup(t_2c_inv(iatom), ri_data%t_2c_inv(1, iatom), group_size, ngroups, para_env)
      END DO

      !the AO based tensors
      CALL create_2c_tensor(ks_t_split(1), dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(ks_t_split(1), ks_t_split(2))

      CALL create_2c_tensor(ks_t_sub(1, 1), dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_AO, ri_data%bsizes_AO, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(ks_t_sub(1, 1), t_2c_ao_tmp(1))

      nspins = SIZE(ks_t_sub, 1)
      nimg = SIZE(ks_t_sub, 2)
      IF (nspins == 2) CALL dbt_create(ks_t_sub(1, 1), ks_t_sub(2, 1))
      DO i_img = 2, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ks_t_sub(1, 1), ks_t_sub(i_spin, i_img))
         END DO
      END DO

      !Finally the HFX potential matrices
      !For now, we do a convoluted things where we go to tensors first, then back to matrices.
      CALL create_2c_tensor(work_sub, dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_RI, ri_data%bsizes_RI, &
                            name="(RI | RI)")
      CALL dbt_create(ri_data%mat_2c_pot(1, 1), work)

      ALLOCATE (dbcsr_pgrid(0:pdims_2d(1) - 1, 0:pdims_2d(2) - 1))
      iproc = 0
      DO i = 0, pdims_2d(1) - 1
         DO j = 0, pdims_2d(2) - 1
            dbcsr_pgrid(i, j) = iproc
            iproc = iproc + 1
         END DO
      END DO

      !We need to have the same exact 2d block dist as the tensors
      ALLOCATE (col_dist(natom), row_dist(natom))
      row_dist(:) = dist1(:)
      col_dist(:) = dist2(:)

      ALLOCATE (RI_blk_size(natom))
      RI_blk_size(:) = ri_data%bsizes_RI(:)

      CALL dbcsr_distribution_new(dbcsr_dist_sub, group=para_env_sub%get_handle(), pgrid=dbcsr_pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(mat_2c_pot(1), dist=dbcsr_dist_sub, name="sub", matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=RI_blk_size, col_blk_size=RI_blk_size)

      DO i_img = 1, nimg
         IF (i_img > 1) CALL dbcsr_create(mat_2c_pot(i_img), template=mat_2c_pot(1))
         CALL dbt_copy_matrix_to_tensor(ri_data%mat_2c_pot(1, i_img), work)
         CALL get_tensor_occupancy(work, nze, occ)
         IF (nze == 0) CYCLE

         CALL copy_2c_to_subgroup(work_sub, work, group_size, ngroups, para_env)
         CALL dbt_copy_tensor_to_matrix(work_sub, mat_2c_pot(i_img))
         CALL dbcsr_filter(mat_2c_pot(i_img), ri_data%filter_eps)
         CALL dbt_clear(work_sub)
      END DO

      CALL dbt_destroy(work)
      CALL dbt_destroy(work_sub)
      CALL dbt_pgrid_destroy(pgrid_2d)
      CALL dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (col_dist, row_dist, RI_blk_size, dbcsr_pgrid)
      CALL timestop(handle)

   END SUBROUTINE get_subgroup_2c_tensors

! **************************************************************************************************
!> \brief copy all required 3c tensors from the main MPI group to the subgroups
!> \param t_3c_int ...
!> \param t_3c_work_2 ...
!> \param t_3c_work_3 ...
!> \param t_3c_apc ...
!> \param t_3c_apc_sub ...
!> \param iatom_to_subgroup ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_subgroup_3c_tensors(t_3c_int, t_3c_work_2, t_3c_work_3, t_3c_apc, t_3c_apc_sub, &
                                      iatom_to_subgroup, group_size, ngroups, para_env, para_env_sub, ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_int, t_3c_work_2, t_3c_work_3
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_apc, t_3c_apc_sub
      TYPE(cp_1d_logical_p_type), DIMENSION(:)           :: iatom_to_subgroup
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_3c_tensors'

      INTEGER                                            :: handle, handle2, i_img, i_RI, i_spin, &
                                                            natom, nblks, nimg, nspins
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            dist1, dist2, dist3
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pdims
      REAL(dp)                                           :: occ
      TYPE(dbt_pgrid_type)                               :: pgrid, pgrid_2d
      TYPE(dbt_type)                                     :: tmp, tmp_sub

      CALL timeset(routineN, handle)

      natom = SIZE(ri_data%bsizes_RI)
      nblks = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*natom))
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks + 1:i_RI*nblks) = ri_data%bsizes_RI_split(:)
      END DO

      !Create the pgrid for the configuration correspoinding to ri_data%t_3c_int_ctr_3
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env_sub, pdims_2d, pgrid_2d, &
                            tensor_dims=[SIZE(bsizes_RI_ext_split), SIZE(ri_data%bsizes_AO_split)])
      CALL dbt_pgrid_destroy(pgrid_2d)

      pdims(1) = pdims_2d(1)
      pdims(2) = 1
      pdims(3) = pdims_2d(2)
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid)

      !Utility atomic size tensors for less messages sent
      CALL create_3c_tensor(tmp, dist1, dist2, dist3, &
                            ri_data%pgrid, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(tmp_sub, dist1, dist2, dist3, &
                            pgrid, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      !Create all required 3c tensors in that configuration
      nimg = SIZE(t_3c_int)

      CALL create_3c_tensor(t_3c_int(1), dist1, dist2, dist3, &
                            pgrid, bsizes_RI_ext_split, ri_data%bsizes_AO_split, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      DO i_img = 2, nimg
         CALL dbt_create(t_3c_int(1), t_3c_int(i_img))
      END DO
      CALL dbt_create(t_3c_int(1), t_3c_work_3(1))
      CALL dbt_create(t_3c_int(1), t_3c_work_3(2))
      CALL dbt_create(t_3c_int(1), t_3c_work_3(3))

      !Finally copy the integrals into the subgroups
      CALL timeset(routineN//"_ints", handle2)
      DO i_img = 1, nimg
         CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_1(1, i_img), nze, occ)
         IF (nze == 0) CYCLE
         CALL dbt_copy(ri_data%t_3c_int_ctr_1(1, i_img), tmp, order=[2, 1, 3])
         CALL copy_3c_to_subgroup(tmp_sub, tmp, group_size, ngroups, para_env)
         CALL dbt_copy(tmp_sub, t_3c_int(i_img), move_data=.TRUE.)
         CALL dbt_filter(t_3c_int(i_img), ri_data%filter_eps)
      END DO
      CALL timestop(handle2)
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_destroy(tmp)
      CALL dbt_destroy(tmp_sub)

      !Do the same for the t_3c_ctr_2 configuration
      pdims(1) = 1
      pdims(2) = pdims_2d(1)
      pdims(3) = pdims_2d(2)
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid)

      CALL create_3c_tensor(t_3c_work_2(1), dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO_split, bsizes_RI_ext_split, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(2))
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(3))

      !Utility atomic size tensors for less messages sent
      CALL create_3c_tensor(tmp, dist1, dist2, dist3, &
                            ri_data%pgrid, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO | RI AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(tmp_sub, dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO | RI AO)")
      DEALLOCATE (dist1, dist2, dist3)

      !Finally copy data from t_3c_apc to the subgroups
      nspins = SIZE(t_3c_apc, 1)
      CALL timeset(routineN//"_apc", handle2)
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            !TODO: necessary to switch to atomic block sizes? Could it be faster otherwise?
            CALL dbt_create(t_3c_work_2(1), t_3c_apc_sub(i_spin, i_img))
            CALL get_tensor_occupancy(t_3c_apc(i_spin, i_img), nze, occ)
            IF (nze == 0) CYCLE
            CALL dbt_copy(t_3c_apc(i_spin, i_img), tmp, move_data=.TRUE.)
            CALL copy_3c_to_subgroup(tmp_sub, tmp, group_size, ngroups, para_env, iatom_to_subgroup, 1)
            CALL dbt_copy(tmp_sub, t_3c_apc_sub(i_spin, i_img), move_data=.TRUE.)
            CALL dbt_filter(t_3c_apc_sub(i_spin, i_img), ri_data%filter_eps)
         END DO
      END DO
      CALL timestop(handle2)
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_destroy(tmp)
      CALL dbt_destroy(tmp_sub)

      CALL timestop(handle)

   END SUBROUTINE get_subgroup_3c_tensors

! **************************************************************************************************
!> \brief Get the sparsity pattern related to the non-symmetric AO basis overlap neighbor list
!> \param pattern ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_sparsity_pattern(pattern, ri_data, qs_env)
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: pattern
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: iatom, j_img, jatom, mj_img, natom, nimg
      INTEGER, DIMENSION(3)                              :: cell_j
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c

      NULLIFY (nl_2c, nl_iterator, kpoints, cell_to_index, dft_control, index_to_cell, para_env)

      CALL get_qs_env(qs_env, kpoints=kpoints, dft_control=dft_control, para_env=para_env, natom=natom)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell, sab_nl=nl_2c)

      nimg = ri_data%nimg
      pattern(:, :, :) = 0

      !We use the symmetric nl for all images that have an opposite cell
      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         mj_img = get_opp_index(j_img, qs_env)
         IF (mj_img > nimg .OR. mj_img < 1) CYCLE

         IF (ri_data%present_images(j_img) == 0) CYCLE

         pattern(iatom, jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      !If there is no opposite cell present, then we take into account the non-symmetric nl
      CALL get_kpoint_info(kpoints, sab_nl_nosym=nl_2c)

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         mj_img = get_opp_index(j_img, qs_env)
         IF (mj_img .LE. nimg .AND. mj_img > 0) CYCLE

         IF (ri_data%present_images(j_img) == 0) CYCLE

         pattern(iatom, jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL para_env%sum(pattern)

      !If the opposite image is considered, then there is no need to compute diagonal twice
      DO j_img = 2, nimg
         DO iatom = 1, natom
            IF (pattern(iatom, iatom, j_img) .NE. 0) THEN
               mj_img = get_opp_index(j_img, qs_env)
               IF (mj_img > nimg .OR. mj_img < 1) CYCLE
               pattern(iatom, iatom, mj_img) = 0
            END IF
         END DO
      END DO

      ! -1 => unoccupied, 0 => occupied
      pattern(:, :, :) = pattern(:, :, :) - 1

   END SUBROUTINE get_sparsity_pattern

! **************************************************************************************************
!> \brief Distribute the iatom, jatom, b_img triplet over the subgroupd to spread the load
!>        the group id for each triplet is passed as the value of sparsity_pattern(i, j, b),
!>        with -1 being an unoccupied block
!> \param iatom_to_subgroup ...
!> \param sparsity_pattern ...
!> \param ngroups ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_sub_dist(iatom_to_subgroup, sparsity_pattern, ngroups, ri_data)
      TYPE(cp_1d_logical_p_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: iatom_to_subgroup
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: sparsity_pattern
      INTEGER, INTENT(IN)                                :: ngroups
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: b_img, i, iatom, igroup, jatom, lb, &
                                                            natom, nimg
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: max_group_per_at
      REAL(dp)                                           :: cost, cost_tot, fac
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: bins, cost_per_at

      natom = SIZE(sparsity_pattern, 2)
      nimg = SIZE(sparsity_pattern, 3)

      ALLOCATE (bins(ngroups))
      bins = 0.0_dp

      !To avoid unnecessary data replication accross the subgroups, we want to have a limited number
      !of subgroup with the data of a given iatom. At the minimum, all groups have 1 atom
      ALLOCATE (cost_per_at(natom))
      cost_per_at = 0.0_dp
      cost_tot = 0.0_dp
      DO b_img = 1, nimg
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (sparsity_pattern(iatom, jatom, b_img) == -1) CYCLE
               cost = REAL(ri_data%bsizes_AO(iatom)*ri_data%bsizes_AO(jatom), dp)
               cost_per_at(iatom) = cost_per_at(iatom) + cost
               cost_tot = cost_tot + cost
            END DO
         END DO
      END DO

      ALLOCATE (max_group_per_at(natom))
      lb = MAX(1, CEILING(REAL(ngroups, dp)/REAL(natom, dp)))
      max_group_per_at = MAX(lb, ngroups/ri_data%n_mem)
      DO iatom = 1, natom
         fac = cost_per_at(iatom)/cost_tot*REAL(natom, dp)
         max_group_per_at(iatom) = MIN(ngroups, CEILING(REAL(max_group_per_at(iatom), dp)*fac))
      END DO

      DO iatom = 1, natom
         NULLIFY (iatom_to_subgroup(iatom)%array)
         ALLOCATE (iatom_to_subgroup(iatom)%array(ngroups))
         iatom_to_subgroup(iatom)%array(:) = .FALSE.
         DO i = 1, max_group_per_at(iatom)
            igroup = MINLOC(bins, 1) - 1
            bins(igroup + 1) = bins(igroup + 1) + 1.0_dp
            iatom_to_subgroup(iatom)%array(igroup + 1) = .TRUE.
         END DO
      END DO

      bins = 0.0_dp
      DO b_img = 1, nimg
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (sparsity_pattern(iatom, jatom, b_img) == -1) CYCLE
               igroup = MINLOC(bins, 1, MASK=iatom_to_subgroup(iatom)%array) - 1

               !Use cost information from previous SCF if available
               IF (ri_data%kp_cost(iatom, jatom, b_img) > EPSILON(0.0_dp)) THEN
                  cost = ri_data%kp_cost(iatom, jatom, b_img)
               ELSE
                  cost = REAL(ri_data%bsizes_AO(iatom)*ri_data%bsizes_AO(jatom), dp)
               END IF
               bins(igroup + 1) = bins(igroup + 1) + cost
               sparsity_pattern(iatom, jatom, b_img) = igroup
            END DO
         END DO
      END DO

   END SUBROUTINE get_sub_dist

! **************************************************************************************************
!> \brief A routine that determines the extend of the KP RI-HFX periodic images, including for the
!>        extension of the RI basis
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_kp_and_ri_images(ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kp_and_ri_images'

      INTEGER :: cell_j(3), cell_k(3), handle, i_img, iatom, ikind, j_img, jatom, jcell, katom, &
         kcell, kp_index_lbounds(3), kp_index_ubounds(3), natom, nimg, nkind, pcoord(3), pdims(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_AO_1, dist_AO_2, dist_RI, &
                                                            nRI_per_atom, present_img, RI_cells
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(dp)                                           :: dij, dik, image_range, RI_range, rij(3), &
                                                            rik(3)
      REAL(dp), DIMENSION(:, :), POINTER                 :: pgf_radius
      TYPE(dbt_type)                                     :: t_dummy
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_cart_type)                                 :: mp_comm_t3c
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, dist_2d, nl_2c, nl_iterator, dft_control, particle_set, kpoints, &
               para_env, cell_to_index, pgf_radius)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, distribution_2d=dist_2d, &
                      dft_control=dft_control, particle_set=particle_set, kpoints=kpoints, &
                      para_env=para_env, natom=natom)
      nimg = dft_control%nimages
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      kp_index_lbounds = LBOUND(cell_to_index)
      kp_index_ubounds = UBOUND(cell_to_index)

      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)

      !Determine the range for contributing periodic images, and for the RI basis extension
      ri_data%kp_RI_range = 0.0_dp
      ri_data%kp_image_range = 0.0_dp
      DO ikind = 1, nkind
         !TODO: experiment with reducing the RI range (not only in reproducing supercell energy, but also rel energies)
         CALL init_interaction_radii_orb_basis(basis_set_AO(ikind)%gto_basis_set, ri_data%eps_pgf_orb)
         CALL get_gto_basis_set(basis_set_AO(ikind)%gto_basis_set, kind_radius=RI_range, pgf_radius=pgf_radius)
         ri_data%kp_RI_range = MAX(RI_range, ri_data%kp_RI_range)

         CALL init_interaction_radii_orb_basis(basis_set_RI(ikind)%gto_basis_set, ri_data%eps_pgf_orb)
         CALL get_gto_basis_set(basis_set_RI(ikind)%gto_basis_set, kind_radius=image_range)
         image_range = 2.0_dp*image_range + cutoff_screen_factor*ri_data%hfx_pot%cutoff_radius
         ri_data%kp_image_range = MAX(image_range, ri_data%kp_image_range)
      END DO

      !For the extent of the KP RI-HFX images, we are limited by the RI-HFX potential in
      !(mu^0 sigma^a|P^0) (P^0|Q^b) (Q^b|nu^b lambda^a+c), if there is no contact between
      !any P^0 and Q^b, then image b does not contribute
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)

      ALLOCATE (present_img(nimg))
      present_img = 0
      ri_data%nimg = 0
      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, r=rij, cell=cell_j)

         dij = NORM2(rij)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         IF (dij > ri_data%kp_image_range) CYCLE

         ri_data%nimg = MAX(j_img, ri_data%nimg)
         present_img(j_img) = 1

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(nl_2c)
      CALL para_env%max(ri_data%nimg)
      IF (ri_data%nimg > nimg) &
         CPABORT("Make sure the smallest exponent of the RI-HFX basis is larger than that of the ORB basis.")

      !Keep track of which images will not contribute, so that can be ignored before calculation
      CALL para_env%sum(present_img)
      ALLOCATE (ri_data%present_images(ri_data%nimg))
      ri_data%present_images = 0
      DO i_img = 1, ri_data%nimg
         IF (present_img(i_img) > 0) ri_data%present_images(i_img) = 1
      END DO

      CALL create_3c_tensor(t_dummy, dist_AO_1, dist_AO_2, dist_RI, &
                            ri_data%pgrid, ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI, &
                            map1=[1, 2], map2=[3], name="(AO AO | RI)")

      CALL dbt_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)
      CALL mp_comm_t3c%create(ri_data%pgrid%mp_comm_2d, 3, pdims)
      CALL distribution_3d_create(dist_3d, dist_AO_1, dist_AO_2, dist_RI, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
      CALL dbt_destroy(t_dummy)

      !For the extension of the RI basis P in (mu^0 sigma^a |P^i), we consider an atom if the distance,
      !between mu^0 and P^i if smaller or equal to the kind radius of mu^0
      CALL build_3c_neighbor_lists(nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, dist_3d, &
                                   ri_data%ri_metric, "HFX_3c_nl", qs_env, op_pos=2, sym_ij=.FALSE., &
                                   own_dist=.TRUE.)

      ALLOCATE (RI_cells(nimg))
      RI_cells = 0

      ALLOCATE (nRI_per_atom(natom))
      nRI_per_atom = 0

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, cell_k=cell_k, rik=rik, cell_j=cell_j, &
                                   iatom=iatom, jatom=jatom, katom=katom)
         dik = NORM2(rik)

         IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
             ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

         jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (jcell > nimg .OR. jcell < 1) CYCLE

         IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
             ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

         kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
         IF (kcell > nimg .OR. kcell < 1) CYCLE

         !TODO: at the moment, this is a descrete yes/no criterion to include an atom in the RI basis
         !      maybe we should put a smooth bump function like in Head-Gordon
         IF (dik > ri_data%kp_RI_range) CYCLE
         RI_cells(kcell) = 1

         IF (jcell == 1 .AND. iatom == jatom) nRI_per_atom(iatom) = nRI_per_atom(iatom) + ri_data%bsizes_RI(katom)
      END DO
      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)
      CALL para_env%sum(RI_cells)
      CALL para_env%sum(nRI_per_atom)

      ri_data%ncell_RI = 0
      ALLOCATE (ri_data%img_to_RI_cell(nimg))
      ri_data%img_to_RI_cell = 0
      DO i_img = 1, nimg
         IF (RI_cells(i_img) > 0) THEN
            ri_data%ncell_RI = ri_data%ncell_RI + 1
            ri_data%img_to_RI_cell(i_img) = ri_data%ncell_RI
         END IF
      END DO
      CALL neighbor_list_3c_destroy(nl_3c)

      !Print some info
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, FMT="(/T3,A,F31.3,A)") &
            "KP-HFX_RI_INFO| RI basis extension radius:", ri_data%kp_RI_range*angstrom, " Ang"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I16,A)") &
            "KP-HFX_RI_INFO| The extended RI bases cover up to ", ri_data%ncell_RI, " unit cells"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I18)") &
            "KP-HFX_RI_INFO| Average number of sgf in extended RI bases:", SUM(nRI_per_atom)/natom
         WRITE (ri_data%unit_nr, FMT="(T3,A,F13.3,A)") &
            "KP-HFX_RI_INFO| Consider all image cells within a radius of ", ri_data%kp_image_range*angstrom, " Ang"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I27/)") &
            "KP-HFX_RI_INFO| Number of image cells considered: ", ri_data%nimg
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_kp_and_ri_images

! **************************************************************************************************
!> \brief A routine that creates tensors structure for rho_ao and 3c_ints in a stacked format for
!>        the efficient contractions of rho_sigma^0,lambda^c * (mu^0 sigam^a | P)
!> \param rho_stack ...
!> \param ints_stack ...
!> \param rho_template ...
!> \param ints_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_stack_tensors(rho_stack, ints_stack, rho_template, ints_template, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: rho_stack, ints_stack
      TYPE(dbt_type), INTENT(INOUT)                      :: rho_template, ints_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i_img, nblks, nblks_3c(3), nimg, &
                                                            pdims_2d(2), pdims_3d(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_stack, dist1, &
                                                            dist2, dist3, dist_stack
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env)

      CALL get_qs_env(qs_env, para_env=para_env)

      nimg = ri_data%nimg
      nblks = SIZE(ri_data%bsizes_AO_split)
      ALLOCATE (bsizes_stack(nimg*nblks))
      DO i_img = 1, nimg
         bsizes_stack((i_img - 1)*nblks + 1:i_img*nblks) = ri_data%bsizes_AO_split(:)
      END DO

      ALLOCATE (dist1(nblks), dist2(nblks), dist_stack(nimg*nblks))
      CALL dbt_get_info(rho_template, proc_dist_1=dist1, proc_dist_2=dist2)
      DO i_img = 1, nimg
         dist_stack((i_img - 1)*nblks + 1:i_img*nblks) = dist1(:)
      END DO

      !First 2c tensor matches the distribution of template
      CALL dbt_distribution_new(t_dist, ri_data%pgrid_2d, dist_stack, dist2)
      CALL dbt_create(rho_stack(1), "RHO_stack", t_dist, [1], [2], bsizes_stack, ri_data%bsizes_AO_split)
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist_stack)

      !Second 2c tensor has optimal pgrid
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env, pdims_2d, pgrid, tensor_dims=[nimg*nblks, nblks])
      CALL create_2c_tensor(rho_stack(2), dist1, dist2, pgrid, bsizes_stack, &
                            ri_data%bsizes_AO_split, name="RHO_stack")
      DEALLOCATE (dist1, dist2)
      CALL dbt_pgrid_destroy(pgrid)

      CALL dbt_get_info(ints_template, nblks_total=nblks_3c)
      ALLOCATE (dist1(nblks_3c(1)), dist2(nblks_3c(2)), dist3(nblks_3c(3)))
      ALLOCATE (dist_stack(nimg*nblks_3c(3)), bsizes_RI_ext(nblks_3c(2)))
      CALL dbt_get_info(ints_template, proc_dist_1=dist1, proc_dist_2=dist2, &
                        proc_dist_3=dist3, blk_size_2=bsizes_RI_ext)
      DO i_img = 1, nimg
         dist_stack((i_img - 1)*nblks_3c(3) + 1:i_img*nblks_3c(3)) = dist3(:)
      END DO

      !First 3c tensor matches the distribution of template
      CALL dbt_distribution_new(t_dist, ri_data%pgrid_1, dist1, dist2, dist_stack)
      CALL dbt_create(ints_stack(1), "ints_stack", t_dist, [1, 2], [3], ri_data%bsizes_AO_split, &
                      bsizes_RI_ext, bsizes_stack)
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist3, dist_stack)

      !Second 3c tensor has optimal pgrid
      pdims_3d = 0
      CALL dbt_pgrid_create(para_env, pdims_3d, pgrid, tensor_dims=[nblks_3c(1), nblks_3c(2), nimg*nblks_3c(3)])
      CALL create_3c_tensor(ints_stack(2), dist1, dist2, dist3, pgrid, ri_data%bsizes_AO_split, &
                            bsizes_RI_ext, bsizes_stack, [1, 2], [3], name="ints_stack")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbt_pgrid_destroy(pgrid)

   END SUBROUTINE get_stack_tensors

! **************************************************************************************************
!> \brief Fill the stack 3c integrals (mu^0 sigma^a | P) tensors with a_img according to the order in ac_pairs
!> \param ints_stack ...
!> \param t_3c_ints ...
!> \param ac_pairs ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE fill_ints_stack(ints_stack, t_3c_ints, ac_pairs, ri_data)
      TYPE(dbt_type), INTENT(INOUT)                      :: ints_stack
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_ints
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: ac_pairs
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: a_img, i_img, ind(3), nblks, nimg
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbt_iterator_type)                            :: iter

      !We loop over the a images from the ac_pairs, then copy the 3c ints to the correct spot in
      !in the stack tensor (corresponding to pair index). Distributions match by construction
      nimg = ri_data%nimg
      nblks = SIZE(ri_data%bsizes_AO_split)

      DO i_img = 1, nimg
         IF (ANY(ac_pairs(:, i_img) == 0)) CYCLE
         a_img = ac_pairs(1, i_img)

!$OMP PARALLEL DEFAULT(NONE) SHARED(a_img,i_img,t_3c_ints,ints_stack,nblks) PRIVATE(iter,ind,blk,found)
         CALL dbt_iterator_start(iter, t_3c_ints(a_img))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)
            CALL dbt_get_block(t_3c_ints(a_img), ind, blk, found)
            IF (.NOT. found) CYCLE

            CALL dbt_put_block(ints_stack, [ind(1), ind(2), (i_img - 1)*nblks + ind(3)], SHAPE(blk), blk)
            DEALLOCATE (blk)
         END DO
         CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      END DO !i_img
      CALL dbt_finalize(ints_stack)

   END SUBROUTINE fill_ints_stack

! **************************************************************************************************
!> \brief Fill the stack rho_ao_sigma^0,lambda^c tensors with c_img according to the order in ac_pairs
!> \param rho_stack ...
!> \param rho_ao_t ...
!> \param ac_pairs ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE fill_rho_stack(rho_stack, rho_ao_t, ac_pairs, ri_data)
      TYPE(dbt_type), INTENT(INOUT)                      :: rho_stack
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: rho_ao_t
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: ac_pairs
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: c_img, i_img, ind(2), nblks, nimg
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(dbt_iterator_type)                            :: iter

      !We loop over the a images from the ac_pairs, then copy the 3c ints to the correct spot in
      !in the stack tensor (corresponding to pair index). Distributions match by construction
      nimg = ri_data%nimg
      nblks = SIZE(ri_data%bsizes_AO_split)

      DO i_img = 1, nimg
         IF (ANY(ac_pairs(:, i_img) == 0)) CYCLE
         c_img = ac_pairs(2, i_img)

!$OMP PARALLEL DEFAULT(NONE) SHARED(c_img,i_img,rho_ao_t,rho_stack,nblks) PRIVATE(iter,ind,blk,found)
         CALL dbt_iterator_start(iter, rho_ao_t(c_img))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)
            CALL dbt_get_block(rho_ao_t(c_img), ind, blk, found)
            IF (.NOT. found) CYCLE

            CALL dbt_put_block(rho_stack, [(i_img - 1)*nblks + ind(1), ind(2)], SHAPE(blk), blk)
            DEALLOCATE (blk)
         END DO
         CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      END DO !i_img
      CALL dbt_finalize(rho_stack)

   END SUBROUTINE fill_rho_stack
END MODULE
