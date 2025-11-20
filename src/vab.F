
! **************************************************************************************************
!> \brief Set occupation of molecular orbitals
!> \par History
!>      - module to calculate V_AB based on Farazdel, Dupuis, Clementi, Avram,JACS 1990
!> \author  Pavan Behara
! **************************************************************************************************

MODULE vab
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE hfx_types,                       ONLY: hfx_create
   USE input_constants,                 ONLY: xc_none
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sync
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt,&
                                              hertz,&
                                              kcalmol,&
                                              kjmol,&
                                              wavenumbers
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE vab_qs_ks_methods,               ONLY: qs_ks_update_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'vab'

   PUBLIC :: calculate_vab, &
             write_real_1dmatrix_to, write_real_2dmatrix_to

CONTAINS

! **************************************************************************************************
!> \brief
!>
!>
!> \param
!> \param
!> \date
!> \author
!> \version
! **************************************************************************************************
   SUBROUTINE calculate_vab(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_vab', routineP = moduleN//':'//routineN

      INTEGER :: after, before, err_info, group, i_rep, iw, lwork, n_rep, nao_a_alfa, nao_a_beta, &
         nao_b_alfa, nao_b_beta, natom, nblkcols_total, nblkrows_total, nmo_a_alfa, &
         nmo_a_beta, nmo_b_alfa, nmo_b_beta, nforce_eval, nfullcols_total, &
         nfullrows_total, nimgs, nspins, val_chk
      LOGICAL                                            :: debug_print, do_hfx, first_iter_exc, &
                                                            gapw_xc, haa_test, kclass_k2_check, &
                                                            omit_headers, use_ref_cell
      REAL(kind=dp) :: a_coeff, b_coeff, c_coeff, d_prod_alfa, d_prod_beta, delta, detU_alfa, &
         detU_beta, detVt_alfa, detVt_beta, DH_prod, DT_prod, E_lower, E_upper, haa, haa_k2_klass, &
         haa_klass, haa_mx, hab, hab_klass_k2, habs_ab, half_sp, hbb, hbb_k2_klass, hbb_klass, &
         hbb_mx, hxx, k2__ab, K2_with_asym, K2_with_sym, sab, sdel, V_ab, val_chk1
      REAL(kind=dp), DIMENSION(:), POINTER               :: diag_alfa, diag_beta, work
      REAL(kind=dp), DIMENSION(:, :), POINTER :: A_alfa, A_beta, A_hat_alfa, A_hat_beta, B_alfa, &
         B_beta, B_hat_alfa, B_hat_beta, D_alfa, D_beta, matrix_s_real, matrix_tmp1, matrix_tmp2, &
         P_alfa, P_beta, PA_hat_alfa, PA_hat_beta, PB_hat_alfa, PB_hat_beta, Pp_asym_alfa, &
         Pp_asym_beta, Pp_sym_alfa, Pp_sym_beta, T_alfa, T_beta, tmp1, tmp2, U_alfa, U_beta, &
         Vt_alfa, Vt_beta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp, fm_struct_tmp2
      TYPE(cp_fm_type), POINTER :: matrix_s_fm, mo_coeff_a_alfa, mo_coeff_a_beta, mo_coeff_b_alfa, &
         mo_coeff_b_beta, P_alfa_fm, P_alfa_fm_asym, P_alfa_fm_sym, P_beta_fm, P_beta_fm_asym, &
         P_beta_fm_sym, PA_hat_alfa_fm, PA_hat_beta_fm, PB_hat_alfa_fm, PB_hat_beta_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_modified_a, &
                                                            density_modified_asym, &
                                                            density_modified_b, &
                                                            density_modified_sym, matrix_s, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos_a, mos_b
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy_a, energy_b
      TYPE(qs_environment_type), POINTER                 :: qs_env_a, qs_env_b
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: force_env_section, hfx_presence, &
                                                            hfx_section, hfx_tmp, mixed_section, &
                                                            root_section

      debug_print = .FALSE.

      NULLIFY (qs_env_a, qs_env_b)
      CPASSERT(ASSOCIATED(force_env))
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section, &
                         root_section=root_section)
      mixed_section => section_vals_get_subs_vals(force_env_section, "MIXED")
!Taking the molecular orbitals output of two states A and B
      nforce_eval = SIZE(force_env%sub_force_env)
      WRITE (124, *) "NFORCE_EVAL is ", nforce_eval
      logger => cp_get_default_logger()
      qs_env_a => force_env%sub_force_env(1)%force_env%qs_env
      qs_env_b => force_env%sub_force_env(2)%force_env%qs_env
      NULLIFY (mos_a, mos_b, energy_a, energy_b, ks_env, blacs_env, dft_control, para_env)
      CALL get_qs_env(qs_env=qs_env_a, mos=mos_a, energy=energy_a)
      CALL get_qs_env(qs_env=qs_env_b, mos=mos_b, energy=energy_b, &
                      ks_env=ks_env, blacs_env=blacs_env, &
                      dft_control=dft_control, para_env=para_env)
      iw = cp_print_key_unit_nr(logger, qs_env_a%input, "DFT%SCF%PRINT%ITERATION_INFO", &
                                extension=".scfLog")

      CALL get_mo_set(mos_a(1)%mo_set, mo_coeff=mo_coeff_a_alfa, nmo=nmo_a_alfa, nao=nao_a_alfa)
      CALL get_mo_set(mos_a(2)%mo_set, mo_coeff=mo_coeff_a_beta, nmo=nmo_a_beta, nao=nao_a_beta)

      CALL get_mo_set(mos_b(1)%mo_set, mo_coeff=mo_coeff_b_alfa, nmo=nmo_b_alfa, nao=nao_b_alfa)
      CALL get_mo_set(mos_b(2)%mo_set, mo_coeff=mo_coeff_b_beta, nmo=nmo_b_beta, nao=nao_b_beta)

      IF((nao_a_alfa .NE. nao_a_beta) .OR. (nao_b_alfa .NE. nao_b_beta) .OR. (nao_a_alfa .NE. nao_b_alfa)) &
                   CALL cp_abort(__LOCATION__, "Mismatch in basis sets of initial and final state MOs")

      !overlap matrix = matrix_s
      CALL get_ks_env(ks_env=ks_env, matrix_s=matrix_s)
      IF (debug_print) THEN
         CALL cp_fm_write_formatted(mo_coeff_a_alfa, iw, "mo_coeff_a_alfa")
         CALL cp_fm_write_formatted(mo_coeff_a_beta, iw, "mo_coeff_a_beta")
         CALL cp_fm_write_formatted(mo_coeff_b_alfa, iw, "mo_coeff_b_alfa")
         CALL cp_fm_write_formatted(mo_coeff_b_beta, iw, "mo_coeff_b_beta")
      END IF
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nfullrows_total, &
                          nfullcols_total=nfullcols_total, nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nfullrows_total, ncol_global=nfullcols_total, &
                               nrow_block=nblkrows_total, ncol_block=nblkcols_total, &
                               context=blacs_env, para_env=para_env)
      CALL cp_fm_create(matrix_s_fm, matrix_struct=fm_struct_tmp)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, matrix_s_fm)

      ALLOCATE (A_hat_alfa(nao_a_alfa, nmo_a_alfa))
      ALLOCATE (B_hat_alfa(nao_b_alfa, nmo_b_alfa))
      ALLOCATE (A_beta(nao_a_beta, nmo_a_beta))
      ALLOCATE (B_beta(nao_b_beta, nmo_b_beta))
      ALLOCATE (A_alfa(nao_a_alfa, nmo_a_alfa))
      ALLOCATE (B_alfa(nao_b_alfa, nmo_b_alfa))
      ALLOCATE (matrix_s_real(nao_a_alfa, nao_b_alfa))

!!!!!!!!!!!!ALFA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (matrix_tmp1(nmo_b_alfa, nao_b_alfa))
      ALLOCATE (D_alfa(nmo_a_alfa, nmo_b_alfa))

      CALL cp_fm_get_submatrix(mo_coeff_a_alfa, A_alfa)
      CALL cp_fm_get_submatrix(mo_coeff_b_alfa, B_alfa)
      CALL cp_fm_get_submatrix(matrix_s_fm, matrix_s_real)

      CALL dgemm("T", "N", nmo_b_alfa, nao_a_alfa, nao_b_alfa, 1.0_dp, &
                 B_alfa, nao_b_alfa, matrix_s_real, nao_a_alfa, 0.0_dp, matrix_tmp1, nmo_b_alfa)
      CALL dgemm("N", "N", nmo_b_alfa, nmo_a_alfa, nao_b_alfa, 1.0_dp, &
                 matrix_tmp1, nmo_b_alfa, A_alfa, nao_a_alfa, 0.0_dp, D_alfa, nmo_a_alfa)

      ALLOCATE (U_alfa(SIZE(D_alfa, 1), SIZE(D_alfa, 1)))
      ALLOCATE (Vt_alfa(SIZE(D_alfa, 2), SIZE(D_alfa, 2)))
      ALLOCATE (diag_alfa(MIN(SIZE(D_alfa, 1), SIZE(D_alfa, 2))))

      lwork = MAX(1, 3*MIN(SIZE(D_alfa, 1), SIZE(D_alfa, 2)) + &
                  MAX(SIZE(D_alfa, 1), SIZE(D_alfa, 2)), &
                  5*MIN(SIZE(D_alfa, 1), SIZE(D_alfa, 2)))
      ALLOCATE (work(MAX(1, lwork)))

      !!!matrix, svalues, mrow, ncol, uvec, vtvec

      CALL dgesvd('A', 'A', SIZE(D_alfa, 1), SIZE(D_alfa, 2), D_alfa, &
                  SIZE(D_alfa, 1), diag_alfa, U_alfa, SIZE(D_alfa, 1), Vt_alfa, &
                  SIZE(D_alfa, 2), work, lwork, err_info)

      CALL dgemm('N', 'T', nao_a_alfa, nmo_a_alfa, nmo_a_alfa, 1.0_dp, &
                 A_alfa, nao_a_alfa, Vt_alfa, nmo_a_alfa, 0.0_dp, A_hat_alfa, nao_a_alfa)

      CALL dgemm('N', 'N', nao_b_alfa, nmo_b_alfa, nmo_b_alfa, 1.0_dp, &
                 B_alfa, nao_b_alfa, U_alfa, nmo_b_alfa, 0.0_dp, B_hat_alfa, nao_b_alfa)

      ALLOCATE (P_alfa(nao_a_alfa, nao_b_alfa))
      ALLOCATE (T_alfa(nmo_a_alfa, nmo_b_alfa))
      ALLOCATE (tmp1(nao_a_alfa, nmo_a_alfa))
      T_alfa = 0.0_dp
      CALL product_except_self(diag_alfa, T_alfa)
      IF (debug_print) THEN
         CALL write_real_1dmatrix_to(diag_alfa, 124, "diag_alfa")
         CALL write_real_2dmatrix_to(U_alfa, 124, "U_alfa")
         CALL write_real_2dmatrix_to(Vt_alfa, 124, "Vt_alfa")
         CALL write_real_2dmatrix_to(A_hat_alfa, 124, "A_hat_alfa")
         CALL write_real_2dmatrix_to(B_hat_alfa, 124, "B_hat_alfa")
         CALL write_real_2dmatrix_to(T_alfa, 124, "T_alfa")
      END IF
      CALL dgemm('N', 'N', nao_a_alfa, SIZE(T_alfa, 2), nmo_a_alfa, 1.0_dp, &
                 A_hat_alfa, nao_a_alfa, T_alfa, SIZE(T_alfa, 2), 0.0_dp, tmp1, nao_a_alfa)
      CALL dgemm('N', 'T', nao_a_alfa, nao_b_alfa, nmo_a_alfa, 1.0_dp, &
                 tmp1, nao_a_alfa, B_hat_alfa, nao_b_alfa, 0.0_dp, P_alfa, nao_a_alfa)

!!!!!!!!!!!!! BETA  !!!!!!!!!!!!!!
      ALLOCATE (A_hat_beta(nao_a_beta, nmo_a_beta))
      ALLOCATE (B_hat_beta(nao_b_beta, nmo_b_beta))
      ALLOCATE (matrix_tmp2(nmo_b_beta, nao_b_beta))
      ALLOCATE (D_beta(nmo_a_beta, nmo_b_beta))

      CALL cp_fm_get_submatrix(mo_coeff_a_beta, A_beta)
      CALL cp_fm_get_submatrix(mo_coeff_b_beta, B_beta)

      CALL dgemm("T", "N", nmo_b_beta, nao_a_beta, nao_b_beta, 1.0_dp, &
                 B_beta, nao_b_beta, matrix_s_real, nao_a_beta, 0.0_dp, matrix_tmp2, nmo_b_beta)
      CALL dgemm("N", "N", nmo_b_beta, nmo_a_beta, nao_b_beta, 1.0_dp, &
                 matrix_tmp2, nmo_b_beta, A_beta, nao_a_beta, 0.0_dp, D_beta, nmo_a_beta)

      ALLOCATE (U_beta(SIZE(D_beta, 1), SIZE(D_beta, 1)))
      ALLOCATE (Vt_beta(SIZE(D_beta, 2), SIZE(D_beta, 2)))
      ALLOCATE (diag_beta(MIN(SIZE(D_beta, 1), SIZE(D_beta, 2))))

      lwork = MAX(1, 3*MIN(SIZE(D_beta, 1), SIZE(D_beta, 2)) + &
                  MAX(SIZE(D_beta, 1), SIZE(D_beta, 2)), &
                  5*MIN(SIZE(D_beta, 1), SIZE(D_beta, 2)))
      DEALLOCATE (work)
      ALLOCATE (work(MAX(1, lwork)))

      !!!matrix, svalues, mrow, ncol, uvec, vtvec

      CALL dgesvd('A', 'A', SIZE(D_beta, 1), SIZE(D_beta, 2), D_beta, &
                  SIZE(D_beta, 1), diag_beta, U_beta, SIZE(D_beta, 1), Vt_beta, &
                  SIZE(D_beta, 2), work, lwork, err_info)

      CALL dgemm('N', 'T', nao_a_beta, nmo_a_beta, nmo_a_beta, 1.0_dp, &
                 A_beta, nao_a_beta, Vt_beta, nmo_a_beta, 0.0_dp, A_hat_beta, nao_a_beta)

      CALL dgemm('N', 'N', nao_b_beta, nmo_b_beta, nmo_b_beta, 1.0_dp, &
                 B_beta, nao_b_beta, U_beta, nmo_b_beta, 0.0_dp, B_hat_beta, nao_b_beta)

      ALLOCATE (P_beta(nao_a_beta, nao_b_beta))
      ALLOCATE (T_beta(nmo_a_beta, nmo_b_beta))
      ALLOCATE (tmp2(nao_a_beta, nmo_a_beta))

      CALL product_except_self(diag_beta, T_beta)
      IF (debug_print) THEN
         CALL write_real_1dmatrix_to(diag_beta, 124, "diag_beta")
         CALL write_real_2dmatrix_to(U_beta, 124, "U_beta")
         CALL write_real_2dmatrix_to(Vt_beta, 124, "Vt_beta")
         CALL write_real_2dmatrix_to(A_hat_beta, 124, "A_hat_beta")
         CALL write_real_2dmatrix_to(B_hat_beta, 124, "B_hat_beta")
         CALL write_real_2dmatrix_to(T_beta, 124, "T_beta")
      END IF
      CALL dgemm('N', 'N', nao_a_beta, SIZE(T_beta, 2), nmo_a_beta, 1.0_dp, &
                 A_hat_beta, nao_a_beta, T_beta, SIZE(T_beta, 2), 0.0_dp, tmp2, nao_a_beta)
      CALL dgemm('N', 'T', nao_a_beta, nao_b_beta, nmo_a_beta, 1.0_dp, &
                 tmp2, nao_a_beta, B_hat_beta, nao_b_beta, 0.0_dp, P_beta, nao_a_beta)

      d_prod_alfa = PRODUCT(diag_alfa)
      d_prod_beta = PRODUCT(diag_beta)
      DH_prod = d_prod_alfa*d_prod_beta
      CALL determ(U_alfa, detU_alfa)
      CALL determ(U_beta, detU_beta)
      CALL determ(Vt_alfa, detVt_alfa)
      CALL determ(Vt_beta, detVt_beta)
      WRITE (124, *) "determinants are detU_alfa, detU_beta, detVt_alfa, detVt_beta", &
         detU_alfa, detU_beta, detVt_alfa, detVt_beta
      DT_prod = detU_alfa*detU_beta*detVt_alfa*detVt_beta
      sab = DH_prod*DT_prod
      WRITE (124, *) "DH_PROD", DH_prod
      WRITE (124, *) "DT_PROD", DT_prod
      WRITE (124, *) "S_AB", sab
      P_alfa = d_prod_beta*P_alfa
      P_beta = d_prod_alfa*P_beta
      IF (debug_print) THEN
         CALL write_real_2dmatrix_to(P_alfa, 124, "P_alfa")
         CALL write_real_2dmatrix_to(P_beta, 124, "P_beta")
      END IF
!!!!!!!!!! MODIFIED DENSITY IS IN REAL(:,:) MATRICES, HAVE TO CONVERT THEM TO
!CP2K FORMATS
!CREATING A FULL MATRIX STRUCTURE AND THEN COPYING TO DBCSR FORMAT

      CALL dbcsr_get_info(mos_b(1)%mo_set%mo_coeff_b, nfullrows_total=nfullrows_total, &
                          nfullcols_total=nfullcols_total, nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total)
      CALL cp_fm_struct_create(fm_struct_tmp2, nrow_global=nfullrows_total, ncol_global=nfullcols_total, &
                               nrow_block=nblkrows_total, ncol_block=nblkcols_total, &
                               context=blacs_env, para_env=para_env)

      CALL cp_fm_create(P_alfa_fm, matrix_struct=fm_struct_tmp)
      CALL cp_fm_create(P_beta_fm, matrix_struct=fm_struct_tmp)

      CALL cp_fm_set_submatrix(P_alfa_fm, P_alfa)
      CALL cp_fm_set_submatrix(P_beta_fm, P_beta)

      IF (debug_print) THEN
         CALL cp_fm_write_formatted(P_alfa_fm, iw, "P_alfa_fullmatrix")
         CALL cp_fm_write_formatted(P_beta_fm, iw, "P_beta_fullmatrix")
      END IF
!#################  Checking E_A is same before and after transformation  #################
      NULLIFY (my_cell, my_cell_ref, atomic_kind_set, particle_set, &
               qs_kind_set, subsys)
      kclass_k2_check = .FALSE.
      hxx = 0.0_dp
      first_iter_exc = .FALSE.
      hfx_presence => section_vals_get_subs_vals(qs_env_b%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_presence, explicit=do_hfx)
      CALL get_ks_env(ks_env, subsys=subsys)
      CALL qs_subsys_get(subsys, &
                         cell=my_cell, &
                         cell_ref=my_cell_ref, &
                         use_ref_cell=use_ref_cell, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)

      IF (.NOT. do_hfx) THEN
         first_iter_exc = .TRUE.
         CALL section_vals_val_set(mixed_section, "VAB%HF%FRACTION", &
                                   r_val=1.0_dp)
         CALL section_vals_val_set(qs_env_b%input, "DFT%XC%XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                   i_val=xc_none)
         CALL section_vals_val_set(qs_env_b%input, "DFT%XC%HF%FRACTION", &
                                   r_val=1.0_dp)
         hfx_tmp => section_vals_get_subs_vals(mixed_section, "VAB%HF")
         CALL section_vals_get(hfx_tmp, n_repetition=n_rep)
         CALL section_vals_val_get(mixed_section, "VAB%HF%INTERACTION_POTENTIAL%POTENTIAL_TYPE", i_val=val_chk)
         !WRITE(79,*) "HF potential type in vab input", val_chk
         !WRITE(79,*) "repetition under vab", n_rep
         hfx_section => section_vals_get_subs_vals(qs_env_b%input, "DFT%XC%HF")
         CALL section_vals_get(hfx_section, n_repetition=n_rep)
         !WRITE(79,*) "repetition under DFT%XC", n_rep
         DO i_rep = 1, n_rep
            NULLIFY (hfx_section)
            CALL section_vals_duplicate(hfx_tmp, hfx_section, i_rep_start=i_rep, i_rep_end=i_rep)
         END DO
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=val_chk1)
         !WRITE(79,*) "HF fraction", val_chk1
         CALL section_vals_val_get(hfx_section, "INTERACTION_POTENTIAL%POTENTIAL_TYPE", i_val=val_chk)
         !WRITE(79,*) "HF potential type in xc hf section", val_chk
         !FLUSH(79)
      ELSE
         CALL section_vals_val_set(qs_env_b%input, "DFT%XC%XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                   i_val=xc_none)
         CALL section_vals_val_set(qs_env_b%input, "DFT%XC%HF%FRACTION", &
                                   r_val=1.0_dp)
         hfx_section => section_vals_get_subs_vals(qs_env_b%input, "DFT%XC%HF")
      END IF
      natom = SIZE(particle_set)
      IF (.NOT. do_hfx) THEN
         CALL hfx_create(qs_env_b%x_data, para_env, hfx_section, natom, atomic_kind_set, &
                         qs_kind_set, dft_control, my_cell)
      END IF
      NULLIFY (density_modified_a)
      nspins = dft_control%nspins
      nimgs = dft_control%nimages
      gapw_xc = dft_control%qs_control%gapw_xc

      IF (iw > 0) WRITE (iw, *) "nspins, nimgs, gapw_xc", nspins, nimgs, gapw_xc, NEW_LINE('A')
      IF (iw > 0) WRITE (iw, *) "####### CHECKS ON A and B #######"
      IF (iw > 0) WRITE (iw, *) " E_A before:", energy_a%total
      IF (iw > 0) WRITE (iw, *) " E_B before:", energy_b%total
      IF (iw > 0) WRITE (iw, *) "####### - #######"
      haa = 10000_dp
      hbb = 10000_dp
      ALLOCATE (PA_hat_alfa(nao_a_alfa, nao_a_alfa))
      ALLOCATE (PA_hat_beta(nao_a_beta, nao_a_beta))
      CALL dgemm('N', 'T', nao_a_alfa, nao_a_alfa, nmo_a_alfa, 1.0_dp, &
                 A_hat_alfa, nao_a_alfa, A_hat_alfa, nao_a_alfa, 0.0_dp, PA_hat_alfa, nao_a_alfa)
      CALL dgemm('N', 'T', nao_a_beta, nao_a_beta, nmo_a_beta, 1.0_dp, &
                 A_hat_beta, nao_a_beta, A_hat_beta, nao_a_beta, 0.0_dp, PA_hat_beta, nao_a_beta)
      WRITE (124, *) "############ CHEKCING E_A ############"
      IF (debug_print) THEN
         CALL write_real_2dmatrix_to(PA_hat_alfa, 124, "PA_hat_alfa")
         CALL write_real_2dmatrix_to(PA_hat_beta, 124, "PA_hat_beta")
      END IF

      CALL cp_fm_create(PA_hat_alfa_fm, matrix_struct=fm_struct_tmp)
      CALL cp_fm_create(PA_hat_beta_fm, matrix_struct=fm_struct_tmp)

      CALL cp_fm_set_submatrix(PA_hat_alfa_fm, PA_hat_alfa)
      CALL cp_fm_set_submatrix(PA_hat_beta_fm, PA_hat_beta)

      IF (debug_print) THEN
         CALL cp_fm_write_formatted(PA_hat_alfa_fm, iw, "PA_hat_alfa_fullmatrix")
         CALL cp_fm_write_formatted(PA_hat_beta_fm, iw, "PA_hat_beta_fullmatrix")
      END IF

      CALL get_qs_env(qs_env_b, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)

      ALLOCATE (density_modified_a(2))
      CALL dbcsr_init_p(density_modified_a(1)%matrix)
      CALL dbcsr_create(density_modified_a(1)%matrix, &
                        template=rho_ao(1)%matrix, &
                        name="PA_hat_alfa", &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_init_p(density_modified_a(2)%matrix)
      CALL dbcsr_create(density_modified_a(2)%matrix, &
                        template=rho_ao(2)%matrix, &
                        name="PA_hat_beta", &
                        matrix_type=dbcsr_type_symmetric)

      CALL copy_fm_to_dbcsr(PA_hat_alfa_fm, density_modified_a(1)%matrix)
      CALL copy_fm_to_dbcsr(PA_hat_beta_fm, density_modified_a(2)%matrix)
      rho_ao_kp(1:nspins, 1:1) => density_modified_a(1:nspins)
      rho_ao => density_modified_a
      CALL qs_rho_set(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)
      CALL qs_rho_update_rho(rho, qs_env=qs_env_b)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE., potential_changed=.TRUE.)
      CALL qs_ks_update_qs_env(qs_env_b, just_energy=.TRUE., &
                               calculate_forces=.FALSE., hxx=hxx, &
                               split_energy=.TRUE., first_iter_exc=first_iter_exc)
      haa = hxx
      haa_mx = hxx - energy_b%ex
      IF (kclass_k2_check) THEN
         CALL K2AB(PA_hat_alfa, PA_hat_beta, haa_k2_klass)
      ELSE
         haa_k2_klass = 0.0_dp
      END IF
      haa_klass = haa_mx + 0.5_dp*haa_k2_klass
      IF (iw > 0) WRITE (iw, *) "CHECK ON A IS DONE, E_A after passing modified PA_hat density:", energy_b%total
      IF (kclass_k2_check) THEN
         ALLOCATE (Pp_sym_alfa(nao_a_alfa, nao_b_alfa))
         ALLOCATE (Pp_asym_alfa(nao_a_alfa, nao_b_alfa))
         ALLOCATE (Pp_sym_beta(nao_a_alfa, nao_b_alfa))
         ALLOCATE (Pp_asym_beta(nao_a_alfa, nao_b_alfa))
         Pp_sym_alfa = 0.5_dp*(P_alfa + TRANSPOSE(P_alfa))
         Pp_asym_alfa = 0.5_dp*(P_alfa - TRANSPOSE(P_alfa))
         Pp_sym_beta = 0.5_dp*(P_beta + TRANSPOSE(P_beta))
         Pp_asym_beta = 0.5_dp*(P_beta - TRANSPOSE(P_beta))
         CALL K2AB(Pp_sym_alfa, Pp_sym_beta, K2_with_sym)
         CALL K2AB(Pp_asym_alfa, Pp_asym_beta, K2_with_asym)
         WRITE (45, *) "%%%% K2_AB with KLASS %%%%"
         WRITE (45, *) "K2_with_sym ", K2_with_sym
         WRITE (45, *) "K2_with_asym ", K2_with_asym
         WRITE (45, *) "Total K2 ", K2_with_sym + K2_with_asym
         WRITE (45, *) "---"
         WRITE (45, *) "Multiplying with prefactor 0.5*dtprod/dhprod "
         WRITE (45, *) "---"
         WRITE (45, *) "K2_with_sym ", 0.5_dp*DT_prod*K2_with_sym/DH_prod
         WRITE (45, *) "K2_with_asym ", 0.5_dp*DT_prod*K2_with_asym/DH_prod
         WRITE (45, *) "Total K2 ", 0.5_dp*DT_prod*(K2_with_sym + K2_with_asym)/DH_prod
         WRITE (45, *) "---- End of K2 check ----"
         DEALLOCATE (Pp_sym_alfa, Pp_asym_alfa, Pp_sym_beta, Pp_asym_beta)
      END IF

!#################  Checking E_B is same before and after transformation  #################
      NULLIFY (density_modified_b, rho_ao, rho_ao_kp)
      ALLOCATE (PB_hat_alfa(nao_b_alfa, nao_b_alfa))
      ALLOCATE (PB_hat_beta(nao_b_beta, nao_b_beta))
      CALL dgemm('N', 'T', nao_b_alfa, nao_b_alfa, nmo_b_alfa, 1.0_dp, &
                 B_hat_alfa, nao_b_alfa, B_hat_alfa, nao_b_alfa, 0.0_dp, PB_hat_alfa, nao_b_alfa)
      CALL dgemm('N', 'T', nao_b_beta, nao_b_beta, nmo_b_beta, 1.0_dp, &
                 B_hat_beta, nao_b_beta, B_hat_beta, nao_b_beta, 0.0_dp, PB_hat_beta, nao_b_beta)
      WRITE (124, *) "############ CHEKCING E_B ############"
      IF (debug_print) THEN
         CALL write_real_2dmatrix_to(PB_hat_alfa, 124, "PB_hat_alfa")
         CALL write_real_2dmatrix_to(PB_hat_beta, 124, "PB_hat_beta")
      END IF
      CALL cp_fm_create(PB_hat_alfa_fm, matrix_struct=fm_struct_tmp)
      CALL cp_fm_create(PB_hat_beta_fm, matrix_struct=fm_struct_tmp)

      CALL cp_fm_set_submatrix(PB_hat_alfa_fm, PB_hat_alfa)
      CALL cp_fm_set_submatrix(PB_hat_beta_fm, PB_hat_beta)

      IF (debug_print) THEN
         CALL cp_fm_write_formatted(PB_hat_alfa_fm, iw, "PB_hat_alfa_fullmatrix")
         CALL cp_fm_write_formatted(PB_hat_beta_fm, iw, "PB_hat_beta_fullmatrix")
      END IF
      CALL get_qs_env(qs_env_b, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)

      ALLOCATE (density_modified_b(2))
      CALL dbcsr_init_p(density_modified_b(1)%matrix)
      CALL dbcsr_create(density_modified_b(1)%matrix, &
                        template=rho_ao(1)%matrix, &
                        name="PB_hat_alfa", &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_init_p(density_modified_b(2)%matrix)
      CALL dbcsr_create(density_modified_b(2)%matrix, &
                        template=rho_ao(2)%matrix, &
                        name="PB_hat_beta", &
                        matrix_type=dbcsr_type_symmetric)
      CALL copy_fm_to_dbcsr(PB_hat_alfa_fm, density_modified_b(1)%matrix)
      CALL copy_fm_to_dbcsr(PB_hat_beta_fm, density_modified_b(2)%matrix)
      rho_ao_kp(1:nspins, 1:1) => density_modified_b(1:nspins)
      rho_ao => density_modified_b
      CALL qs_rho_set(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)
      CALL qs_rho_update_rho(rho, qs_env=qs_env_b)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE., potential_changed=.TRUE.)
      CALL qs_ks_update_qs_env(qs_env_b, just_energy=.TRUE., &
                               calculate_forces=.FALSE., hxx=hxx, split_energy=.TRUE.)
      hbb = hxx
      hbb_mx = hxx - energy_b%ex
      IF (kclass_k2_check) THEN
         CALL K2AB(PB_hat_alfa, PB_hat_beta, hbb_k2_klass)
      ELSE
         hbb_k2_klass = 0.0_dp
      END IF
      hbb_klass = hbb_mx + 0.5_dp*hbb_k2_klass

      IF (iw > 0) WRITE (iw, *) "CHECK ON B IS DONE, E_B after passing modified PB_hat density:", energy_b%total
!#################
!#################  Doing anti-symmetric part for k2__ab
!#################
      haa_test = .FALSE.
      NULLIFY (density_modified_asym, rho_ao, rho_ao_kp)

      CALL get_qs_env(qs_env_b, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)

      ALLOCATE (density_modified_asym(2))
      CALL dbcsr_init_p(density_modified_asym(1)%matrix)
      CALL dbcsr_create(density_modified_asym(1)%matrix, &
                        template=rho_ao(1)%matrix, &
                        name="P_AB anti-symmetric alfa", &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_init_p(density_modified_asym(2)%matrix)
      CALL dbcsr_create(density_modified_asym(2)%matrix, &
                        template=rho_ao(2)%matrix, &
                        name="P_AB anti-symmetric beta", &
                        matrix_type=dbcsr_type_antisymmetric)

      CALL cp_fm_create(P_alfa_fm_asym, matrix_struct=fm_struct_tmp)
      CALL cp_fm_create(P_beta_fm_asym, matrix_struct=fm_struct_tmp)
      IF (haa_test) THEN
         CALL cp_fm_transpose(PB_hat_alfa_fm, P_alfa_fm_asym)
         CALL cp_fm_transpose(PB_hat_beta_fm, P_beta_fm_asym)
         CALL cp_fm_scale_and_add(-0.5_dp, P_alfa_fm_asym, 0.5_dp, PB_hat_alfa_fm)
         CALL cp_fm_scale_and_add(-0.5_dp, P_beta_fm_asym, 0.5_dp, PB_hat_beta_fm)
      ELSE
         CALL cp_fm_transpose(P_alfa_fm, P_alfa_fm_asym)
         CALL cp_fm_transpose(P_beta_fm, P_beta_fm_asym)
         CALL cp_fm_scale_and_add(-0.5_dp, P_alfa_fm_asym, 0.5_dp, P_alfa_fm)
         CALL cp_fm_scale_and_add(-0.5_dp, P_beta_fm_asym, 0.5_dp, P_beta_fm)
      END IF
      CALL copy_fm_to_dbcsr(P_alfa_fm_asym, density_modified_asym(1)%matrix)
      CALL copy_fm_to_dbcsr(P_beta_fm_asym, density_modified_asym(2)%matrix)

      rho_ao_kp(1:nspins, 1:1) => density_modified_asym(1:nspins)
      rho_ao => density_modified_asym
      CALL qs_rho_set(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)
      CALL section_vals_val_get(qs_env_b%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      CALL section_vals_val_get(qs_env_b%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
      IF (debug_print) THEN
         after = MIN(MAX(after, 1), 16)
         before = 2
         CALL cp_dbcsr_write_sparse_matrix(density_modified_asym(1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
         CALL cp_dbcsr_write_sparse_matrix(rho_ao(1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
         CALL cp_dbcsr_write_sparse_matrix(rho_ao_kp(1, 1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)

         CALL cp_dbcsr_write_sparse_matrix(density_modified_asym(2)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
      END IF
      CALL qs_rho_update_rho(rho, qs_env=qs_env_b)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE., potential_changed=.TRUE.)
      IF (iw > 0) WRITE (iw, *) "DOING ANTI-SYMMETRIC PART OF VAB"
      CALL qs_ks_update_qs_env(qs_env_b, just_energy=.TRUE., &
                               calculate_forces=.FALSE., split_energy=.TRUE.)

      k2__ab = energy_b%ex
      WRITE (120, *) "k2__ab from anti-symmetric part is", k2__ab

!#################
!#################  Doing symmetric part
!#################

      NULLIFY (density_modified_sym, rho_ao, rho_ao_kp)

      CALL get_qs_env(qs_env_b, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)

      ALLOCATE (density_modified_sym(2))
      CALL dbcsr_init_p(density_modified_sym(1)%matrix)
      CALL dbcsr_create(density_modified_sym(1)%matrix, &
                        template=rho_ao(1)%matrix, &
                        name="P_AB symmetric alfa", &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_init_p(density_modified_sym(2)%matrix)
      CALL dbcsr_create(density_modified_sym(2)%matrix, &
                        template=rho_ao(2)%matrix, &
                        name="P_AB symmetric beta", &
                        matrix_type=dbcsr_type_symmetric)

      CALL cp_fm_create(P_alfa_fm_sym, matrix_struct=fm_struct_tmp)
      CALL cp_fm_create(P_beta_fm_sym, matrix_struct=fm_struct_tmp)
      IF (haa_test) THEN
         CALL cp_fm_transpose(PB_hat_alfa_fm, P_alfa_fm_sym)
         CALL cp_fm_transpose(PB_hat_beta_fm, P_beta_fm_sym)
         CALL cp_fm_scale_and_add(0.5_dp, P_alfa_fm_sym, 0.5_dp, PB_hat_alfa_fm)
         CALL cp_fm_scale_and_add(0.5_dp, P_beta_fm_sym, 0.5_dp, PB_hat_beta_fm)
      ELSE
         CALL cp_fm_transpose(P_alfa_fm, P_alfa_fm_sym)
         CALL cp_fm_transpose(P_beta_fm, P_beta_fm_sym)
         CALL cp_fm_scale_and_add(0.5_dp, P_alfa_fm_sym, 0.5_dp, P_alfa_fm)
         CALL cp_fm_scale_and_add(0.5_dp, P_beta_fm_sym, 0.5_dp, P_beta_fm)
      END IF
      CALL copy_fm_to_dbcsr(P_alfa_fm_sym, density_modified_sym(1)%matrix)
      CALL copy_fm_to_dbcsr(P_beta_fm_sym, density_modified_sym(2)%matrix)

      rho_ao_kp(1:nspins, 1:1) => density_modified_sym(1:nspins)
      rho_ao => density_modified_sym
      CALL qs_rho_set(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp)
      CALL section_vals_val_get(qs_env_b%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      CALL section_vals_val_get(qs_env_b%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
      IF (debug_print) THEN
         after = MIN(MAX(after, 1), 16)
         before = 2
         CALL cp_dbcsr_write_sparse_matrix(density_modified_sym(1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
         CALL cp_dbcsr_write_sparse_matrix(rho_ao(1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
         CALL cp_dbcsr_write_sparse_matrix(rho_ao_kp(1, 1)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)

         CALL cp_dbcsr_write_sparse_matrix(density_modified_sym(2)%matrix, before, after, qs_env_b, para_env, &
                                           output_unit=iw, omit_headers=omit_headers)
      END IF
      CALL qs_rho_update_rho(rho, qs_env=qs_env_b)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE., potential_changed=.TRUE.)
      IF (iw > 0) WRITE (iw, *) "DOING SYMMETRIC PART OF VAB"
      CALL qs_ks_update_qs_env(qs_env_b, just_energy=.TRUE., &
                               calculate_forces=.FALSE., calc_vab=.TRUE., dtprod=DT_prod, &
                               dhprod=DH_prod, k2ab=k2__ab, hxx=hxx, haa=haa, hbb=hbb, split_energy=.TRUE.)

!#################
!#################  Done with AB part
!#################

!#################  VAB with KLASS #################
      IF (kclass_k2_check) THEN
         CALL K2AB(P_alfa, P_beta, hab_klass_k2)
         hab = hxx + 0.5_dp*DT_prod*hab_klass_k2/DH_prod
         haa = haa_klass
         hbb = hbb_klass

         V_ab = (1.0_dp/(1.0_dp - sab*sab))*SQRT(0.25_dp*(haa - hbb)*(haa - hbb) - (haa + hbb)*hab*sab + &
                                                 haa*hbb*sab*sab + hab*hab)
         habs_ab = hab*sab
         a_coeff = 1.0_dp - sab*sab
         b_coeff = -(haa + hbb - habs_ab - habs_ab)
         c_coeff = haa*hbb - hab*hab
         delta = b_coeff*b_coeff - (4.0_dp)*a_coeff*c_coeff
         IF (delta .LT. 0.0) STOP 'delta CANNOT BE negative, check b^2-4ac in vab.F Line 733'
         sdel = SQRT(delta)
         E_upper = (-b_coeff + sdel)/(2.0_dp*a_coeff)
         E_lower = (-b_coeff - sdel)/(2.0_dp*a_coeff)
         half_sp = sdel/(2.0_dp*a_coeff)

         IF (iw > 0) WRITE (iw, *) "            "
         IF (iw > 0) WRITE (iw, *) "            "
         IF (iw > 0) WRITE (iw, *) "########################"
         IF (iw > 0) WRITE (iw, *) "#########  ET From KLASS calculated K2 #########"
         IF (iw > 0) WRITE (iw, *) "------------------------"
         IF (iw > 0) WRITE (iw, *) "HAA KLASS <a|H|a>", haa
         IF (iw > 0) WRITE (iw, *) "HBB KLASS <b|H|b>", hbb
         IF (iw > 0) WRITE (iw, *) "K2AA_KLASS ", haa_k2_klass
         IF (iw > 0) WRITE (iw, *) "K2BB_KLASS ", hbb_k2_klass
         IF (iw > 0) WRITE (iw, *) "            "
         IF (iw > 0) WRITE (iw, *) "KLASS <a|K2|b>", 0.5_dp*DT_prod*hab_klass_k2/DH_prod
         IF (iw > 0) WRITE (iw, *) "HAB KLASS <a|H|b>", hab
         IF (iw > 0) WRITE (iw, *) "SAB", sab
         IF (iw > 0) WRITE (iw, *) "VAB KLASS", V_ab
         IF (iw > 0) WRITE (iw, *) " VAB in different units "
         IF (iw > 0) WRITE (iw, *) "------------------------"
         IF (iw > 0) WRITE (iw, *) "      ", V_ab*evolt, "eV"
         IF (iw > 0) WRITE (iw, *) "      ", V_ab*wavenumbers, "1/cm"
         IF (iw > 0) WRITE (iw, *) "      ", V_ab*kcalmol, "kcal/mol"
         IF (iw > 0) WRITE (iw, *) "      ", V_ab*kjmol, "kJ/mol"
         IF (iw > 0) WRITE (iw, *) "      ", V_ab*hertz, "Hertz"
         IF (iw > 0) WRITE (iw, *) "------------------------"
         IF (iw > 0) WRITE (iw, *) "------------------------"
         IF (iw > 0) WRITE (iw, *) "E+", E_upper
         IF (iw > 0) WRITE (iw, *) "E-", E_lower
         IF (iw > 0) WRITE (iw, *) "(E+ - E-)/2", half_sp
         IF (iw > 0) WRITE (iw, *) "      ", half_sp*evolt, "eV"
         IF (iw > 0) WRITE (iw, *) "      ", half_sp*wavenumbers, "1/cm"
         IF (iw > 0) WRITE (iw, *) "      ", half_sp*kcalmol, "kcal/mol"
         IF (iw > 0) WRITE (iw, *) "      ", half_sp*kjmol, "kJ/mol"
         IF (iw > 0) WRITE (iw, *) "      ", half_sp*hertz, "Hertz"
         IF (iw > 0) WRITE (iw, *) "------------------------"
         IF (iw > 0) WRITE (iw, *) "------------------------"
      END IF

      group = para_env%group
      CALL mp_sync(group)
!#################
!#################
!#################

      !DEALLOCATE(matrix_tmp1, matrix_tmp2, diag_alfa, diag_beta, T_alfa, T_beta, D_alfa, U_alfa, Vt_alfa, &
      !           D_beta, U_beta, Vt_beta, A_alfa, A_beta, B_alfa, B_beta, A_hat_alfa, B_hat_alfa, A_hat_beta, &
      !           B_hat_beta, P_alfa, P_beta, PA_hat_alfa, PA_hat_beta, PB_hat_alfa, PB_hat_beta, tmp1, tmp2, &
      !           matrix_s_real)
      !CALL cp_fm_release(mo_coeff_a_alfa)
      !CALL cp_fm_release(mo_coeff_a_beta)
      !CALL cp_fm_release(mo_coeff_b_alfa)
      !CALL cp_fm_release(mo_coeff_b_beta)
      !CALL cp_fm_release(matrix_s_fm)
      !CALL cp_fm_release(P_alfa_fm)
      !CALL cp_fm_release(P_beta_fm)
      !CALL cp_fm_release(PA_hat_alfa_fm)
      !CALL cp_fm_release(PA_hat_beta_fm)
      !CALL cp_fm_release(PB_hat_alfa_fm)
      !CALL cp_fm_release(PB_hat_beta_fm)
      !CALL cp_fm_release(P_alfa_fm_sym)
      !CALL cp_fm_release(P_beta_fm_sym)
      !CALL cp_fm_release(P_alfa_fm_asym)
      !CALL cp_fm_release(P_beta_fm_asym)
      !CALL cp_fm_struct_release(fm_struct_tmp)
      !CALL cp_fm_struct_release(fm_struct_tmp2)
      !CALL dbcsr_release(density_modified_a(1)%matrix)
      !CALL dbcsr_release(density_modified_a(2)%matrix)
      !DEALLOCATE(density_modified_a)
      !CALL dbcsr_release(density_modified_b(1)%matrix)
      !CALL dbcsr_release(density_modified_b(2)%matrix)
      !DEALLOCATE(density_modified_b)
      !CALL dbcsr_release(density_modified_sym(1)%matrix)
      !CALL dbcsr_release(density_modified_sym(2)%matrix)
      !DEALLOCATE(density_modified_sym)
      !CALL dbcsr_release(density_modified_asym(1)%matrix)
      !CALL dbcsr_release(density_modified_asym(2)%matrix)
      !DEALLOCATE(density_modified_asym)
   END SUBROUTINE calculate_vab

   SUBROUTINE klass(i, j, k, l, goto_case)
      INTEGER                                            :: i, j, k, l
      INTEGER, INTENT(OUT)                               :: goto_case

      CHARACTER(len=*), PARAMETER :: routineN = 'klass', routineP = moduleN//':'//routineN

      INTEGER                                            :: kount
      LOGICAL                                            :: IEQJ, IEQK, IEQL, JEQK, JEQL, KEQL

      kount = 0
      IEQJ = i .EQ. j
      KEQL = k .EQ. l
      IEQK = i .EQ. k
      JEQL = j .EQ. l
      IEQL = i .EQ. l
      JEQK = j .EQ. k

      IF (IEQJ .AND. KEQL) THEN
         goto_case = 5
         IF (IEQK) goto_case = 1
         RETURN
      END IF
      IF (IEQJ) kount = kount + 1
      IF (KEQL) kount = kount + 1
      IF (IEQK) kount = kount + 1
      IF (JEQL) kount = kount + 1
      IF (IEQL) kount = kount + 1
      IF (JEQK) kount = kount + 1
      IF (kount .EQ. 3) THEN
         goto_case = 2
         IF (IEQJ) goto_case = 3
         RETURN
      END IF
      IF (kount .EQ. 2) THEN
         goto_case = 4
         RETURN
      END IF
      IF (kount .GT. 3 .OR. kount .LT. 0) THEN
         goto_case = 10
         RETURN
      END IF
      IF (IEQJ) THEN
         goto_case = 6
         RETURN
      END IF
      IF (KEQL) THEN
         goto_case = 7
         RETURN
      END IF
      goto_case = 8
      RETURN
   END SUBROUTINE klass

   SUBROUTINE K2AB(P_alfa, P_beta, K2_AB)
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: P_alfa, P_beta
      REAL(kind=dp), INTENT(OUT)                         :: K2_AB

      CHARACTER(len=*), PARAMETER :: routineN = 'K2AB', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, klss, l, nints, stat
      REAL(KIND=dp)                                      :: val, val2

      K2_AB = 0.0_dp
      nints = 0
      OPEN (15, file="2elec_int.txt", status='OLD', action='READ', iostat=stat)
      IF (stat /= 0) STOP 'ERROR: FILE 2elec_int.txt NOT FOUND'
      READ (15, *)
      DO
         !READ(15, FORMAT="I12, I12, I12, I12, E16.8", iostat=stat) i, j, k, l, val
         READ (15, *, iostat=stat) i, j, k, l, val
         !STOP 'in K2_AB'
         IF (stat /= 0) EXIT
         nints = nints + 1
         CALL klass(i, j, k, l, klss)
         WRITE (191, *) "i,j,k,l,val", i, j, k, l, val, klss
         val2 = val + val
         IF (klss .EQ. 1) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j))*val
         END IF
         IF (klss .EQ. 2) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(j, l)*P_alfa(k, i) - P_beta(j, l)*P_beta(k, i))*val2
         END IF
         IF (klss .EQ. 3) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(i, k)*P_alfa(l, j) &
                             - P_beta(i, k)*P_beta(l, j))*val2
         END IF
         IF (klss .EQ. 4) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(j, l)*P_alfa(k, i) - P_beta(j, l)*P_beta(k, i) &
                             - P_alfa(i, k)*P_alfa(l, j) - P_beta(i, k)*P_beta(l, j) &
                             - P_alfa(j, k)*P_alfa(l, i) - P_beta(j, k)*P_beta(l, i))*val
         END IF
         IF (klss .EQ. 5) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) &
                             - P_beta(i, l)*P_beta(k, j))*val2
         END IF
         IF (klss .EQ. 6) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(i, k)*P_alfa(l, j) - P_beta(i, k)*P_beta(l, j))*val2
         END IF
         IF (klss .EQ. 7) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(j, l)*P_alfa(k, i) - P_beta(j, l)*P_beta(k, i))*val2
         END IF
         IF (klss .EQ. 8) THEN
            K2_AB = K2_AB + (-P_alfa(i, l)*P_alfa(k, j) - P_beta(i, l)*P_beta(k, j) &
                             - P_alfa(j, l)*P_alfa(k, i) - P_beta(j, l)*P_beta(k, i) &
                             - P_alfa(i, k)*P_alfa(l, j) - P_beta(i, k)*P_beta(l, j) &
                             - P_alfa(j, k)*P_alfa(l, i) - P_beta(j, k)*P_beta(l, i))*val2
         END IF
         IF (klss .EQ. 10) THEN
            STOP 'KLSS IS WRONG!!!!'
         END IF
      END DO
      CLOSE (15)
      WRITE (45, *) "Number of integrals read:", nints
   END SUBROUTINE K2AB

   SUBROUTINE determ(A, det)
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: A
      REAL(kind=dp), INTENT(OUT)                         :: det

      CHARACTER(len=*), PARAMETER :: routineN = 'determ', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, INFO, j, N
      INTEGER, DIMENSION(:), POINTER                     :: IPIV
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: TMP

      N = SIZE(A, 1)
      ALLOCATE (TMP(N, N))
      ALLOCATE (IPIV(N))
      DO i = 1, N
         DO j = 1, N
            TMP(i, j) = A(i, j)
         END DO
      END DO
      CALL dgetrf(N, N, TMP, N, IPIV, INFO)
      IF (INFO > 0) STOP
      det = 1.0_dp
      DO i = 1, N
         IF (ipiv(i) .NE. i) THEN
            det = -det*TMP(i, i)
         ELSE
            det = det*TMP(i, i)
         ENDIF
      END DO
      DEALLOCATE (TMP, IPIV)
   END SUBROUTINE determ

   SUBROUTINE product_except_self(diag_mat, T_mat)
      REAL(kind=dp), DIMENSION(:), POINTER               :: diag_mat
      REAL(kind=dp), DIMENSION(:, :), INTENT(OUT), &
         POINTER                                         :: T_mat

      CHARACTER(len=*), PARAMETER :: routineN = 'product_except_self', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, N
      REAL(kind=dp)                                      :: r_prod

      T_mat = 0.0_dp
      N = SIZE(diag_mat)
      IF (N == 0) RETURN
      DO i = 1, N
         DO j = 1, N
            IF (i == j) THEN
               T_mat(i, j) = 1.0_dp
            END IF
         END DO
      END DO
      DO i = 2, N
         T_mat(i, i) = T_mat(i - 1, i - 1)*diag_mat(i - 1)
      END DO
      r_prod = 1.0_dp
      DO i = N, 1, -1
         T_mat(i, i) = T_mat(i, i)*r_prod
         r_prod = r_prod*diag_mat(i)
      END DO
   END SUBROUTINE product_except_self

   SUBROUTINE write_real_2dmatrix_to(matrix, write_nr, matrix_name)
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: matrix
      INTEGER                                            :: write_nr
      CHARACTER(len=*)                                   :: matrix_name

      CHARACTER(len=*), PARAMETER :: routineN = 'write_real_2dmatrix_to', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, M

      WRITE (write_nr, *) "----------"
      WRITE (write_nr, *) matrix_name
      WRITE (write_nr, *) "----------"
      M = SIZE(matrix, 1)
      DO i = 1, M
         WRITE (write_nr, *) matrix(i, :)
      END DO
      WRITE (write_nr, *) "----------"
   END SUBROUTINE write_real_2dmatrix_to

   SUBROUTINE write_real_1dmatrix_to(matrix, write_nr, matrix_name)
      REAL(kind=dp), DIMENSION(:), POINTER               :: matrix
      INTEGER                                            :: write_nr
      CHARACTER(len=*)                                   :: matrix_name

      CHARACTER(len=*), PARAMETER :: routineN = 'write_real_1dmatrix_to', &
         routineP = moduleN//':'//routineN

      WRITE (write_nr, *) "----------"
      WRITE (write_nr, *) matrix_name
      WRITE (write_nr, *) "----------"
      WRITE (write_nr, 101) matrix(:)
101   FORMAT(F15.6)
      WRITE (write_nr, *) "----------"
   END SUBROUTINE write_real_1dmatrix_to

   !Following routine is from  from dbcsr_example3.F
   SUBROUTINE random_dist(dist_array, dist_size, nbins)
      INTEGER, DIMENSION(:), INTENT(out), POINTER        :: dist_array
      INTEGER, INTENT(in)                                :: dist_size, nbins

      INTEGER                                            :: i

      ALLOCATE (dist_array(dist_size))
      DO i = 1, dist_size
         dist_array(i) = MODULO(nbins - i, nbins)
      END DO

   END SUBROUTINE random_dist

END MODULE vab
