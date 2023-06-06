!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE ml_ri_perso
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ml_ri_perso'

   PUBLIC :: print_ri_density

CONTAINS

! **************************************************************************************************
!> \brief This is a local modification to CP2K, in which we print quantites related to the RI density,
!>        which is read from a file containing the RI coefficients
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE print_ri_density(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(default_path_length)                     :: mpi_filename
      CHARACTER(default_string_length)                   :: basis_type, filename, header, my_pos_cube
      INTEGER                                            :: i, natom, nsgf_tot, unit_nr
      LOGICAL                                            :: append_cube, do_fft, do_print, &
                                                            file_exists, mpi_io, print_energy, &
                                                            response_only
      REAL(dp)                                           :: ehartree, eqmmm
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: density_coeffs
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: dummy_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rho_g, rho_r, v_hartree_g, v_hartree_r
      TYPE(pw_type), POINTER                             :: rho_core
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: print_section

      NULLIFY (para_env, qs_kind_set, print_section, dummy_fm, atomic_kind_set, cell, dft_control, &
               particle_set, pw_env, auxbas_pw_pool, logger, particles, subsys, poisson_env, energy)

      print_section => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%RI_DENSITY")
      CALL section_vals_get(print_section, explicit=do_print)
      IF (.NOT. do_print) RETURN

      logger => cp_get_default_logger()

      !First we need to be able to read the RI coefficients
      CALL get_qs_env(qs_env, para_env=para_env, natom=natom, qs_kind_set=qs_kind_set, cell=cell, &
                      atomic_kind_set=atomic_kind_set, dft_control=dft_control, pw_env=pw_env, &
                      particle_set=particle_set, subsys=subsys, energy=energy)
      CALL qs_subsys_get(subsys, particles=particles)

      !TODO: generalize
      CPASSERT(dft_control%nspins == 1)

      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%RI_DENSITY%BASIS_TYPE", c_val=basis_type)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf_tot, basis_type=TRIM(basis_type))

      do_fft = .FALSE.
      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%RI_DENSITY%FOURIER_TRANSFORM", l_val=do_fft)

      unit_nr = -1
      file_exists = .FALSE.
      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%RI_DENSITY%RI_COEFFS_FILE", c_val=filename)
      IF (para_env%ionode) THEN
         INQUIRE (FILE=filename, EXIST=file_exists)
         IF (.NOT. file_exists) CPABORT("Trying to read non-existing RI_COEFFS file")

         CALL open_file(file_name=TRIM(filename), file_action="READ", file_form="FORMATTED", &
                        file_position="REWIND", file_status="OLD", unit_number=unit_nr)
      END IF

      ALLOCATE (density_coeffs(nsgf_tot))
      density_coeffs = 0.0_dp
      IF (unit_nr > 0) THEN
         READ (unit_nr, *) header
         DO i = 1, nsgf_tot
            READ (unit_nr, *) density_coeffs(i)
         END DO
      END IF
      CALL mp_sum(density_coeffs, para_env%group)

      IF (para_env%ionode) THEN
         CALL close_file(unit_number=unit_nr)
      END IF

      !Then we need the density on the grid
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_r, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL calculate_wavefunction(dummy_fm, 1, rho_r, rho_g, atomic_kind_set, qs_kind_set, cell, &
                                  dft_control, particle_set, pw_env, TRIM(basis_type), density_coeffs)

      !Print the electronic density to file
      IF (BTEST(cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                           "DFT%PRINT%RI_DENSITY%E_DENSITY_CUBE"), cp_p_file)) THEN

         append_cube = section_get_lval(qs_env%input, "DFT%PRINT%RI_DENSITY%E_DENSITY_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF

         filename = "RI_ELECTRON_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%RI_DENSITY%E_DENSITY_CUBE", &
                                        extension=".cube", middle_name=TRIM(filename), &
                                        file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                        fout=mpi_filename)

         IF (do_fft) THEN
            CALL pw_zero(rho_r)
            CALL pw_transfer(rho_g, rho_r)
         END IF

         CALL cp_pw_to_cube(rho_r, unit_nr, "RI ELECTRON DENSITY", particles=particles, &
                            stride=section_get_ivals(qs_env%input, "DFT%PRINT%RI_DENSITY%E_DENSITY_CUBE%STRIDE"), &
                            mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, qs_env%input, &
                                           "DFT%PRINT%RI_DENSITY%E_DENSITY_CUBE", mpi_io=mpi_io)
      END IF

      !Print the electrostatic potential to a cube file
      IF (BTEST(cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                           "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE"), cp_p_file)) THEN

         CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_r, use_data=REALDATA3D, in_space=REALSPACE)

         response_only = section_get_lval(qs_env%input, "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE%RESPONSE_DENSITY_ONLY")

         NULLIFY (rho_core)
         CALL get_qs_env(qs_env, rho_core=rho_core)
         IF (.NOT. response_only) CALL pw_axpy(rho_core, rho_g)

         CALL pw_poisson_solve(poisson_env, rho_g, ehartree, v_hartree_g)
         CALL pw_transfer(v_hartree_g, v_hartree_r)

         append_cube = section_get_lval(qs_env%input, "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF

         filename = "RI_V_ELECTROSTATIC"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE", &
                                        extension=".cube", middle_name=TRIM(filename), &
                                        file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                        fout=mpi_filename)

         CALL cp_pw_to_cube(v_hartree_r, unit_nr, "ELECTROSTATIC POTENTIAL DUE TO RI DENSITY", particles=particles, &
                            stride=section_get_ivals(qs_env%input, "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE%STRIDE"), &
                            mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, qs_env%input, &
                                           "DFT%PRINT%RI_DENSITY%V_ELECTROSTATIC_CUBE", mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_r)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_g)
      END IF

      !Printing the details of the electrostatic energy
      print_energy = .FALSE.
      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%RI_DENSITY%PRINT_ELECTROSTATIC_ENERGY", l_val=print_energy)

      IF (print_energy) THEN

         unit_nr = cp_logger_get_default_io_unit(logger)

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(/,A,/,A)") &
               " !-----------------------------------------------------------------------------!", &
               " Electrostatic energy from external RI density coefficients [a.u.]:"
         END IF

         !The QM electrostatic total energy (Hartree energy electronic density + ions)
         CALL pw_transfer(rho_r, rho_g)
         NULLIFY (rho_core)
         CALL get_qs_env(qs_env, rho_core=rho_core)
         CALL pw_axpy(rho_core, rho_g)
         CALL pw_poisson_solve(poisson_env, rho_g, ehartree)

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(A,F26.16)") &
               " Hartree energy (from QM electronic density and ions):", ehartree
         END IF

         IF (qs_env%qmmm) THEN
            !The QM/MM electrostatic energy (from QM electron)
            eqmmm = pw_integral_ab(rho_r, qs_env%ks_qmmm_env%v_qmmm_rspace%pw)

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(A, F32.16)") &
                  " QM/MM electrostatic energy (from QM electrons):", eqmmm
            END IF

            !The QM/MM electrostatic energy (from QM ions)
            CALL pw_transfer(rho_core, rho_r)
            eqmmm = pw_integral_ab(rho_r, qs_env%ks_qmmm_env%v_qmmm_rspace%pw)

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(A, F37.16)") &
                  " QM/MM electrostatic energy (from QM ions):", eqmmm
            END IF

         END IF

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(A,/)") &
               " !-----------------------------------------------------------------------------!"
         END IF
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g)

   END SUBROUTINE print_ri_density

END MODULE ml_ri_perso
