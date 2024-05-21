!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   
! **************************************************************************************************
MODULE cp_dbcsr_multiplication
   USE mathlib, ONLY: lcm, gcd
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_cfm_types, ONLY: cp_cfm_type
   USE dbcsr_api, ONLY: dbcsr_distribution_get, &
                        dbcsr_convert_sizes_to_offsets, dbcsr_add, &
                        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, &
                        dbcsr_deallocate_matrix, &
                        dbcsr_desymmetrize, dbcsr_distribution_new, &
                        dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_matrix_type, &
                        dbcsr_iterator_type, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
                        dbcsr_iterator_start, dbcsr_iterator_stop, &
                        dbcsr_multiply, dbcsr_norm, dbcsr_p_type, dbcsr_release, &
                        dbcsr_reserve_all_blocks, dbcsr_scale, dbcsr_type, &
                        dbcsr_valid_index, dbcsr_verify_matrix, &
                        dbcsr_distribution_type, dbcsr_distribution_release, &
                        dbcsr_norm_frobenius, &
                        dbcsr_type_antisymmetric, dbcsr_type_complex_8, dbcsr_type_no_symmetry, dbcsr_type_real_8, &
                        dbcsr_type_symmetric, dbcsr_type_real_8
   USE dbt_api, ONLY: dbt_create, dbt_destroy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_type, &
                      dbt_contract
   USE cp_fm_basic_linalg, ONLY: cp_fm_gemm
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_info, &
                          cp_fm_release, &
                          cp_fm_to_fm, &
                          cp_fm_type
   USE message_passing, ONLY: mp_para_env_type
   USE distribution_2d_types, ONLY: distribution_2d_get, &
                                    distribution_2d_type
   USE kinds, ONLY: dp, default_string_length, int_8
   USE message_passing, ONLY: mp_comm_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_multiplication'

   PUBLIC :: cp_dbcsr_multiply


CONTAINS

! **************************************************************************************************
!> \brief A test, offloading dbcsr_multiply to dbt_contract, to see
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                first_row, last_row, first_column, last_column, first_k, last_k, &
                                retain_sparsity, filter_eps, flop)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                          :: beta
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                      :: first_row, last_row, first_column, &
                                                            last_column, first_k, last_k
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL         :: flop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(1)                              :: c1, c2, nc1, nc2
      TYPE(dbt_type)                                     :: tensor_a, tensor_b, tensor_c

      CALL timeset(routineN, handle)

      IF (PRESENT(first_row) .OR. PRESENT(last_row) .OR. PRESENT(first_column) &
          .OR. PRESENT(last_column) .OR. PRESENT(first_k) .OR. PRESENT(last_k)) THEN

         !Simply call dbcsr_multiply when these options are present
         !TODO: might be able to use dbt_contract bounds for this
         CALL dbcsr_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                             first_row, last_row, first_column, last_column, first_k, last_k, &
                             retain_sparsity, filter_eps, flop) 

      ELSE IF (dbcsr_get_data_type(matrix_a) .NE. dbcsr_type_real_8) THEN

         !We only offload double precision floats to DBT
         CALL dbcsr_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                             first_row, last_row, first_column, last_column, first_k, last_k, &
                             retain_sparsity, filter_eps, flop)

      ELSE

         CALL dbt_create(matrix_a, tensor_a)
         CALL dbt_create(matrix_b, tensor_b)
         CALL dbt_create(matrix_c, tensor_c)

         CALL dbt_copy_matrix_to_tensor(matrix_a, tensor_a)
         CALL dbt_copy_matrix_to_tensor(matrix_b, tensor_b)

         c1(1) = 2
         nc1(1) = 1
         IF (transa == "T") THEN !TODO: use DBCSR transpose type
            c1(1) = 1
            nc1(1) = 2
         END IF

         c2(1) = 1
         nc2(1) = 2
         IF (transb == "T") THEN
            c2(1) = 2
            nc2(1) = 1
         END IF

         CALL dbt_contract(alpha, tensor_a, tensor_b, beta, tensor_c, &
                           contract_1=c1, notcontract_1=nc1, &
                           contract_2=c2, notcontract_2=nc2, &
                           map_1=[1], map_2=[2], filter_eps=filter_eps, &
                           flop=flop)

         CALL dbt_copy_tensor_to_matrix(tensor_c, matrix_c)

         CALL dbt_destroy(tensor_a)
         CALL dbt_destroy(tensor_b)
         CALL dbt_destroy(tensor_c)

      END IF
      CALL timestop(handle)
   END SUBROUTINE cp_dbcsr_multiply

END MODULE cp_dbcsr_multiplication
