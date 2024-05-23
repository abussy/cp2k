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
                        dbcsr_deallocate_matrix, dbcsr_get_stored_coordinates, &
                        dbcsr_desymmetrize, dbcsr_distribution_new, dbcsr_clear,&
                        dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_matrix_type, &
                        dbcsr_iterator_type, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
                        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_has_symmetry, &
                        dbcsr_multiply, dbcsr_norm, dbcsr_p_type, dbcsr_release, &
                        dbcsr_reserve_all_blocks, dbcsr_scale, dbcsr_type, dbcsr_reserve_blocks,&
                        dbcsr_valid_index, dbcsr_verify_matrix, dbcsr_finalize, &
                        dbcsr_distribution_type, dbcsr_distribution_release, &
                        dbcsr_norm_frobenius, dbcsr_transpose, dbcsr_get_block_p, &
                        dbcsr_type_antisymmetric, dbcsr_type_complex_8, dbcsr_type_no_symmetry, dbcsr_type_real_8, &
                        dbcsr_type_symmetric, dbcsr_type_real_8, dbcsr_put_block
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
   USE message_passing, ONLY: mp_para_env_type, mp_cart_type, mp_comm_type
   USE distribution_2d_types, ONLY: distribution_2d_get, &
                                    distribution_2d_type
   USE kinds, ONLY: dp, default_string_length, int_8
   USE dbm_api, ONLY: dbm_type, dbm_distribution_new, dbm_distribution_obj, dbm_distribution_release, &
                      dbm_iterator, dbm_iterator_stop, dbm_iterator_next_block, dbm_iterator_start, &
                      dbm_iterator_blocks_left, dbm_release, dbm_create, dbm_multiply, dbm_get_block_p, &
                      dbm_put_block, dbm_clear, dbm_finalize, dbm_get_stored_coordinates, dbm_get_num_blocks

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

      INTEGER                                            :: handle, group, pdims(2)
      TYPE(dbm_type)                                     :: dbm_a, dbm_b, dbm_c
      LOGICAL                                            :: dbm_transa, dbm_transb
      TYPE(mp_comm_type)                                 :: mp_comm
      TYPE(mp_cart_type)                                 :: dbm_comm

      CALL timeset(routineN, handle)

      IF (PRESENT(first_row) .OR. PRESENT(last_row) .OR. PRESENT(first_column) &
          .OR. PRESENT(last_column) .OR. PRESENT(first_k) .OR. PRESENT(last_k)) THEN

         !Simply call dbcsr_multiply when these options are present
         CALL dbcsr_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                             first_row, last_row, first_column, last_column, first_k, last_k, &
                             retain_sparsity, filter_eps, flop) 

      ELSE

         CALL dbcsr_get_info(matrix_a, group=group)
         CALL mp_comm%set_handle(group)

         pdims(:) = 0
         CALL dbm_comm%create(mp_comm, 2, pdims)

         CALL create_dbm_from_dbcsr(matrix_a, dbm_a, dbm_comm)
         CALL create_dbm_from_dbcsr(matrix_b, dbm_b, dbm_comm)
         CALL create_dbm_from_dbcsr(matrix_c, dbm_c, dbm_comm)

         CALL copy_dbcsr_to_dbm(matrix_a, dbm_a)
         CALL copy_dbcsr_to_dbm(matrix_b, dbm_b)
         !TODO: might be more efficient to do a dbcsr_add at the end, and avoid desymmetrizing 
         !      or do a dbcsr_scale and a copy with summation 
         IF (beta > EPSILON(0.0_dp) .OR. PRESENT(retain_sparsity)) CALL copy_dbcsr_to_dbm(matrix_c, dbm_c)

         dbm_transa = transa == dbcsr_transpose
         dbm_transb = transb == dbcsr_transpose

         CALL dbm_multiply(dbm_transa, dbm_transb, alpha, dbm_a, dbm_b, beta, dbm_c, &
                           retain_sparsity=retain_sparsity, filter_eps=filter_eps, flop=flop)

         CALL copy_dbm_to_dbcsr(dbm_c, matrix_c)

         CALL dbm_release(dbm_a)
         CALL dbm_release(dbm_b)
         CALL dbm_release(dbm_c)

         CALL dbm_comm%free()

      END IF
      CALL timestop(handle)
   END SUBROUTINE cp_dbcsr_multiply

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE create_dbm_from_dbcsr(dbcsr_mat, dbm_mat, dbm_comm)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbm_type), INTENT(INOUT)                      :: dbm_mat
      TYPE(mp_cart_type), INTENT(IN)                     :: dbm_comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbm_from_dbcsr'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: row_dist, col_dist, row_blk_size, col_blk_size
      TYPE(dbm_distribution_obj)                         :: dbm_dist

      NULLIFY(row_dist, col_dist, row_blk_size, col_blk_size)

      CALL timeset(routineN, handle)

      !Get distribution and block size from DBCSR matrix, create a DBM matrix with the same
      CALL dbcsr_get_info(dbcsr_mat, row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                          proc_row_dist=row_dist, proc_col_dist=col_dist)

      CALL dbm_distribution_new(dbm_dist, dbm_comm, row_dist, col_dist)

      CALL dbm_create(dbm_mat, "dbm_mat", dbm_dist, row_blk_size, col_blk_size)

      CALL dbm_distribution_release(dbm_dist)

      CALL timestop(handle)

   END SUBROUTINE create_dbm_from_dbcsr

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE copy_dbm_to_dbcsr(dbm_mat, dbcsr_mat)
      TYPE(dbm_type), INTENT(INOUT)                      :: dbm_mat
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat

      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbm_to_dbcsr'

      TYPE(dbm_iterator)                                 :: iter
      INTEGER                                            :: row, col, iblk, nblk, handle
      REAL(dp), DIMENSION(:, :), POINTER                 :: blk
      LOGICAL                                            :: has_symmetry
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: row_idx, col_idx 

      NULLIFY(blk)

      CALL timeset(routineN, handle)

      CALL dbcsr_clear(dbcsr_mat)
      has_symmetry = dbcsr_has_symmetry(dbcsr_mat)
      nblk = dbm_get_num_blocks(dbm_mat)
      ALLOCATE(row_idx(nblk), col_idx(nblk))
      iblk = 0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(dbm_mat,dbcsr_mat,has_symmetry,iblk,nblk,row_idx,col_idx) &
!$OMP PRIVATE(iter,row,col,blk)
      CALL dbm_iterator_start(iter, dbm_mat)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, row=row, column=col)
         CALL dbm_get_block_p(dbm_mat, row, col, blk)

         IF (.NOT. ASSOCIATED(blk)) CYCLE
         IF (has_symmetry) THEN
            IF (checker_tr(row, col)) CYCLE                                              
            IF (row > col) CALL swap(row, col)   
         END IF

!$OMP CRITICAL
         iblk = iblk + 1
         row_idx(iblk) = row
         col_idx(iblk) = col
!$OMP END CRITICAL

      END DO
      CALL dbm_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_reserve_blocks(dbcsr_mat, row_idx(:iblk), col_idx(:iblk))
      CALL dbcsr_finalize(dbcsr_mat)
      
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(dbm_mat,dbcsr_mat,has_symmetry) &
!$OMP PRIVATE(iter,row,col,blk)
      CALL dbm_iterator_start(iter, dbm_mat)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, row=row, column=col)
         CALL dbm_get_block_p(dbm_mat, row, col, blk)
         
         IF (has_symmetry .AND. checker_tr(row, col)) CYCLE
         IF (.NOT. ASSOCIATED(blk)) CYCLE
         IF (has_symmetry .AND. row > col) THEN
            CALL dbcsr_put_block(dbcsr_mat, col, row, TRANSPOSE(blk))
         ELSE
            CALL dbcsr_put_block(dbcsr_mat, row, col, blk)
         END IF
      END DO
      CALL dbm_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE copy_dbm_to_dbcsr
 
! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE copy_dbcsr_to_dbm(dbcsr_mat, dbm_mat)
      TYPE(dbcsr_type), TARGET, INTENT(IN)               :: dbcsr_mat
      TYPE(dbm_type), INTENT(INOUT)                      :: dbm_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbm_from_dbcsr'

      TYPE(dbcsr_iterator_type)                          :: iter
      INTEGER                                            :: row, col, handle
      REAL(dp), DIMENSION(:, :), POINTER                 :: blk
      TYPE(dbcsr_type), POINTER                          :: dbcsr_desym

      NULLIFY(blk, dbcsr_desym)

      CALL timeset(routineN, handle)

      IF (dbcsr_has_symmetry(dbcsr_mat)) THEN
         ALLOCATE (dbcsr_desym)
         CALL dbcsr_desymmetrize(dbcsr_mat, dbcsr_desym)
      ELSE
         dbcsr_desym => dbcsr_mat
      END IF

      
      CALL dbcsr_iterator_start(iter, dbcsr_desym)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row=row, column=col, block=blk)

         IF (.NOT. ASSOCIATED(blk)) CYCLE

         CALL dbm_put_block(dbm_mat, row, col, blk)
      END DO
      CALL dbcsr_iterator_stop(iter)

      IF (dbcsr_has_symmetry(dbcsr_mat)) THEN
         CALL dbcsr_release(dbcsr_desym)
         DEALLOCATE (dbcsr_desym)
      END IF

      CALL dbm_finalize(dbm_mat)

      CALL timestop(handle)

   END SUBROUTINE copy_dbcsr_to_dbm

! ************************************************************************************************** 
! ************************************************************************************************** 
   ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)                                      
      INTEGER, INTENT(IN)                                :: row, column                              
      LOGICAL                                            :: transpose                                
                                                                                                     
      transpose = BTEST(column + row, 0) .EQV. column .GE. row                                       
                                                                                                     
   END FUNCTION checker_tr

! ************************************************************************************************** 
! ************************************************************************************************** 
   ELEMENTAL SUBROUTINE swap(a, b)                                                                   
      INTEGER, INTENT(INOUT)                             :: a, b                                     
      INTEGER                                            :: tmp                                      
                                                                                                     
      tmp = a                                                                                        
      a = b                                                                                          
      b = tmp                                                                                        
   END SUBROUTINE swap

END MODULE cp_dbcsr_multiplication
