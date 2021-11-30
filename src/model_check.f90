module perturb_subs
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  implicit none

  integer, parameter :: nvars = 6
  character(len=500), dimension(nvars), parameter :: model_names = &
    (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", &
                            "reg1/vsh", "reg1/eta", "reg1/rho"/)
  integer, parameter :: vpv_idx = 1, vph_idx = 2, vsv_idx = 3, &
                        vsh_idx = 4, eta_idx = 5, rho_idx = 6

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars) :: ref_model, &
                                                                          new_model
  ! 6 parameter perturbation + 5 extra perturbation
  ! don't change the order unless you know what you are doing
  character(len=500), dimension(11), parameter :: perturb_names = &
    (/character(len=500) :: "reg1/dvpvvpv", "reg1/dvphvph", "reg1/dvsvvsv", &
                            "reg1/dvshvsh", "reg1/detaeta", "reg1/drhorho", &
                            "reg1/dvp",     "reg1/dvs",     "reg1/dvp_vs_ratio", &
                            "reg1/dvsh_vsv_ratio_1", "reg1/dvsh_vsv_ratio_2"/)
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars+5) :: perturb_model

  contains
  subroutine get_sys_args(solver_file)
    use global_var, only : myrank
    use global_var, only : exit_mpi
    implicit none
    character(len=500), intent(in) :: solver_file

    call getarg(1, solver_file)

    if (trim(solver_file) == '') then
      call exit_mpi('Usage: xmodel_check /path/to/solver_data.bp')
    endif

    if (myrank == 0) then
      print*, "Solver file: ", trim(solver_file)
    
    endif
  end subroutine get_sys_args

  ! subroutine calculate_other_perturbation()

  !   use global_var, only : ONE_THIRD, TWO_THIRDS

  !   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: dvp, dvs, &
  !     dvp_vs_ratio, dvsv_vsh_ratio_1, dvsv_vsh_ratio_2

  !   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: vp_0, vp_1, &
  !                                                                 vs_0, vs_1, &
  !                                                                 tmpArr0, tmpArr1

  !   dvsv_vsh_ratio_1 = log( (new_model(:,:,:,:,vsv_idx)/new_model(:,:,:,:,vsh_idx)) / &
  !                           (ref_model(:,:,:,:,vsv_idx)/ref_model(:,:,:,:,vsh_idx)) )

  !   vp_0 = ref_model(:, :, :, :, vpv_idx)
  !   vp_1 = new_model(:, :, :, :, vpv_idx)
  !   dvp = log(vp_1 / vp_0)

  !   vs_0 = sqrt(TWO_THIRDS * ref_model(:, :, :, :, vsv_idx) ** 2 + &
  !               ONE_THIRD  * ref_model(:, :, :, :, vsh_idx) ** 2)
  !   vs_1 = sqrt(TWO_THIRDS * new_model(:, :, :, :, vsv_idx) ** 2 + &
  !               ONE_THIRD  * new_model(:, :, :, :, vsh_idx) ** 2)
  !   dvs = log(vs_1 / vs_0)

  !   dvp_vs_ratio = log( (vp_1/vs_1) / (vp_0/vs_0))

  !   tmpArr0 = (ref_model(:,:,:,:,vsh_idx) - ref_model(:,:,:,:,vsv_idx)) / vs_0
  !   tmpArr1 = (new_model(:,:,:,:,vsh_idx) - new_model(:,:,:,:,vsv_idx)) / vs_1
  !   dvsv_vsh_ratio_2 = tmpArr1 - tmpArr0

  !   perturb_model(:,:,:,:,7) = dvp
  !   perturb_model(:,:,:,:,8) = dvs
  !   perturb_model(:,:,:,:,9) = dvp_vs_ratio
  !   perturb_model(:,:,:,:,10) = dvsv_vsh_ratio_1
  !   perturb_model(:,:,:,:,11) = dvsv_vsh_ratio_2

  ! end subroutine calculate_other_perturbation

end module perturb_subs

program main

  use mpi
  use adios_read_mod
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  use AdiosIO
  use perturb_subs
  implicit none

  character(len=500) :: solver_file

  integer :: ier, i, j, k, ispec, iglob
  real(kind=CUSTOM_REAL) :: threshold, nglob_total, nspec_total, nglob_local, nspec_local
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: checkarray 
  real(kind=CUSTOM_REAL) :: x,y,z,r

  call init_mpi()

  if(myrank == 0) print*, "mpi done"

  call get_sys_args(solver_file)
 
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                                  "verbose=1", ier)

  ! if(myrank == 0) print*, "Read ref model file"
  ! call read_bp_file_real(ref_model_file, model_names, ref_model)
  ! if(myrank == 0) print*, "Read new model file"
  ! call read_bp_file_real(new_model_file, model_names, new_model)
  ! if(myrank == 0) print*, "Done reading"

  ! call read_bp_file_int(input_solver_file, "reg1/idoubling", idoubling)
  call read_bp_file_int(solver_file, "reg1/ibool", ibool)
  call read_bp_file_real(solver_file, "reg1/x_global", x_glob)
  call read_bp_file_real(solver_file, "reg1/y_global", y_glob)
  call read_bp_file_real(solver_file, "reg1/z_global", z_glob)

  perturb_model(:,:,:,:,1:6) = log(new_model / ref_model)
  nglob_local = 0.0
  nspec_local = 0.0
  threshold = 800.0

  do ispec = 1, NSPEC

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i, j, k, ispec)
          x = x_glob(iglob)
          y = y_glob(iglob)
          z = z_glob(iglob)

          r = sqrt(x**2 + y**2 + z**2)

          if (r > 1-800.0/6731.0) then
            checkarray(i,j,k,ispec) = 1.0
          endif
        
        enddo
      enddo
    enddo
    if ( ANY(checkarray(:,:,:,ispec)==1.0) ) then
      nspec_local = nspec_local + 1
      nglob_local = nglob_local + sum(checkarray(:,:,:,ispec)==1.0)
    endif

  enddo

  call sum_all_all_cr(nglob_local, nglob_total)
  call sum_all_all_cr(nspec_local, nspec_total)


  if(myrank == 0) print*, "Total elements:   ", nspec_total
  if(myrank == 0) print*, "Total GLL points: ", nglob_total
  ! call write_bp_file(perturb_model, perturb_names, "KERNELS_GROUP", outputfn)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

  if(myrank == 0) print*, "Job finished"

end program main
