!******************ACURE LOGICALS*************************!

LOGICAL :: l_acure_bl_nuc                      = .FALSE.
LOGICAL :: l_acure_ait_width                   = .FALSE.
LOGICAL :: l_acure_cloud_ph                    = .FALSE.
LOGICAL :: l_acure_carb_ff_ems                 = .FALSE.
LOGICAL :: l_acure_carb_bb_ems                 = .FALSE.
LOGICAL :: l_acure_carb_res_ems                = .FALSE.
LOGICAL :: l_acure_anth_so2                    = .FALSE.
LOGICAL :: l_acure_carb_ff_diam                = .FALSE.
LOGICAL :: l_acure_carb_bb_diam                = .FALSE.
LOGICAL :: l_acure_carb_res_diam               = .FALSE.
LOGICAL :: l_acure_prim_moc                    = .FALSE.
LOGICAL :: l_acure_prim_so4_diam               = .FALSE.
LOGICAL :: l_acure_sea_spray                   = .FALSE.
LOGICAL :: l_acure_volc_so2                    = .FALSE.
LOGICAL :: l_acure_bvoc_soa                    = .FALSE.
LOGICAL :: l_acure_dms                         = .FALSE.
LOGICAL :: l_acure_dry_dep_ait                 = .FALSE.
LOGICAL :: l_acure_dry_dep_acc                 = .FALSE.
LOGICAL :: l_acure_dry_dep_so2                 = .FALSE.
LOGICAL :: l_acure_kappa_oc                    = .FALSE.
LOGICAL :: l_acure_sig_w                       = .FALSE.
LOGICAL :: l_acure_cloud_ice_thresh            = .FALSE.
LOGICAL :: l_acure_convective_plume_scavenging = .FALSE.
LOGICAL :: l_acure_scav_diam                   = .FALSE.
LOGICAL :: l_acure_bc_ri                       = .FALSE.
LOGICAL :: l_acure_oxidants_oh                 = .FALSE.
LOGICAL :: l_acure_oxidants_o3                 = .FALSE.
LOGICAL :: l_acure_autoconv_exp_lwp            = .FALSE.
LOGICAL :: l_acure_autoconv_exp_nd             = .FALSE.
LOGICAL :: l_acure_pcalc_index                 = .FALSE.


!*****************End ACURE LOGICALS**********************!

!******************ACURE REALS*************************!

REAL :: acure_bl_nuc                      = rmdi
REAL :: acure_ait_width                   = rmdi
REAL :: acure_cloud_ph                    = rmdi
REAL :: acure_carb_ff_ems_eur             = rmdi
REAL :: acure_carb_ff_ems_nam             = rmdi
REAL :: acure_carb_ff_ems_chi             = rmdi
REAL :: acure_carb_ff_ems_asi             = rmdi
REAL :: acure_carb_ff_ems_mar             = rmdi
REAL :: acure_carb_ff_ems_r               = rmdi
REAL :: acure_carb_bb_ems_sam             = rmdi
REAL :: acure_carb_bb_ems_naf             = rmdi
REAL :: acure_carb_bb_ems_saf             = rmdi
REAL :: acure_carb_bb_ems_bnh             = rmdi
REAL :: acure_carb_bb_ems_rnh             = rmdi
REAL :: acure_carb_bb_ems_rsh             = rmdi
REAL :: acure_carb_res_ems_chi            = rmdi
REAL :: acure_carb_res_ems_asi            = rmdi
REAL :: acure_carb_res_ems_afr            = rmdi
REAL :: acure_carb_res_ems_lat            = rmdi
REAL :: acure_carb_res_ems_r              = rmdi
REAL :: acure_anth_so2_chi                = rmdi
REAL :: acure_anth_so2_asi                = rmdi
REAL :: acure_anth_so2_eur                = rmdi
REAL :: acure_anth_so2_nam                = rmdi
REAL :: acure_anth_so2_r                  = rmdi
REAL :: acure_carb_ff_diam                = rmdi
REAL :: acure_carb_bb_diam                = rmdi
REAL :: acure_carb_res_diam               = rmdi
REAL :: acure_prim_moc                    = rmdi
REAL :: acure_prim_so4_diam               = rmdi
REAL :: acure_sea_spray                   = rmdi
REAL :: acure_volc_so2                    = rmdi
REAL :: acure_bvoc_soa                    = rmdi
REAL :: acure_dms                         = rmdi
REAL :: acure_dry_dep_ait                 = rmdi
REAL :: acure_dry_dep_acc                 = rmdi
REAL :: acure_dry_dep_so2                 = rmdi
REAL :: acure_kappa_oc                    = rmdi
REAL :: acure_sig_w                       = rmdi
REAL :: acure_cloud_ice_thresh            = rmdi
REAL :: acure_convective_plume_scavenging = rmdi
REAL :: acure_scav_diam                   = rmdi
REAL :: acure_bc_ri                       = rmdi
REAL :: acure_oxidants_oh                 = rmdi
REAL :: acure_oxidants_o3                 = rmdi
REAL :: acure_autoconv_exp_lwp            = rmdi
REAL :: acure_autoconv_exp_nd             = rmdi
REAL :: acure_pcalc_index                 = rmdi


!*****************End ACURE REALS**********************!

Insert the following at the end of the namelist block starting "NAMELIST/run_ukca/"

         l_acure_bl_nuc, acure_bl_nuc,                            &
         l_acure_ait_width, acure_ait_width,                      &
         l_acure_cloud_ph, acure_cloud_ph,                        &
         l_acure_carb_ff_ems, acure_carb_ff_ems_eur,              &
         acure_carb_ff_ems_nam, acure_carb_ff_ems_chi,            &
         acure_carb_ff_ems_asi, acure_carb_ff_ems_mar,            &
         acure_carb_ff_ems_r, l_acure_carb_bb_ems,                &
         acure_carb_bb_ems_sam, acure_carb_bb_ems_naf,            &
         acure_carb_bb_ems_saf, acure_carb_bb_ems_bnh,            &
         acure_carb_bb_ems_rnh, acure_carb_bb_ems_rsh,            &
         l_acure_carb_res_ems, acure_carb_res_ems_chi,            &
         acure_carb_res_ems_asi, acure_carb_res_ems_afr,          &
         acure_carb_res_ems_lat, acure_carb_res_ems_r,            &
         l_acure_anth_so2, acure_anth_so2_chi,                    &
         acure_anth_so2_asi, acure_anth_so2_eur,                  &
         acure_anth_so2_nam, acure_anth_so2_r                     &
         l_acure_carb_ff_diam, acure_carb_ff_diam                 &
         l_acure_carb_bb_diam, acure_carb_bb_diam,                &
         l_acure_carb_res_diam, acure_carb_res_diam,              &
         l_acure_prim_moc, acure_prim_moc,                        &
         l_acure_prim_so4_diam, acure_prim_so4_diam,              &
         l_acure_sea_spray, acure_sea_spray,                      &
         l_acure_volc_so2, acure_volc_so2,                        &
         l_acure_bvoc_soa, acure_bvoc_soa,                        &
         l_acure_dms, acure_dms,                                  &
         l_acure_dry_dep_ait, acure_dry_dep_ait,                  &
         l_acure_dry_dep_acc, acure_dry_dep_acc,                  &
         l_acure_dry_dep_so2, acure_dry_dep_so2,                  &
         l_acure_kappa_oc, acure_kappa_oc,                        &
         l_acure_sig_w, acure_sig_w,                              &
         l_acure_cloud_ice_thresh, acure_cloud_ice_thresh,        &
         l_acure_convective_plume_scavenging,                     &
         acure_convective_plume_scavenging,                       &
         l_acure_scav_diam, acure_scav_diam,                      &
         l_acure_bc_ri, acure_bc_ri,                              &
         l_acure_oxidants_oh, acure_oxidants_oh,                  &
         l_acure_oxidants_o3, acure_oxidants_o3,                  &
         l_acure_autoconv_exp_lwp, acure_autoconv_exp_lwp,        &
         l_acure_autoconv_exp_nd, acure_autoconv_exp_nd,          &
         l_acure_pcalc_index


Insert the following before the line "CALL umPrint('- - - - - - end of namelist - - - - - -', &)"

WRITE(lineBuffer,'(A,L1)')'  l_acure_bl_nuc = ',l_acure_bl_nuc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_bl_nuc = ',acure_bl_nuc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_ait_width = ',l_acure_ait_width
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_ait_width = ',acure_ait_width
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_cloud_ph = ',l_acure_cloud_ph
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_cloud_ph = ',acure_cloud_ph
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_ff_ems = ',l_acure_carb_ff_ems
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_eur = ',acure_carb_ff_ems_eur
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_nam = ',acure_carb_ff_ems_nam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_chi = ',acure_carb_ff_ems_chi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_asi = ',acure_carb_ff_ems_asi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_mar = ',acure_carb_ff_ems_mar
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_ems_r = ',acure_carb_ff_ems_r
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_bb_ems = ',l_acure_carb_bb_ems
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_sam = ',acure_carb_bb_ems_sam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_naf = ',acure_carb_bb_ems_naf
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_saf = ',acure_carb_bb_ems_saf
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_bnh = ',acure_carb_bb_ems_bnh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_rnh = ',acure_carb_bb_ems_rnh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_ems_rsh = ',acure_carb_bb_ems_rsh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_res_ems = ',l_acure_carb_res_ems
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_ems_chi = ',acure_carb_res_ems_chi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_ems_asi = ',acure_carb_res_ems_asi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_ems_afr = ',acure_carb_res_ems_afr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_ems_lat = ',acure_carb_res_ems_lat
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_ems_r = ',acure_carb_res_ems_r
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_anth_so2 = ',l_acure_anth_so2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_anth_so2_chi = ',acure_anth_so2_chi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_anth_so2_asi = ',acure_anth_so2_asi
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_anth_so2_eur = ',acure_anth_so2_eur
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_anth_so2_nam = ',acure_anth_so2_nam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_anth_so2_r = ',acure_anth_so2_r
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_ff_diam = ',l_acure_carb_ff_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_ff_diam = ',acure_carb_ff_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_bb_diam = ',l_acure_carb_bb_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_bb_diam = ',acure_carb_bb_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_carb_res_diam = ',l_acure_carb_res_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_carb_res_diam = ',acure_carb_res_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_prim_moc = ',l_acure_prim_moc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_prim_moc = ',acure_prim_moc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_prim_so4_diam = ',l_acure_prim_so4_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_prim_so4_diam = ',acure_prim_so4_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_sea_spray = ',l_acure_sea_spray
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_sea_spray = ',acure_sea_spray
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_volc_so2 = ',l_acure_volc_so2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_volc_so2 = ',acure_volc_so2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_bvoc_soa = ',l_acure_bvoc_soa
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_bvoc_soa = ',acure_bvoc_soa
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_dms = ',l_acure_dms
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_dms = ',acure_dms
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_dry_dep_ait = ',l_acure_dry_dep_ait
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_dry_dep_ait = ',acure_dry_dep_ait
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_dry_dep_acc = ',l_acure_dry_dep_acc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_dry_dep_acc = ',acure_dry_dep_acc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_dry_dep_so2 = ',l_acure_dry_dep_so2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_dry_dep_so2 = ',acure_dry_dep_so2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_kappa_oc = ',l_acure_kappa_oc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_kappa_oc = ',acure_kappa_oc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_sig_w = ',l_acure_sig_w
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_sig_w = ',acure_sig_w
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_cloud_ice_thresh = ', &
                                                        l_acure_cloud_ice_thresh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_cloud_ice_thresh = ',acure_cloud_ice_thresh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_convective_plume_scavenging = ', &
                                             l_acure_convective_plume_scavenging
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_convective_plume_scavenging = ', &
                                               acure_convective_plume_scavenging
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_scav_diam = ',l_acure_scav_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_scav_diam = ',acure_scav_diam
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_bc_ri = ',l_acure_bc_ri
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_bc_ri = ',acure_bc_ri
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_oxidants_oh = ',l_acure_oxidants_oh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_oxidants_oh = ',acure_oxidants_oh
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_oxidants_o3 = ',l_acure_oxidants_o3
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_oxidants_o3 = ',acure_oxidants_o3
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_autoconv_exp_lwp = ', &
                                                        l_acure_autoconv_exp_lwp
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_autoconv_exp_lwp = ',acure_autoconv_exp_lwp
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_autoconv_exp_nd = ',l_acure_autoconv_exp_nd
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' acure_autoconv_exp_nd = ',acure_autoconv_exp_nd
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')'  l_acure_pcalc_index = ',l_acure_pcalc_index
CALL umPrint(lineBuffer,src='ukca_option_mod')


Insert the following in "TYPE my_namelist"

  REAL :: acure_bl_nuc
  REAL :: acure_ait_width
  REAL :: acure_cloud_ph
  REAL :: acure_carb_ff_ems_eur
  REAL :: acure_carb_ff_ems_nam
  REAL :: acure_carb_ff_ems_chi
  REAL :: acure_carb_ff_ems_asi
  REAL :: acure_carb_ff_ems_mar
  REAL :: acure_carb_ff_ems_r
  REAL :: acure_carb_bb_ems_sam
  REAL :: acure_carb_bb_ems_naf
  REAL :: acure_carb_bb_ems_saf
  REAL :: acure_carb_bb_ems_bnh
  REAL :: acure_carb_bb_ems_rnh
  REAL :: acure_carb_bb_ems_rsh
  REAL :: acure_carb_res_ems_chi
  REAL :: acure_carb_res_ems_asi
  REAL :: acure_carb_res_ems_afr
  REAL :: acure_carb_res_ems_lat
  REAL :: acure_carb_res_ems_r
  REAL :: acure_anth_so2_chi
  REAL :: acure_anth_so2_asi
  REAL :: acure_anth_so2_eur
  REAL :: acure_anth_so2_nam
  REAL :: acure_anth_so2_r
  REAL :: acure_carb_ff_diam
  REAL :: acure_carb_bb_diam
  REAL :: acure_carb_res_diam
  REAL :: acure_prim_moc
  REAL :: acure_prim_so4_diam
  REAL :: acure_sea_spray
  REAL :: acure_volc_so2
  REAL :: acure_bvoc_soa
  REAL :: acure_dms
  REAL :: acure_dry_dep_ait
  REAL :: acure_dry_dep_acc
  REAL :: acure_dry_dep_so2
  REAL :: acure_kappa_oc
  REAL :: acure_sig_w
  REAL :: acure_cloud_ice_thresh
  REAL :: acure_convective_plume_scavenging
  REAL :: acure_scav_diam
  REAL :: acure_bc_ri
  REAL :: acure_oxidants_oh
  REAL :: acure_oxidants_o3
  REAL :: acure_autoconv_exp_lwp
  REAL :: acure_autoconv_exp_nd


  LOGICAL :: l_acure_bl_nuc
  LOGICAL :: l_acure_ait_width
  LOGICAL :: l_acure_cloud_ph
  LOGICAL :: l_acure_carb_ff_ems
  LOGICAL :: l_acure_carb_bb_ems
  LOGICAL :: l_acure_carb_res_ems
  LOGICAL :: l_acure_anth_so2
  LOGICAL :: l_acure_carb_ff_diam
  LOGICAL :: l_acure_carb_bb_diam
  LOGICAL :: l_acure_carb_res_diam
  LOGICAL :: l_acure_prim_moc
  LOGICAL :: l_acure_prim_so4_diam
  LOGICAL :: l_acure_sea_spray
  LOGICAL :: l_acure_volc_so2
  LOGICAL :: l_acure_bvoc_soa
  LOGICAL :: l_acure_dms
  LOGICAL :: l_acure_dry_dep_ait
  LOGICAL :: l_acure_dry_dep_acc
  LOGICAL :: l_acure_dry_dep_so2
  LOGICAL :: l_acure_kappa_oc
  LOGICAL :: l_acure_sig_w
  LOGICAL :: l_acure_cloud_ice_thresh
  LOGICAL :: l_acure_convective_plume_scavenging
  LOGICAL :: l_acure_scav_diam
  LOGICAL :: l_acure_bc_ri
  LOGICAL :: l_acure_oxidants_oh
  LOGICAL :: l_acure_oxidants_o3
  LOGICAL :: l_acure_autoconv_exp_lwp
  LOGICAL :: l_acure_autoconv_exp_nd
  LOGICAL :: l_acure_pcalc_index

After "IF (mype == 0) THEN", insert

  my_nml % acure_bl_nuc           = acure_bl_nuc
  my_nml % acure_ait_width        = acure_ait_width
  my_nml % acure_cloud_ph         = acure_cloud_ph
  my_nml % acure_carb_ff_ems_eur  = acure_carb_ff_ems_eur
  my_nml % acure_carb_ff_ems_nam  = acure_carb_ff_ems_nam
  my_nml % acure_carb_ff_ems_chi  = acure_carb_ff_ems_chi
  my_nml % acure_carb_ff_ems_asi  = acure_carb_ff_ems_asi
  my_nml % acure_carb_ff_ems_mar  = acure_carb_ff_ems_mar
  my_nml % acure_carb_ff_ems_r    = acure_carb_ff_ems_r
  my_nml % acure_carb_bb_ems_sam  = acure_carb_bb_ems_sam
  my_nml % acure_carb_bb_ems_naf  = acure_carb_bb_ems_naf
  my_nml % acure_carb_bb_ems_saf  = acure_carb_bb_ems_saf
  my_nml % acure_carb_bb_ems_bnh  = acure_carb_bb_ems_bnh
  my_nml % acure_carb_bb_ems_rnh  = acure_carb_bb_ems_rnh
  my_nml % acure_carb_bb_ems_rsh  = acure_carb_bb_ems_rsh
  my_nml % acure_carb_res_ems_chi = acure_carb_res_ems_chi
  my_nml % acure_carb_res_ems_asi = acure_carb_res_ems_asi
  my_nml % acure_carb_res_ems_afr = acure_carb_res_ems_afr
  my_nml % acure_carb_res_ems_lat = acure_carb_res_ems_lat
  my_nml % acure_carb_res_ems_r   = acure_carb_res_ems_r
  my_nml % acure_anth_so2_chi     = acure_anth_so2_chi
  my_nml % acure_anth_so2_asi     = acure_anth_so2_asi
  my_nml % acure_anth_so2_eur     = acure_anth_so2_eur
  my_nml % acure_anth_so2_nam     = acure_anth_so2_nam
  my_nml % acure_anth_so2_r       = acure_anth_so2_r
  my_nml % acure_carb_ff_diam     = acure_carb_ff_diam
  my_nml % acure_carb_bb_diam     = acure_carb_bb_diam
  my_nml % acure_carb_res_diam    = acure_carb_res_diam
  my_nml % acure_prim_moc         = acure_prim_moc
  my_nml % acure_prim_so4_diam    = acure_prim_so4_diam
  my_nml % acure_sea_spray        = acure_sea_spray
  my_nml % acure_volc_so2         = acure_volc_so2
  my_nml % acure_bvoc_soa         = acure_bvoc_soa
  my_nml % acure_dms              = acure_dms
  my_nml % acure_dry_dep_ait      = acure_dry_dep_ait
  my_nml % acure_dry_dep_acc      = acure_dry_dep_acc
  my_nml % acure_dry_dep_so2      = acure_dry_dep_so2
  my_nml % acure_kappa_oc         = acure_kappa_oc
  my_nml % acure_sig_w            = acure_sig_w
  my_nml % acure_cloud_ice_thresh = acure_cloud_ice_thresh
  my_nml % acure_convective_plume_scavenging = acure_convective_plume_scavenging
  my_nml % acure_scav_diam        = acure_scav_diam
  my_nml % acure_bc_ri            = acure_bc_ri
  my_nml % acure_oxidants_oh      = acure_oxidants_oh
  my_nml % acure_oxidants_o3      = acure_oxidants_o3
  my_nml % acure_autoconv_exp_lwp = acure_autoconv_exp_lwp
  my_nml % acure_autoconv_exp_nd  = acure_autoconv_exp_nd


  my_nml % l_acure_bl_nuc           = l_acure_bl_nuc
  my_nml % l_acure_ait_width        = l_acure_ait_width
  my_nml % l_acure_cloud_ph         = l_acure_cloud_ph
  my_nml % l_acure_carb_ff_ems      = l_acure_carb_ff_ems
  my_nml % l_acure_carb_bb_ems      = l_acure_carb_bb_ems
  my_nml % l_acure_carb_res_ems     = l_acure_carb_res_ems
  my_nml % l_acure_anth_so2         = l_acure_anth_so2
  my_nml % l_acure_carb_ff_diam     = l_acure_carb_ff_diam
  my_nml % l_acure_carb_bb_diam     = l_acure_carb_bb_diam
  my_nml % l_acure_carb_res_diam    = l_acure_carb_res_diam
  my_nml % l_acure_prim_moc         = l_acure_prim_moc
  my_nml % l_acure_prim_so4_diam    = l_acure_prim_so4_diam
  my_nml % l_acure_sea_spray        = l_acure_sea_spray
  my_nml % l_acure_volc_so2         = l_acure_volc_so2
  my_nml % l_acure_bvoc_soa         = l_acure_bvoc_soa
  my_nml % l_acure_dms              = l_acure_dms
  my_nml % l_acure_dry_dep_ait      = l_acure_dry_dep_ait
  my_nml % l_acure_dry_dep_acc      = l_acure_dry_dep_acc
  my_nml % l_acure_dry_dep_so2      = l_acure_dry_dep_so2
  my_nml % l_acure_kappa_oc         = l_acure_kappa_oc
  my_nml % l_acure_sig_w            = l_acure_sig_w
  my_nml % l_acure_cloud_ice_thresh = l_acure_cloud_ice_thresh
  my_nml % l_acure_convective_plume_scavenging = l_acure_convective_plume_scavenging
  my_nml % l_acure_scav_diam        = l_acure_scav_diam
  my_nml % l_acure_bc_ri            = l_acure_bc_ri
  my_nml % l_acure_oxidants_oh      = l_acure_oxidants_oh
  my_nml % l_acure_oxidants_o3      = l_acure_oxidants_o3
  my_nml % l_acure_autoconv_exp_lwp = l_acure_autoconv_exp_lwp
  my_nml % l_acure_autoconv_exp_nd  = l_acure_autoconv_exp_nd
  my_nml % l_acure_pcalc_index      = l_acure_pcalc_index

After "IF (mype /= 0) THEN", insert

  acure_bl_nuc           = my_nml % acure_bl_nuc
  acure_ait_width        = my_nml % acure_ait_width
  acure_cloud_ph         = my_nml % acure_cloud_ph
  acure_carb_ff_ems_eur  = my_nml % acure_carb_ff_ems_eur
  acure_carb_ff_ems_nam  = my_nml % acure_carb_ff_ems_nam
  acure_carb_ff_ems_chi  = my_nml % acure_carb_ff_ems_chi
  acure_carb_ff_ems_asi  = my_nml % acure_carb_ff_ems_asi
  acure_carb_ff_ems_mar  = my_nml % acure_carb_ff_ems_mar
  acure_carb_ff_ems_r    = my_nml % acure_carb_ff_ems_r
  acure_carb_bb_ems_sam  = my_nml % acure_carb_bb_ems_sam
  acure_carb_bb_ems_naf  = my_nml % acure_carb_bb_ems_naf
  acure_carb_bb_ems_saf  = my_nml % acure_carb_bb_ems_saf
  acure_carb_bb_ems_bnh  = my_nml % acure_carb_bb_ems_bnh
  acure_carb_bb_ems_rnh  = my_nml % acure_carb_bb_ems_rnh
  acure_carb_bb_ems_rsh  = my_nml % acure_carb_bb_ems_rsh
  acure_carb_res_ems_chi = my_nml % acure_carb_res_ems_chi
  acure_carb_res_ems_asi = my_nml % acure_carb_res_ems_asi
  acure_carb_res_ems_afr = my_nml % acure_carb_res_ems_afr
  acure_carb_res_ems_lat = my_nml % acure_carb_res_ems_lat
  acure_carb_res_ems_r   = my_nml % acure_carb_res_ems_r
  acure_anth_so2_chi     = my_nml % acure_anth_so2_chi
  acure_anth_so2_asi     = my_nml % acure_anth_so2_asi
  acure_anth_so2_eur     = my_nml % acure_anth_so2_eur
  acure_anth_so2_nam     = my_nml % acure_anth_so2_nam
  acure_anth_so2_r       = my_nml % acure_anth_so2_r
  acure_carb_ff_diam     = my_nml % acure_carb_ff_diam
  acure_carb_bb_diam     = my_nml % acure_carb_bb_diam
  acure_carb_res_diam    = my_nml % acure_carb_res_diam
  acure_prim_moc         = my_nml % acure_prim_moc
  acure_prim_so4_diam    = my_nml % acure_prim_so4_diam
  acure_sea_spray        = my_nml % acure_sea_spray
  acure_volc_so2         = my_nml % acure_volc_so2
  acure_bvoc_soa         = my_nml % acure_bvoc_soa
  acure_dms              = my_nml % acure_dms
  acure_dry_dep_ait      = my_nml % acure_dry_dep_ait
  acure_dry_dep_acc      = my_nml % acure_dry_dep_acc
  acure_dry_dep_so2      = my_nml % acure_dry_dep_so2
  acure_kappa_oc         = my_nml % acure_kappa_oc
  acure_sig_w            = my_nml % acure_sig_w
  acure_cloud_ice_thresh = my_nml % acure_cloud_ice_thresh
  acure_convective_plume_scavenging = my_nml % acure_convective_plume_scavenging
  acure_scav_diam        = my_nml % acure_scav_diam
  acure_bc_ri            = my_nml % acure_bc_ri
  acure_oxidants_oh      = my_nml % acure_oxidants_oh
  acure_oxidants_o3      = my_nml % acure_oxidants_o3
  acure_autoconv_exp_lwp = my_nml % acure_autoconv_exp_lwp
  acure_autoconv_exp_nd  = my_nml % acure_autoconv_exp_nd


  l_acure_bl_nuc           = my_nml % l_acure_bl_nuc
  l_acure_ait_width        = my_nml % l_acure_ait_width
  l_acure_cloud_ph         = my_nml % l_acure_cloud_ph
  l_acure_carb_ff_ems      = my_nml % l_acure_carb_ff_ems
  l_acure_carb_bb_ems      = my_nml % l_acure_carb_bb_ems
  l_acure_carb_res_ems     = my_nml % l_acure_carb_res_ems
  l_acure_anth_so2         = my_nml % l_acure_anth_so2
  l_acure_carb_ff_diam     = my_nml % l_acure_carb_ff_diam
  l_acure_carb_bb_diam     = my_nml % l_acure_carb_bb_diam
  l_acure_carb_res_diam    = my_nml % l_acure_carb_res_diam
  l_acure_prim_moc         = my_nml % l_acure_prim_moc
  l_acure_prim_so4_diam    = my_nml % l_acure_prim_so4_diam
  l_acure_sea_spray        = my_nml % l_acure_sea_spray
  l_acure_volc_so2         = my_nml % l_acure_volc_so2
  l_acure_bvoc_soa         = my_nml % l_acure_bvoc_soa
  l_acure_dms              = my_nml % l_acure_dms
  l_acure_dry_dep_ait      = my_nml % l_acure_dry_dep_ait
  l_acure_dry_dep_acc      = my_nml % l_acure_dry_dep_acc
  l_acure_dry_dep_so2      = my_nml % l_acure_dry_dep_so2
  l_acure_kappa_oc         = my_nml % l_acure_kappa_oc
  l_acure_sig_w            = my_nml % l_acure_sig_w
  l_acure_cloud_ice_thresh = my_nml % l_acure_cloud_ice_thresh
  l_acure_convective_plume_scavenging = my_nml % l_acure_convective_plume_scavenging
  l_acure_scav_diam        = my_nml % l_acure_scav_diam
  l_acure_bc_ri            = my_nml % l_acure_bc_ri
  l_acure_oxidants_oh      = my_nml % l_acure_oxidants_oh
  l_acure_oxidants_o3      = my_nml % l_acure_oxidants_o3
  l_acure_autoconv_exp_lwp = my_nml % l_acure_autoconv_exp_lwp
  l_acure_autoconv_exp_nd  = my_nml % l_acure_autoconv_exp_nd
  l_acure_pcalc_index      = my_nml % l_acure_pcalc_index
