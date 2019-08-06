import rose.upgrade
import re
import sys

class UpgradeError(Exception):

    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__



class vn111_tXXXX(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #XXXX by <author>."""

    BEFORE_TAG = "vn11.1"
    AFTER_TAG = "vn11.1_tXXXX"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        return config, self.reports

class vn111_t4751(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #4751 by christophersymonds."""

    BEFORE_TAG = "vn11.1"
    AFTER_TAG = "vn11.1_t4751"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration with ACURE PPE variables."""

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_bl_nuc"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_bl_nuc"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_ait_width"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_ait_width"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_cloud_ph"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_cloud_ph"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_ff_ems"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_1"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_3"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_4"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_5"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_ems_region_6"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_bb_ems"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_1"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_3"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_4"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_5"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_ems_region_6"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_res_ems"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_ems_region_1"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_ems_region_2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_ems_region_3"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_ems_region_4"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_ems_region_5"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_ff_diam"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_ff_diam"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_bb_diam"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_diam"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_res_diam"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_res_diam"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_carb_bb_vertical"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_carb_bb_vertical"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_prim_moc"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_prim_moc"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_prim_so4_diam"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_prim_so4_diam"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_sea_spray"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_sea_spray"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_anth_so2"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anth_so2_region_1"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anth_so2_region_2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anth_so2_region_3"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anth_so2_region_4"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anth_so2_region_5"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_anthropogenic_so2_emission_height"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_anthropogenic_so2_emission_height"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_volc_so2"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_volc_so2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_bvoc_soa"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_bvoc_soa"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_dms"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_dms"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_dry_dep_ait"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_dry_dep_ait"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_dry_dep_acc"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_dry_dep_acc"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_dry_dep_so2"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_dry_dep_so2"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_deposition_rate_substep"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_deposition_rate_substep"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_kappa_oc"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_kappa_oc"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_sig_w"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_sig_w"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_rain_frac"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_rain_frac"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_cloud_ice_thresh"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_cloud_ice_thresh"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_convective_plume_scavenging"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_convective_plume_scavenging"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_bc_ri"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_bc_ri"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_oxidants"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_oxidants"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_in_cloud_so2_oxidisation"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_in_cloud_so2_oxidisation"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_dp_cor_strat"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_dp_cor_strat"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_cloud_drop_spectral_dispersion"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_cloud_drop_spectral_dispersion"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_autoconversion"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_autoconversion"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_ar"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_ar"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_ice_width"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_ice_width"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_mp_dz_scal"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "acure_mp_dz_scal"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_acure_pcalc_index"],".false.")

        self.add_setting(config,["namelist:run_radiation",
                                 "l_acure_bparam"],".false.")

        self.add_setting(config,["namelist:run_radiation",
                                 "l_acure_two_d_fsd_factor"],".false.")

        self.add_setting(config,["namelist:run_cloud",
                                 "l_acure_dbsdtbs_turb_0"],".false.")

        self.add_setting(config,["namelist:run_stochastic",
                                 "l_acure_m_ci"],".false.")

        self.add_setting(config,["namelist:run_stochastic",
                                 "l_acure_a_ent_1_rp"],".false.")

        self.add_setting(config,["namelist:run_precip",
                                 "l_acure_c_r_correl"],".false.")

        self.add_setting(config,["namelist:run_precip",
                                 "l_acure_ai"],".false.")

        return config, self.reports


