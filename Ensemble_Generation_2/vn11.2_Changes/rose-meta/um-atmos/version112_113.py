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



class vn112_tXXXX(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #XXXX by <author>."""

    BEFORE_TAG = "vn11.2"
    AFTER_TAG = "vn11.2_tXXXX"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        return config, self.reports


class vn112_t4751(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #4751 by christophersymonds."""

    BEFORE_TAG = "vn11.2"
    AFTER_TAG = "vn11.2_t4751"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration with ACURE PPE variables."""

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bl_nuc"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bl_nuc"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_ageing"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_ageing"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_accumulation_width"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_accumulation_width"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_aitken_width"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_aitken_width"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_cloudph"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_cloudph"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_ff_ems"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_ff_ems"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bb_ems"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bb_ems"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bf_ems"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bf_ems"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_ff_diameter"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_ff_diameter"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bb_diameter"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bb_diameter"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bf_diameter"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bf_diameter"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_su_parfrac"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_su_parfrac"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_su_diameter"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_su_diameter"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_ss_flux"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_ss_flux"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_anth_so2"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_anth_so2"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_volc_so2"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_volc_so2"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_bio_soa"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_bio_soa"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dms_flux"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dms_flux"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_aitken_ddrydepaer"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_aitken_ddrydepaer"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_accum_ddrydepaer"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_accum_ddrydepaer"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_drydep_so2"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_drydep_so2"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_kappa_oc"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_kappa_oc"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_sigw"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_sigw"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dust_ems"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dust_ems"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_rain_frac"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_rain_frac"],"0.0"

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_rr_thres"],".false."

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_rr_thres"],"0.0"

        return config, self.reports
