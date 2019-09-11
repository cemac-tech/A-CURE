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

class vn112_t4751(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #4751 by christophersymonds."""

    BEFORE_TAG = "vn11.2"
    AFTER_TAG = "vn11.2_t4751"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration with ACURE PPE variables."""

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_01"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_01"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_02"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_02"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_03"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_03"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_04"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_04"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_05"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_05"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_06"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_06"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_07"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_07"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_08"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_08"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_09"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_09"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_10"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_10"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_11"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_11"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_12"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_12"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_13"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_13"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_14"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_14"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_15"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_15"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_16"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_16"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_17"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_17"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_18"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_18"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_19"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_19"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_20"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_20"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_21"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_21"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_22"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_22"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_23"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_23"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_24"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_24"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_25"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_25"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_26"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_26"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_27"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_27"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_28"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_28"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_29"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_29"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_30"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_30"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_31"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_31"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_32"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_32"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_33"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_33"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_34"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_34"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_35"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_35"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_36"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_36"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_37"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_37"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_38"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_38"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_39"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_39"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_40"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_40"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_41"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_41"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_42"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_42"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_43"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_43"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_44"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_44"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_45"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_45"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_46"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_46"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_47"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_47"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_48"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_48"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_49"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_49"],"0.0")

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_acure_dummy_50"],".false.")

        self.add_setting(config,["namelist:run_ukca",
                                 "ukca_acure_dummy_50"],"0.0")

        return config, self.reports
