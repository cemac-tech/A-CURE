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

# class vn112_tXXXX(rose.upgrade.MacroUpgrade):

#    """Upgrade macro for ticket #XXXX by <author>."""

#    BEFORE_TAG = "vn11.2"
#    AFTER_TAG = "vn11.2_tXXXX"
#
#    def upgrade(self, config, meta_config=None):
#        """Upgrade a UM runtime app configuration."""
#        # Input your macro commands here
#        return config, self.reports
#

class vn112_t4653(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4653 by christophersymonds."""

    BEFORE_TAG = "vn11.2"
    AFTER_TAG = "vn11.2_t4653"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_aeros_volc_so2"],".false.")

        self.add_setting(config,["namelist:run_ukca",
				 "ukca_aeros_volc_so2"],"1.0")

        return config, self.reports

