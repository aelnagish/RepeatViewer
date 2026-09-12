# Collect .py source so the Windows runtime hook can patch
# `del obj` in scipy.stats._distn_infrastructure (Python 3.12.0 + PyInstaller).
module_collection_mode = "pyz+py"
