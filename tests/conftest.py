import sys
from unittest.mock import MagicMock

sys.modules['deeplc'] = MagicMock()
sys.modules['deeplc.feat_extractor'] = MagicMock()
sys.modules['deeplc.deeplc'] = MagicMock()
