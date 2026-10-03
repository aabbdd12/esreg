"""The data folder: data/ beside python/ in the public repository
(replication/data), or data/ at the root of the development project."""
import os

_here = os.path.dirname(os.path.abspath(__file__))
_cands = [os.path.join(_here, os.pardir, "data"),
          os.path.join(_here, os.pardir, os.pardir, "data")]
DATA = os.path.normpath(next((d for d in _cands if os.path.isdir(d)), _cands[0])) + os.sep
