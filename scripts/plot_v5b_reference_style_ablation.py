"""
plot_v5b_reference_style_ablation.py
====================================

Legacy entry point -- delegates to the unified early ablation/plume sequence
plotter (fixed coordinate frame, continuous scene).

Run:
    python scripts\\plot_v5b_reference_style_ablation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))

if __name__ == "__main__":
    from plot_v5b_unified_ablation_sequence import main
    raise SystemExit(main())
