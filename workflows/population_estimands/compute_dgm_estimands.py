"""H1-standard population estimands: archive verification or explicit integration."""
from pathlib import Path
import sys
PACKAGE = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(PACKAGE / 'src/dgm'))
from compute_dgm_estimands import main
if __name__ == '__main__':
    main()
