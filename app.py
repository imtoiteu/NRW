"""NRW Material Parameter Extraction Laboratory - Streamlit entry point.

    export NRW_APP_PASSWORD='...'
    streamlit run app.py --server.address 0.0.0.0 --server.port 8501

All logic lives in src/nrw/web (UI) and src/nrw (tested numerical core).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))  # works without `pip install -e .`

from nrw.web.ui import main  # noqa: E402

main()
