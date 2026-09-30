"""Session-scoped state for the Streamlit app.

Uploaded data and results live only in the user's own Streamlit session (server memory);
they are never written to a web-served location and are never shared between sessions
(no global st.cache_* is used for user data). ``clear_session_data`` removes them.
"""
from __future__ import annotations

from collections import OrderedDict

DATA_KEYS = ("data", "input_name", "input_bytes", "input_summary", "demo_meta", "results", "current_key",
             "comparison", "mc_key", "exports", "fig_memo",
             "data2", "input2_name", "input2_bytes", "input2_summary", "demo2_meta")
MAX_CACHED_RESULTS = 6


def init(ss) -> None:
    ss.setdefault("authenticated", False)
    ss.setdefault("results", OrderedDict())


def clear_session_data(ss) -> None:
    """Forget uploaded data, results and prepared exports (keeps the login)."""
    for k in DATA_KEYS:
        if k in ss:
            del ss[k]
    ss["results"] = OrderedDict()


def logout(ss) -> None:
    clear_session_data(ss)
    ss["authenticated"] = False


def remember(ss, key: str, outcome) -> None:
    """Session-local memo of analysis outcomes (bounded)."""
    res = ss.setdefault("results", OrderedDict())
    res[key] = outcome
    res.move_to_end(key)
    while len(res) > MAX_CACHED_RESULTS:
        res.popitem(last=False)


def recall(ss, key: str):
    res = ss.get("results")
    if res and key in res:
        res.move_to_end(key)
        return res[key]
    return None
