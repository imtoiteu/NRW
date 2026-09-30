"""Access control for the web application.

* The password is read ONLY from the environment variable ``NRW_APP_PASSWORD`` (never from
  the repository). If it is unset or empty the application refuses to serve any data page.
* Comparison uses ``hmac.compare_digest`` on SHA-256 digests (constant time, no length leak).
* Failed attempts are rate-limited per client IP (in-process, shared by all sessions):
  after ``MAX_FAILURES`` failures within ``WINDOW_S`` seconds the IP is locked out for
  ``LOCKOUT_S`` seconds; every failure also costs a fixed delay.
* The password itself is never stored in session state; only a boolean "authenticated".

Limitation (documented in docs/web_app.md): over plain http the password travels in clear
text. For anything beyond short IP testing use an SSH tunnel or a TLS reverse proxy.
"""
from __future__ import annotations

import hashlib
import hmac
import os
import threading
import time

ENV_VAR = "NRW_APP_PASSWORD"
MIN_LENGTH = 8
MAX_FAILURES = 5
WINDOW_S = 300.0
LOCKOUT_S = 300.0
FAIL_DELAY_S = 1.0

_lock = threading.Lock()
_failures: dict[str, list[float]] = {}
_locked_until: dict[str, float] = {}


def configured_password() -> str | None:
    pw = os.environ.get(ENV_VAR, "")
    return pw if pw else None


def password_problem() -> str | None:
    """Return a human-readable configuration problem, or None if the setup is acceptable."""
    pw = configured_password()
    if pw is None:
        return (f"The environment variable {ENV_VAR} is not set. The administrator must start the "
                f"application with e.g.  export {ENV_VAR}='<strong password>'  (see docs/web_app.md).")
    if len(pw) < MIN_LENGTH:
        return f"{ENV_VAR} is shorter than {MIN_LENGTH} characters; choose a stronger password."
    return None


def _digest(x: str) -> bytes:
    return hashlib.sha256(x.encode("utf-8")).digest()


def check_password(candidate: str, expected: str | None = None) -> bool:
    exp = configured_password() if expected is None else expected
    if not exp or candidate is None:
        return False
    return hmac.compare_digest(_digest(candidate), _digest(exp))


def locked_out(client: str, now: float | None = None) -> float:
    """Seconds of lockout remaining for ``client`` (0 if not locked)."""
    now = time.time() if now is None else now
    with _lock:
        return max(0.0, _locked_until.get(client, 0.0) - now)


def register_attempt(client: str, ok: bool, now: float | None = None, sleep: bool = True) -> None:
    now = time.time() if now is None else now
    with _lock:
        if ok:
            _failures.pop(client, None)
            return
        fails = [t for t in _failures.get(client, []) if now - t < WINDOW_S] + [now]
        _failures[client] = fails
        if len(fails) >= MAX_FAILURES:
            _locked_until[client] = now + LOCKOUT_S
            _failures[client] = []
    if sleep:
        time.sleep(FAIL_DELAY_S)


def attempt_login(candidate: str, client: str, expected: str | None = None, sleep: bool = True):
    """Returns (ok, message)."""
    rem = locked_out(client)
    if rem > 0:
        return False, f"Too many failed attempts. Try again in {int(rem) + 1} s."
    ok = check_password(candidate, expected)
    register_attempt(client, ok, sleep=sleep)
    return (True, "") if ok else (False, "Incorrect password.")


def _reset_for_tests():
    with _lock:
        _failures.clear()
        _locked_until.clear()
