"""Import shims for the two dependencies that do not import cleanly everywhere.

``antpack``
    Version 0.4 moved the numbering tools behind a license key and the package
    ``__init__`` imports its Qt GUI (``qt_material``), so a headless install raises
    ``ModuleNotFoundError`` before any numbering class is reachable. We pin
    ``antpack==0.3.8.6.2`` (pre-license); we stub the GUI module so either version
    imports, and surface a clear error if the installed version is license-gated.

the upstream numbering library
    The four modules the prescreen needs are vendored in :mod:`prescreen._proteintyper`,
    so the package is self-contained: nothing has to be on disk beside it. Setting
    ``$PROTEINTYPER_SRC`` to a checkout's ``src`` directory switches to the upstream
    library instead, which is how you check the two have not drifted. That path cannot
    simply import the package — its ``__init__`` pulls in dependencies the prescreen does
    not carry, and the submodules use relative imports — so it registers a synthetic
    parent pointing at the tree, letting the relative imports resolve without executing
    the real ``__init__``.
"""

from __future__ import annotations

import importlib
import os
import sys
import types

PROTEINTYPER_SRC_ENV = "PROTEINTYPER_SRC"
#: The vendored copy, used unless $PROTEINTYPER_SRC points somewhere else.
_VENDORED = f"{__package__}._proteintyper"


def _stub_qt_material() -> None:
    if "qt_material" in sys.modules:
        return
    try:
        importlib.import_module("qt_material")
    except ModuleNotFoundError:
        stub = types.ModuleType("qt_material")
        stub.apply_stylesheet = lambda *a, **k: None
        stub.list_themes = lambda *a, **k: []
        sys.modules["qt_material"] = stub


def _quiet_loguru() -> None:
    """The numbering modules log one DEBUG line per rejected sequence; silence it.

    Set ``PRESCREEN_VERBOSE=1`` to keep their logging.
    """
    if os.environ.get("PRESCREEN_VERBOSE"):
        return
    try:
        from loguru import logger
    except ModuleNotFoundError:
        return
    logger.disable("proteintyper_lib")
    logger.disable(_VENDORED)


def antpack_annotators():
    """Return ``(SingleChainAnnotator, PairedChainAnnotator)``.

    Raises
    ------
    RuntimeError
        If the installed antpack is license-gated (0.4+ without a key installed).
    """
    _stub_qt_material()
    try:
        from antpack import PairedChainAnnotator, SingleChainAnnotator
    except Exception as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "antpack is required for antibody-format classification. "
            "Install the pin used in production: pip install antpack==0.3.8.6.2"
        ) from exc
    return SingleChainAnnotator, PairedChainAnnotator


def proteintyper(module: str):
    """Import one of the numbering / paratope modules the prescreen depends on.

    The vendored copy in :mod:`prescreen._proteintyper` is used unless
    ``$PROTEINTYPER_SRC`` points at the upstream library checkout's ``src`` directory,
    which switches to that copy — the way to check the two have not drifted. An
    *installed* copy of the upstream library is deliberately not preferred: which code
    produced a verdict should depend on one explicit variable, not on what happens to be
    on the path.

    Args:
        module: ``cdr_novelty``, ``scaffold_cdr``, ``scaffold_search`` or
            ``novelty_scales``.
    """
    _stub_qt_material()
    _quiet_loguru()

    src = os.environ.get(PROTEINTYPER_SRC_ENV)
    if src:
        pkg_dir = os.path.join(src, "proteintyper_lib")
        if not os.path.isdir(pkg_dir):
            raise RuntimeError(
                f"${PROTEINTYPER_SRC_ENV} is set but {pkg_dir!r} does not exist; "
                "unset it to use the copy vendored in prescreen._proteintyper"
            )
        if "proteintyper_lib" not in sys.modules:
            pkg = types.ModuleType("proteintyper_lib")
            pkg.__path__ = [pkg_dir]
            sys.modules["proteintyper_lib"] = pkg
        return importlib.import_module(f"proteintyper_lib.{module}")

    return importlib.import_module(f"{_VENDORED}.{module}")
