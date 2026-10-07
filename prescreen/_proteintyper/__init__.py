"""The numbering and paratope modules the prescreen shares with its upstream library.

The prescreen does not classify antibodies itself: the IMGT numbering, the CDR identity
metric and the scaffold paratope masks all come from the upstream library this package was
split out of. Sharing that code is the point — a prescreen that numbered sequences
differently from the novelty scales it sits beside would give a different answer about the
same design.

Four modules are copied verbatim from upstream (commit ``a02fae8``):

``novelty_scales``   the published cut points (only ``CDR_DE_NOVO_CUTOFF`` is read here)
``cdr_novelty``      antpack IMGT numbering, ``extract_cdrs``, ``cdr_identity``,
                     ``compare_regions``
``scaffold_cdr``     affibody / monobody / DARPin references, masks, framework identity
``scaffold_search``  MMseqs2 search of the scaffold families

Only one thing is changed, and it is marked ``vendored:``: the two branches that recover a
paratope from a *structure* still import the upstream ``paratope`` module and now raise a
clear error when it is absent. Both are reachable only when a caller passes
``structure_path``; the prescreen is sequence-only and never does (see the "no structure
arm" blind spot in ``PRESCREEN_SPEC.md``).

The copy keeps the prescreen deployable on its own — before it, importing the package
anywhere without the upstream source on disk failed at the first ``classify()`` call, which
is a long way from the import that caused it. To run against the upstream copy instead (to
check the two have not drifted), set ``$PROTEINTYPER_SRC`` to its ``src`` directory;
:func:`prescreen.compat.proteintyper` prefers it when it is set.

Refreshing the copy is a file copy plus re-applying the two marked guards; the modules
import nothing else from upstream.
"""
