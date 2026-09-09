PYTHON ?= python

.PHONY: check

check:
	$(PYTHON) -m compileall -q scripts tests
	$(PYTHON) tests/test_repo_tree.py
	$(PYTHON) tests/test_validation_workbook.py
	$(PYTHON) tests/test_temporal_outputs.py
	$(PYTHON) tests/test_corpus_manifests.py
