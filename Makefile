.PHONY: install test lint gate release-validate clean

install:
	pip install -e ".[dev]"

test:
	pytest -q

lint:
	ruff check src tests scripts

gate:
	@if python scripts/validate_release.py data/releases/smoke; then \
		echo "FAIL: validator accepted planted errors"; exit 1; \
	else \
		echo "OK: validator rejected planted errors"; \
	fi

release-validate:
	python scripts/validate_release.py data/releases/$(REL)

clean:
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache .ruff_cache *.egg-info build dist
