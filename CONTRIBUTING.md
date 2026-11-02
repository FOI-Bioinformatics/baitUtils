# Contributing to baitUtils

Thank you for your interest in contributing to baitUtils! This document provides guidelines for contributing to the project.

## 🚀 Quick Start for Contributors

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/baitUtils.git
   cd baitUtils
   ```
3. **Set up development environment**:
   ```bash
   conda create -n baitutils-dev -c conda-forge -c bioconda python=3.12 pblat minimap2 bedtools pybedtools viennarna
   conda activate baitutils-dev
   pip install -e ".[dev]"
   ```
4. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## 🧪 Development Guidelines

### Code Style
- Follow PEP 8 Python style guidelines
- Use type hints where appropriate
- Add docstrings to all functions and classes
- Keep functions focused and modular

### Testing
- Write tests for new behaviour; `tests/conftest.py` provides a small dataset
  with hand-derived expected values, fake `pblat` and `minimap2` scripts and a
  `run_cli` fixture for end-to-end runs
- Ensure all existing tests pass (the suite uses pytest fixtures):
  ```bash
  python -m pytest tests/
  ```
- Tests that need bedtools, pblat or ViennaRNA skip when the tool is absent;
  CI runs them in its `external-tools` job
- CI runs Python 3.12 and 3.13 on Ubuntu and macOS and reports coverage to
  Codecov; keep new code covered by tests

### Linting
CI fails on these checks:
```bash
ruff check baitUtils/ --select=E9,F63,F7,F82,F401,F811,F841
flake8 baitUtils/ --select=E9,F63,F7,F82
```

### Documentation
- Update docstrings and help text
- Add examples for new commands or features
- Update README.md if adding new functionality or changing options and outputs
- Update CLAUDE.md with technical details
- Add an entry to the "Unreleased" section of CHANGELOG.md

## 📝 Contribution Types

### Bug Reports
- Open a GitHub issue
- Include detailed reproduction steps
- Provide system information and versions
- Include sample data if possible

### Feature Requests
- Describe the use case and motivation
- Provide examples of desired functionality
- Consider backward compatibility

### Code Contributions
- Follow the development workflow above
- Include comprehensive tests
- Update documentation
- Ensure all CI checks pass

## 🔧 Architecture Overview

Each command (`stats`, `plot`, `map`, `check`, `fill`, `evaluate`,
`compare`) has a module with `add_arguments(parser)` and `main(args)`;
`__main__.py` dispatches to them. CLAUDE.md maps every command and module.

### Key Components
- `mapping_utils.py` - pblat and minimap2 runners and the single PSL/PAF parser
- `coverage_stats.py` - Per-base coverage, breadth, depth and exports
- `gap_analysis.py`, `reference_analyzer.py` - Gap and reference features
- `quality_scorer.py` - Quality score in 0-1
- `design.py` - Colours, fonts and chart style shared by every figure and report
- `coverage_viz.py`, `interactive_plots.py`, `comparative_visualizations.py` - Figures
- `report_generator.py`, `comparative_report_generator.py` - HTML reports

### Testing Structure
- Unit tests in `tests/` directory, run with pytest
- End-to-end command runs on the fixture dataset (`tests/conftest.py`)
- Regression tests per audit phase (`tests/test_phase*.py`)

## 🎯 Contribution Areas

We welcome contributions in:
- **New statistical methods** for oligo evaluation
- **Additional visualization types** for coverage analysis
- **Performance optimizations** for large datasets
- **New export formats** and reporting options
- **Documentation improvements** and examples
- **Bug fixes** and code quality improvements

## 📋 Pull Request Process

1. **Update documentation** for any user-facing changes
2. **Add or update tests** to maintain coverage
3. **Run the full test suite** locally
4. **Leave the version alone**; maintainers set it in `_version.py` when releasing
5. **Write clear commit messages**
6. **Submit pull request** with detailed description

### PR Template
```
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature  
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
- [ ] Added/updated unit tests
- [ ] All tests pass locally
- [ ] Tested with sample data

## Documentation
- [ ] Updated docstrings
- [ ] Updated README if needed
- [ ] Updated help text

## Checklist
- [ ] Code follows project style guidelines
- [ ] Self-review completed
- [ ] Backward compatibility maintained
```

## ❓ Questions?

- Open an issue for questions about contributing
- Email: andreas.sjodin@gmail.com

Thank you for contributing to baitUtils! 🧬