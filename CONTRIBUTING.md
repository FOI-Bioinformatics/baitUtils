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
   conda create -n baitutils-dev python=3.12
   conda activate baitutils-dev
   conda install -c bioconda pblat bedtools
   pip install -e .
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
- Write unit tests for new functionality
- Ensure all existing tests pass:
  ```bash
  python -m unittest discover tests/
  ```
- Test coverage should be maintained above 80%

### Documentation
- Update docstrings and help text
- Add examples for new commands or features
- Update README.md if adding new functionality
- Update CLAUDE.md with technical details

## 📝 Contribution Types

### Bug Reports
- Use the GitHub issue template
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

### Core Modules
- **Phase 1**: `evaluate.py` - Individual oligo set analysis
- **Phase 2**: Enhanced evaluation with interactive reporting  
- **Phase 3**: `compare.py` - Multi-set comparative analysis

### Key Components
- `coverage_stats.py` - Coverage analysis engine
- `quality_scorer.py` - Quality assessment system
- `*_visualizations.py` - Plotting and visualization
- `*_report_generator.py` - HTML report generation

### Testing Structure
- Unit tests in `tests/` directory
- Mock-based testing for file I/O
- Integration tests for command workflows

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
4. **Update version** in `_version.py` if needed
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
- Join discussions in GitHub Discussions
- Email: andreas.sjodin@gmail.com

Thank you for contributing to baitUtils! 🧬