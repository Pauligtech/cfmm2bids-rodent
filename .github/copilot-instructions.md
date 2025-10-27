# GitHub Copilot Instructions for cfmm2bids

## Code Quality Tools

This repository uses several tools to maintain code quality and consistency:

### 1. Ruff (Python Linting and Formatting)

Ruff is a fast Python linter and formatter that combines the functionality of multiple tools (flake8, black, isort, etc.).

**Configuration**: `pyproject.toml`

**Usage**:
```bash
# Check for linting issues
ruff check .

# Automatically fix issues
ruff check --fix .

# Check formatting
ruff format --check .

# Apply formatting
ruff format .
```

**In your code**: 
- Follow PEP 8 style guidelines
- Use double quotes for strings
- Line length is 88 characters
- Import statements should be sorted and organized

### 2. Snakefmt (Snakemake Formatting)

Snakefmt formats Snakemake files (Snakefile, *.smk) for consistency.

**Usage**:
```bash
# Check formatting
snakefmt --check --compact-diff .

# Apply formatting
snakefmt .
```

**In your code**:
- Snakemake files should follow snakefmt conventions
- Use consistent indentation and spacing

### 3. Pre-commit Hooks

Pre-commit hooks automatically run checks before commits.

**Setup**:
```bash
# Install pre-commit (if not already installed)
pip install pre-commit

# Install the git hooks
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

**Hooks configured**:
- `ruff` - Python linting with auto-fix
- `ruff-format` - Python formatting
- `snakefmt` - Snakemake formatting
- `trailing-whitespace` - Remove trailing whitespace
- `end-of-file-fixer` - Ensure files end with newline
- `check-yaml` - Validate YAML syntax
- `check-added-large-files` - Prevent committing large files
- `check-merge-conflict` - Check for merge conflict markers

### 4. GitHub Actions CI

Two GitHub Actions workflows run automatically on push and pull requests:

**Lint Workflow** (`.github/workflows/lint.yml`):
- Runs `ruff check` and `ruff format --check`
- Ensures Python code meets linting and formatting standards

**Snakeformat Workflow** (`.github/workflows/snakeformat.yml`):
- Runs `snakefmt --check`
- Ensures Snakemake files are properly formatted

## Development Workflow

1. **Before committing**:
   ```bash
   # Run all pre-commit hooks
   pre-commit run --all-files
   ```

2. **Fix issues automatically**:
   ```bash
   # Fix Python linting issues
   ruff check --fix .
   
   # Format Python code
   ruff format .
   
   # Format Snakemake files
   snakefmt .
   ```

3. **Commit your changes**:
   - Pre-commit hooks will run automatically
   - If hooks fail, fix the issues and commit again

4. **Push to GitHub**:
   - GitHub Actions will run the lint and snakeformat checks
   - PR will show check status

## Copilot Coding Guidelines

When generating or modifying code:

### Python Files
- Follow PEP 8 conventions
- Use type hints where appropriate
- Write docstrings for functions and classes
- Keep functions focused and modular
- Use meaningful variable names
- Imports should be organized: standard library, third-party, local

### Snakemake Files
- Follow snakefmt formatting
- Use descriptive rule names
- Include input, output, log, and resources where appropriate
- Add comments to explain complex logic
- Keep rules focused on single tasks

### General
- Write clear commit messages
- Update documentation when changing functionality
- Test your changes before committing
- Run formatters and linters before pushing

## Quick Reference

```bash
# Check everything
ruff check . && ruff format --check . && snakefmt --check .

# Fix everything
ruff check --fix . && ruff format . && snakefmt .

# Pre-commit on all files
pre-commit run --all-files

# Update pre-commit hooks
pre-commit autoupdate
```
