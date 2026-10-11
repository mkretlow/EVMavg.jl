# Contributing to EVMavg.jl

Thank you for considering contributing to EVMavg.jl! This document provides guidelines for contributing to the package.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:
- A clear description of the problem
- Steps to reproduce the bug
- Expected vs. actual behavior
- Your Julia version and OS
- A minimal code example demonstrating the issue

### Suggesting Enhancements

Enhancement suggestions are welcome! Please open an issue describing:
- The proposed feature or enhancement
- Why it would be useful
- How it might work

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following the coding standards below
3. **Add tests** for any new functionality
4. **Update documentation** as needed
5. **Ensure all tests pass** by running `julia --project -e 'using Pkg; Pkg.test()'`
6. **Submit a pull request** with a clear description of the changes

## Development Setup

```julia
# Clone the repository
git clone https://github.com/mkretlow/EVMavg.jl.git
cd EVMavg.jl

# Activate the project environment
julia --project

# Install dependencies
using Pkg
Pkg.instantiate()

# Run tests
Pkg.test()
```

## Coding Standards

### Style Guidelines

- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 92 characters (prefer 80 when possible)
- Use descriptive variable names
- Add docstrings to all public functions

### Documentation

All public functions must have docstrings following this format:

```julia
"""
    function_name(arg1, arg2; kwarg=default)

Brief description of what the function does.

# Arguments
- `arg1::Type`: Description of arg1
- `arg2::Type`: Description of arg2
- `kwarg::Type`: Description of keyword argument (default: `default`)

# Returns
- `ReturnType`: Description of return value

# Examples
```julia
result = function_name(value1, value2)
```

# Throws
- `ErrorType`: When this error occurs
"""
```

### Testing

- All new functionality must include tests
- Tests should cover:
  - Normal operation
  - Edge cases
  - Error conditions
  - Input validation
- Aim for high code coverage (>90%)

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in present tense (e.g., "Add", "Fix", "Update")
- Reference issue numbers when applicable

Good commit messages:
- ✓ "Add support for asymmetric uncertainties"
- ✓ "Fix confidence level calculation for small samples"
- ✓ "Update documentation for evm function"

Poor commit messages:
- ✗ "Fixed stuff"
- ✗ "Update"
- ✗ "WIP"

## Project Structure

```
EVMavg.jl/
├── src/
│   └── EVMavg.jl          # Main module file
├── test/
│   └── runtests.jl        # Test suite
├── examples/
│   └── basic_usage.jl     # Example scripts
├── Project.toml           # Package metadata
├── README.md              # Package overview
├── LICENSE                # License file
└── CONTRIBUTING.md        # This file
```

## Priority Areas for Contribution

We especially welcome contributions in these areas:

1. **Visualization Tools**: Functions to plot probability densities and weights
2. **Additional Examples**: Real-world applications from various fields
3. **Performance Optimization**: Improving computation speed for large datasets
4. **Documentation**: Expanding tutorials and use cases (e.g. a Documenter.jl site)

## Testing Your Changes

Before submitting a pull request:

```julia
# Run the test suite
using Pkg
Pkg.test("EVMavg")

# Check for performance regressions (if applicable)
using BenchmarkTools
@benchmark evm(data, errors)
```

## Code Review Process

1. All pull requests will be reviewed by maintainers
2. Reviewers may request changes or improvements
3. Once approved, your PR will be merged into `main`
4. Your contribution will be acknowledged in the release notes

## Questions?

If you have questions about contributing, feel free to:
- Open an issue on GitHub
- Contact the maintainer: Mike Kretlow <mike@kretlow.de>

## License

By contributing to EVMavg.jl, you agree that your contributions will be licensed under the MIT License.

Thank you for helping improve EVMavg.jl!
