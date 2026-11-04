# banger

A fast, simple math calculator, symbolic engine, and terminal grapher written in Go. Zero external dependencies.

## run

Interactive mode:

```sh
go run .
```

Single expression from the command line:

```sh
go run . "int(e^(-x^2), x, -inf, inf)"
```

Run tests:

```sh
go test ./...
```

## syntax and session

- Implicit multiplication: `2x`, `3pi`, `2(x + 1)`, `(x + 1)(x - 1)`, `2sin(x)`
- Right-associative powers: `2^3^2`, `2^-3`, `**`
- Postfix factorials: `5!`, `7!!`
- Base literals: `0xFF`, `0b1010`, `0o77`
- Variables and custom functions: `a = 5`, `f(x) = x^3 - 3x`, `g(x, y) = sqrt(x^2 + y^2)`
- Prime derivative notation: `f'(x)`, `f''(2)`, `sin'(0)`
- History and last answer: `ans`, `_`, `$1`, `$2`
- Statement chaining on one line: `a = 3; b = 4; hypot(a, b)`
- Automatic exact value hints when applicable (`[= pi^2/6]`, `[= sqrt(pi)]`, `[= 1/3]`)
- Built-in constants: `pi`, `tau`, `e`, `phi`, `i` (or `j`), `euler`, `catalan`, `inf`, `c_light`, `g_grav`, `planck`, `hbar`, `boltz`, `avogadro`

## commands

- `help` or `help <topic>` (`calc`, `linalg`, `complex`, `num`, `stats`, `special`, `plot`)
- `mode rad` or `mode deg` to switch angle units
- `vars` and `funcs` to list user variables and functions
- `history` or `hist` to view previous results
- `reset` to clear session state
- `clear` or `cls` to clear the terminal
- `exit` or `quit` to leave

## graphing and tables

High-resolution Unicode Braille plotting in the terminal:

- `plot <f1 [, f2, ...]> [xmin xmax [ymin ymax]]`
- `plotparam <x(t)>, <y(t)> [tmin tmax]`
- `plotpolar <r(t)> [tmin tmax]`
- `slopefield <dy/dx> [xmin xmax ymin ymax]`
- `table <f1 [, f2, ...]> [xmin xmax step]`

## math features

### algebra and complex numbers

Full complex support across standard operations and functions (`3 + 4i`, `e^(i*pi) + 1`, `sqrt(-4)`).

- Core: `sqrt`, `cbrt`, `root(n, x)`, `pow(x, y)`, `exp`, `ln`, `log` / `log10`, `log2`, `logb(b, x)`
- Rounding and sign: `abs`, `floor`, `ceil`, `round(x [, places])`, `trunc`, `mod(a, b)`, `sgn`, `heaviside`, `hypot`
- Angles: `deg2rad`, `rad2deg`
- Complex: `re`, `im`, `conj`, `arg` / `phase`, `polar`, `rect(r, theta)`
- Trig and hyperbolic: `sin`, `cos`, `tan`, `sec`, `csc`, `cot`, `asin`, `acos`, `atan`, `atan2`, `asec`, `acsc`, `acot`, `sinh`, `cosh`, `tanh`, `sech`, `csch`, `coth`, `asinh`, `acosh`, `atanh`

### calculus and differential equations

- Symbolic derivatives and integrals: `diff(expr [, x [, order [, pt]]])`, `int(expr [, x])` / `antideriv`
- Symbolic algebra and series: `simplify(expr)`, `expand(expr [, x])`, `at(expr, x, val)`, `taylor(expr [, x [, center [, order]]])`, `maclaurin(expr [, x [, order]])`
- Definite and improper integrals: `int(expr, x, a, b)` using adaptive 15-point Gauss-Kronrod quadrature (supports `-inf`, `inf`, and complex integrands)
- Numerical derivatives and limits: `ndiff(expr, x, pt [, order])`, `lim(expr, x, target [, dir])`
- Sums and products: `sum(expr, k, a, b)` (supports `inf`), `prod(expr, k, a, b)`
- Equation solving: `solve(eq, x [, min, max])` / `roots` (finds all real and complex polynomial roots, or real roots of transcendental equations)
- Differential equations: `rk4(dydx, x, y, x0, y0, xEnd [, steps])`
- Vector calculus: `grad`, `hessian`, `div` / `divergence`, `curl`, `laplacian`

### linear algebra and vectors

Vectors use `[1, 2, 3]` and matrices use `[1, 2; 3, 4]` or `[[1, 2], [3, 4]]` (real or complex entries).

- Operators: `A + B`, `A - B`, `A * B`, `A * v`, `A^k`, `A^-1`, `A'` (conjugate transpose)
- Vectors: `dot`, `cross`, `outer`, `norm(x [, p])`, `unit`, `proj`, `angle`, `dim` / `size` / `len`
- Builders: `eye`, `zeros`, `ones`, `diag`, `hilbert`, `vandermonde`
- Reductions and subspaces: `rref`, `ref`, `rank`, `nullity`, `nullspace` / `kernel`, `colspace`, `rowspace`, `gramschmidt`
- Solvers and decompositions: `det`, `inv`, `pinv`, `adj`, `trace` / `tr`, `trans`, `ctrans`, `linsolve(A, b)`, `eig` / `eigenvals`, `eigvecs`, `charpoly`, `expm`, `cond`

### number theory, combinatorics, and bases

- Primes and divisors: `isprime`, `nextprime`, `prevprime`, `nthprime`, `primepi`, `factor`, `primefactors`, `divisors`, `sigma(n [, k])`, `totient` / `eulerphi`, `mobius`
- Modular arithmetic: `gcd`, `lcm`, `xgcd`, `modpow(b, e, m)`, `modinv(a, m)`, `crt([rems], [mods])`, `znorder(a, n)`, `primroot(n)`, `legendre(a, p)`, `jacobi(a, n)`
- Combinatorics: `fact`, `ncr` / `binom` / `choose`, `npr` / `perm`, `fib`, `lucas`, `catalannum`, `stirling1`, `stirling2`, `bell`, `partition`, `subfact` / `derange`
- Rationals and bases: `frac(x [, maxDen])`, `cf(x [, terms])`, `hex`, `bin`, `oct`, `base(n, b)`

### special functions

- `gamma`, `lgamma`, `digamma` / `psi`, `beta`, `gammainc`, `betainc`
- `zeta(s)` across the complex plane (including `s <= 0`)
- `lambertw(x [, branch])` for `0` and `-1` branches
- `erf`, `erfc`, `erfinv`, `besselj(n, x)`, `bessely(n, x)`, `sinc`

### statistics and probability

Descriptive functions accept a vector `mean([1, 2, 3])` or multiple arguments `mean(1, 2, 3)`.

- Descriptive: `mean` / `avg`, `geomean`, `harmmean`, `rms`, `median`, `mode`, `var`, `pvar`, `std` / `stdev`, `pstd`, `iqr`, `mad`, `skewness`, `kurtosis`, `quantile(X, p)`, `percentile(X, p)`, `min`, `max`, `sum`, `prod`, `sort`, `zscore`
- Two-variable and regression: `cov(X, Y)`, `corr(X, Y)`, `linreg(X, Y)` (returns `[slope, intercept, r^2]`)
- Distributions: `normpdf`, `normcdf`, `norminv`, `tpdf`, `tcdf`, `tinv`, `chipdf`, `chicdf`, `chiinv`, `fpdf`, `fcdf`, `binompmf`, `binomcdf`, `poisspmf`, `poisscdf`, `geompmf`, `geomcdf`, `exppdf`, `expcdf`

## examples

```text
Calculate > diff(x^x, x)
Result    = x^x*(ln(x) + 1)

Calculate > int(e^(-x^2), x, -inf, inf)
Result    = 1.77245385091  [= sqrt(pi)]

Calculate > taylor(sin(x), x, 0, 7)
Result    = x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7

Calculate > solve(x^3 - 6x^2 + 11x - 6 = 0, x)
Result    = [1, 2, 3]

Calculate > eig([0, -1; 1, 0])
Result    = [-i, i]

Calculate > crt([2, 3, 2], [3, 5, 7])
Result    = 23

Calculate > plot sin(x), cos(x) -2*pi 2*pi
Calculate > plotpolar 1 + cos(3t) 0 2*pi
Calculate > slopefield x - y -4 4 -4 4
```
