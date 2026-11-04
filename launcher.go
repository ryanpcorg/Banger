package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strings"
)

func printBanner() {
	println("██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗██████╗ ")
	println("██╔══██╗██╔══██╗████╗  ██║██╔════╝ ██╔════╝██╔══██╗")
	println("██████╔╝███████║██╔██╗ ██║██║  ███╗█████╗  ██████╔╝")
	println("██╔══██╗██╔══██║██║╚██╗██║██║   ██║██╔══╝  ██╔══██╗")
	println("██████╔╝██║  ██║██║ ╚████║╚██████╔╝███████╗██║  ██║")
	println("╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝")
	println("---------------------------------------------------")
	println("The calculator your mum hates and you dad thinks is BS (Ballistically Stupid)")
	println("Type 'help' for CAS/Calculus/LinAlg/Stats/Graphing docs | 'clear' | 'exit'\n")
}

func splitTopLevel(s string, sep byte) []string {
	var parts []string
	depthParen := 0
	depthBrack := 0
	start := 0
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case '(':
			depthParen++
		case ')':
			if depthParen > 0 {
				depthParen--
			}
		case '[':
			depthBrack++
		case ']':
			if depthBrack > 0 {
				depthBrack--
			}
		default:
			if s[i] == sep && depthParen == 0 && depthBrack == 0 {
				parts = append(parts, strings.TrimSpace(s[start:i]))
				start = i + 1
			}
		}
	}
	parts = append(parts, strings.TrimSpace(s[start:]))
	return parts
}

func splitTopLevelWords(s string) []string {
	var words []string
	depthParen := 0
	depthBrack := 0
	start := -1
	for i := 0; i < len(s); i++ {
		ch := s[i]
		switch ch {
		case '(':
			if start < 0 {
				start = i
			}
			depthParen++
		case ')':
			if depthParen > 0 {
				depthParen--
			}
		case '[':
			if start < 0 {
				start = i
			}
			depthBrack++
		case ']':
			if depthBrack > 0 {
				depthBrack--
			}
		case ' ', '\t':
			if depthParen == 0 && depthBrack == 0 {
				if start >= 0 {
					words = append(words, s[start:i])
					start = -1
				}
			}
		default:
			if start < 0 {
				start = i
			}
		}
	}
	if start >= 0 {
		words = append(words, s[start:])
	}
	return words
}

func endsWithOperator(tok string) bool {
	if tok == "" {
		return true
	}
	last := tok[len(tok)-1]
	return last == '+' || last == '-' || last == '*' || last == '/' || last == '^' || last == ',' || last == '='
}

func tryEvalBound(eng *Engine, s string) (float64, bool) {
	nd, err := parseMath(s)
	if err != nil {
		return 0, false
	}
	if dependsOn(nd, "x") || dependsOn(nd, "X") || dependsOn(nd, "t") || dependsOn(nd, "T") || dependsOn(nd, "y") || dependsOn(nd, "theta") {
		return 0, false
	}
	v, err := eng.evalAST(nd, nil, 0)
	if err != nil || !v.IsReal() {
		return 0, false
	}
	r := real(v.S)
	if math.IsNaN(r) || math.IsInf(r, 0) {
		return 0, false
	}
	return r, true
}

func parseTrailingBounds(eng *Engine, body string, count int) (string, []float64, bool) {
	words := splitTopLevelWords(body)
	if len(words) <= count {
		return body, nil, false
	}
	prevTok := words[len(words)-count-1]
	if endsWithOperator(prevTok) {
		return body, nil, false
	}
	vals := make([]float64, count)
	for i := 0; i < count; i++ {
		v, ok := tryEvalBound(eng, words[len(words)-count+i])
		if !ok {
			return body, nil, false
		}
		vals[i] = v
	}
	exprPart := strings.Join(words[:len(words)-count], " ")
	return strings.TrimSpace(exprPart), vals, true
}

func printFormattedResult(val Value) {
	s := val.String()
	lines := strings.Split(s, "\n")
	if len(lines) == 1 {
		fmt.Printf("Result    = %s\n\n", lines[0])
		return
	}
	fmt.Printf("Result    = %s\n", lines[0])
	indent := strings.Repeat(" ", 12)
	for _, ln := range lines[1:] {
		fmt.Printf("%s%s\n", indent, ln)
	}
	fmt.Println()
}

func printHelp(topic string) {
	switch strings.ToLower(strings.TrimSpace(topic)) {
	case "calc", "calculus":
		fmt.Println("CALCULUS & DIFFERENTIAL EQUATIONS:")
		fmt.Println("  diff(expr [, x [, order [, pt]]])  Symbolic derivative: diff(x^x, x), diff(sin(x), x, 3)")
		fmt.Println("  f'(x), f''(2), sin'(0)             Lagrange prime derivative notation")
		fmt.Println("  ndiff(expr, x, pt [, order])       High-order Richardson numerical derivative")
		fmt.Println("  int(expr [, x])                    Symbolic indefinite integral: int(3x^2 + cos(2x), x)")
		fmt.Println("  int(expr, x, a, b)                 Adaptive 15-pt Gauss-Kronrod integral (supports -inf, inf)")
		fmt.Println("  taylor(expr, x, center, order)     Symbolic Taylor series: taylor(sin(x), x, 0, 7)")
		fmt.Println("  maclaurin(expr [, x [, order]])    Maclaurin polynomial around 0")
		fmt.Println("  lim(expr, x, target [, dir])       Limits (+1 right, -1 left, inf): lim(sin(x)/x, x, 0)")
		fmt.Println("  sum(expr, k, a, b)                 Finite & infinite series: sum(1/n^2, n, 1, inf)")
		fmt.Println("  prod(expr, k, a, b)                Finite products: prod(1 - 1/k^2, k, 2, 10)")
		fmt.Println("  solve(eq, x [, min, max])          All complex poly roots or real transcendental roots")
		fmt.Println("  rk4(dydx, x, y, x0, y0, x1 [, n])  4th-order Runge-Kutta ODE solver")
		fmt.Println("  grad, hessian, div, curl, laplacian  Multivariable vector calculus")
		fmt.Println("  simplify(expr), expand(expr, x), at(expr, x, val)")
		fmt.Println()

	case "linalg", "matrix", "vectors":
		fmt.Println("LINEAR ALGEBRA & VECTORS:")
		fmt.Println("  Vectors:  u = [1, 2, 3]            Matrices: A = [1, 2; 3, 4] or [[1, 2], [3, 4]]")
		fmt.Println("  Operators: A + B, A - B, A * B, A * v, A^3, A^-1 (inverse), A' (conjugate transpose)")
		fmt.Println("  dot(u, v), cross(u, v), outer(u, v), norm(v [, p]), unit(v), proj(u, v), angle(u, v)")
		fmt.Println("  det(A), inv(A), pinv(A), adj(A), trace(A), trans(A), ctrans(A), cond(A), expm(A)")
		fmt.Println("  rref(A), ref(A), rank(A), nullity(A), nullspace(A), colspace(A), rowspace(A)")
		fmt.Println("  linsolve(A, b)                     Exact LU or least-squares Ax = b solver")
		fmt.Println("  eig(A), eigvecs(A), charpoly(A)    Complex/real eigenvalues, eigenvectors, char polynomial")
		fmt.Println("  gramschmidt(A), eye(n), zeros(r,c), ones(r,c), diag(v), hilbert(n), vandermonde(v)")
		fmt.Println()

	case "complex":
		fmt.Println("COMPLEX NUMBERS:")
		fmt.Println("  Imaginary unit: i (or literals like 3i, 2.5i). Euler identity: e^(i*pi) + 1")
		fmt.Println("  re(z), im(z), abs(z), arg(z), conj(z), polar(z), rect(r, theta)")
		fmt.Println("  All elementary & special functions (sqrt(-4), ln(-1), sin(1+i), zeta, gamma) support C.")
		fmt.Println()

	case "num", "numtheory":
		fmt.Println("NUMBER THEORY, COMBINATORICS & BASES:")
		fmt.Println("  gcd(...), lcm(...), xgcd(a, b)     Extended GCD [g, x, y] with ax + by = g")
		fmt.Println("  isprime(n), nextprime(n), prevprime(n), nthprime(k), primepi(x)")
		fmt.Println("  factor(n), primefactors(n), divisors(n), sigma(n [, k]), totient(n), mobius(n)")
		fmt.Println("  modpow(b, e, m), modinv(a, m), crt([rems], [mods]), znorder(a, n), primroot(n)")
		fmt.Println("  legendre(a, p), jacobi(a, n)")
		fmt.Println("  n!, n!!, ncr(n, r), npr(n, r), fib(n), lucas(n), catalannum(n)")
		fmt.Println("  stirling1(n, k), stirling2(n, k), bell(n), partition(n), subfact(n)")
		fmt.Println("  frac(x), cf(x), hex(n), bin(n), oct(n), base(n, b), literals 0xFF, 0b1010, 0o77")
		fmt.Println()

	case "stats", "prob":
		fmt.Println("STATISTICS & PROBABILITY DISTRIBUTIONS:")
		fmt.Println("  Descriptive (vector or variadic): mean, geomean, harmmean, rms, median, mode")
		fmt.Println("  var, pvar, std, pstd, iqr, mad, skewness, kurtosis, quantile(X, p), percentile(X, p)")
		fmt.Println("  min, max, sum, prod, sort(X), zscore(X), cov(X, Y), corr(X, Y)")
		fmt.Println("  linreg(X, Y)                       Returns [slope, intercept, R^2]")
		fmt.Println("  Normal:     normpdf(x [,mu,sd]), normcdf(x [,mu,sd]), norminv(p [,mu,sd])")
		fmt.Println("  Student-t:  tpdf(t, df), tcdf(t, df), tinv(p, df)")
		fmt.Println("  Chi-Square: chipdf(x, df), chicdf(x, df), chiinv(p, df)")
		fmt.Println("  F-dist:     fpdf(x, d1, d2), fcdf(x, d1, d2)")
		fmt.Println("  Discrete:   binompmf/cdf(k,n,p), poisspmf/cdf(k,lam), geompmf/cdf(k,p), exppdf/cdf(x,lam)")
		fmt.Println()

	case "special":
		fmt.Println("SPECIAL FUNCTIONS & CONSTANTS:")
		fmt.Println("  gamma(z), lgamma(x), digamma(x), beta(a, b), gammainc(a, x), betainc(a, b, x)")
		fmt.Println("  zeta(s)                            Riemann zeta function across C (incl. s < 0)")
		fmt.Println("  lambertw(x [, branch])             Product log W_0(x) and W_-1(x)")
		fmt.Println("  erf(x), erfc(x), erfinv(x), besselj(n, x), bessely(n, x), sinc(z)")
		fmt.Println("  Constants: pi, tau, e, phi, i, euler, catalan, inf, c_light, g_grav, planck, hbar")
		fmt.Println()

	case "plot", "graph":
		fmt.Println("HIGH-RESOLUTION BRAILLE GRAPHING & TABLES:")
		fmt.Println("  plot <f1 [, f2, ...]> [xmin xmax [ymin ymax]]")
		fmt.Println("    Examples: plot sin(x), cos(x) -2*pi 2*pi")
		fmt.Println("              plot tan(x) -5 5 -6 6")
		fmt.Println("  plotparam <x(t)>, <y(t)> [tmin tmax]")
		fmt.Println("    Example:  plotparam cos(3t), sin(2t) 0 2*pi")
		fmt.Println("  plotpolar <r(t)> [tmin tmax]")
		fmt.Println("    Example:  plotpolar 1 + cos(t) 0 2*pi")
		fmt.Println("  slopefield <dy/dx> [xmin xmax ymin ymax]")
		fmt.Println("    Example:  slopefield x - y -4 4 -4 4")
		fmt.Println("  table <f1 [, f2]> [xmin xmax step]")
		fmt.Println("    Example:  table x^2, sin(x) -2 2 0.5")
		fmt.Println()

	default:
		fmt.Println("BANGER MATHEMATICAL ENGINE - HELP OVERVIEW")
		fmt.Println("  Syntax:     Implicit multiplication (2x, 3pi, 2(x+1), sin(2x)), right-assoc ^, n!, n!!")
		fmt.Println("  Variables:  a = 42, A = [1, 2; 3, 4], ans (last result), $1, $2 (history refs)")
		fmt.Println("  Functions:  f(x) = x^3 - 3x, g(x, y) = sqrt(x^2 + y^2), derivatives f'(x), f''(2)")
		fmt.Println("  Chaining:   a = 3; b = 4; hypot(a, b)")
		fmt.Println()
		fmt.Println("  Session Commands:")
		fmt.Println("    mode [rad|deg]   Switch or view angle mode (default: rad)")
		fmt.Println("    vars / funcs     Inspect defined variables and custom functions")
		fmt.Println("    history          Show numbered calculation history ($1, $2, ...)")
		fmt.Println("    reset            Clear user variables, functions, and history")
		fmt.Println("    clear / cls      Clear terminal screen")
		fmt.Println("    exit / quit      Exit Banger")
		fmt.Println()
		fmt.Println("  Detailed Topic Help:")
		fmt.Println("    help calc        Symbolic/numerical derivatives, integrals, Taylor, limits, ODEs, solve")
		fmt.Println("    help linalg      Vectors, matrices, RREF, eigenvalues, nullspace, linear systems")
		fmt.Println("    help complex     Complex arithmetic, polar/rectangular forms, complex analysis")
		fmt.Println("    help num         Number theory, primes, CRT, modular arithmetic, combinatorics")
		fmt.Println("    help stats       Descriptive statistics, regression, probability distributions")
		fmt.Println("    help special     Zeta, Gamma, Beta, Lambert W, Bessel, Error functions")
		fmt.Println("    help plot        Unicode Braille multi-curve, parametric, polar & slope-field plots")
		fmt.Println()
	}
}

func handleCommand(eng *Engine, input string) bool {
	lower := strings.ToLower(input)

	if lower == "clear" || lower == "cls" {
		clearScreen()
		printBanner()
		return true
	}

	if lower == "help" || strings.HasPrefix(lower, "help ") {
		topic := ""
		if len(input) > 5 {
			topic = input[5:]
		}
		printHelp(topic)
		return true
	}

	if lower == "reset" {
		eng.Reset()
		fmt.Print("Session variables, functions, and history cleared.\n\n")
		return true
	}

	if lower == "mode" || strings.HasPrefix(lower, "mode ") {
		arg := strings.TrimSpace(lower[4:])
		switch arg {
		case "":
			if eng.Mode == ModeRad {
				fmt.Print("Angle mode: rad (radians)\n\n")
			} else {
				fmt.Print("Angle mode: deg (degrees)\n\n")
			}
		case "rad", "radians":
			eng.Mode = ModeRad
			fmt.Print("Angle mode set to radians.\n\n")
		case "deg", "degrees":
			eng.Mode = ModeDeg
			fmt.Print("Angle mode set to degrees.\n\n")
		default:
			fmt.Print("Usage: mode [rad|deg]\n\n")
		}
		return true
	}

	if lower == "vars" {
		names := eng.SortedVarNames()
		if len(names) == 0 {
			fmt.Print("No user variables defined. (Built-in: pi, tau, e, phi, i, euler, catalan, ans)\n\n")
			return true
		}
		for _, k := range names {
			fmt.Printf("  %s = %s\n", k, eng.Vars[k].String())
		}
		fmt.Println()
		return true
	}

	if lower == "funcs" {
		names := eng.SortedFuncNames()
		if len(names) == 0 {
			fmt.Print("No user functions defined. Define one with e.g. f(x) = x^3 - 2x\n\n")
			return true
		}
		for _, k := range names {
			fmt.Printf("  %s\n", eng.Funcs[k].Raw)
		}
		fmt.Println()
		return true
	}

	if lower == "history" || lower == "hist" {
		if len(eng.History) == 0 {
			fmt.Print("History is empty.\n\n")
			return true
		}
		for idx, h := range eng.History {
			fmt.Printf("  $%d: %s  =>  %s\n", idx+1, h.Input, strings.ReplaceAll(h.Val.String(), "\n", " "))
		}
		fmt.Println()
		return true
	}

	if strings.HasPrefix(lower, "plot ") || strings.HasPrefix(lower, "graph ") {
		idx := strings.IndexByte(input, ' ')
		body := strings.TrimSpace(input[idx+1:])
		xMin, xMax := -10.0, 10.0
		yMin, yMax := math.NaN(), math.NaN()

		if rem, b4, ok := parseTrailingBounds(eng, body, 4); ok {
			body = rem
			xMin, xMax, yMin, yMax = b4[0], b4[1], b4[2], b4[3]
		} else if rem, b2, ok := parseTrailingBounds(eng, body, 2); ok {
			body = rem
			xMin, xMax = b2[0], b2[1]
		}

		exprs := splitTopLevel(body, ',')
		out, err := eng.PlotFunctions(exprs, xMin, xMax, yMin, yMax, 64, 18, true)
		if err != nil {
			fmt.Printf("Plot error: %v\n\n", err)
		} else {
			fmt.Println(out)
		}
		return true
	}

	if strings.HasPrefix(lower, "plotparam ") {
		body := strings.TrimSpace(input[len("plotparam "):])
		tMin, tMax := 0.0, 2*math.Pi
		if rem, b2, ok := parseTrailingBounds(eng, body, 2); ok {
			body = rem
			tMin, tMax = b2[0], b2[1]
		}
		parts := splitTopLevel(body, ',')
		if len(parts) != 2 {
			fmt.Print("Usage: plotparam <x(t)>, <y(t)> [tmin tmax]\n\n")
			return true
		}
		out, err := eng.PlotParametric(parts[0], parts[1], tMin, tMax, 64, 18, true)
		if err != nil {
			fmt.Printf("Plot error: %v\n\n", err)
		} else {
			fmt.Println(out)
		}
		return true
	}

	if strings.HasPrefix(lower, "plotpolar ") {
		body := strings.TrimSpace(input[len("plotpolar "):])
		tMin, tMax := 0.0, 2*math.Pi
		if rem, b2, ok := parseTrailingBounds(eng, body, 2); ok {
			body = rem
			tMin, tMax = b2[0], b2[1]
		}
		out, err := eng.PlotPolar(body, tMin, tMax, 64, 18, true)
		if err != nil {
			fmt.Printf("Plot error: %v\n\n", err)
		} else {
			fmt.Println(out)
		}
		return true
	}

	if strings.HasPrefix(lower, "slopefield ") {
		body := strings.TrimSpace(input[len("slopefield "):])
		xMin, xMax, yMin, yMax := -5.0, 5.0, -5.0, 5.0
		if rem, b4, ok := parseTrailingBounds(eng, body, 4); ok {
			body = rem
			xMin, xMax, yMin, yMax = b4[0], b4[1], b4[2], b4[3]
		} else if rem, b2, ok := parseTrailingBounds(eng, body, 2); ok {
			body = rem
			xMin, xMax = b2[0], b2[1]
		}
		out, err := eng.PlotSlopeField(body, xMin, xMax, yMin, yMax, 64, 18, true)
		if err != nil {
			fmt.Printf("Slope field error: %v\n\n", err)
		} else {
			fmt.Println(out)
		}
		return true
	}

	if strings.HasPrefix(lower, "table ") {
		body := strings.TrimSpace(input[len("table "):])
		xMin, xMax, step := -5.0, 5.0, 1.0
		if rem, b3, ok := parseTrailingBounds(eng, body, 3); ok {
			body = rem
			xMin, xMax, step = b3[0], b3[1], b3[2]
		} else if rem, b2, ok := parseTrailingBounds(eng, body, 2); ok {
			body = rem
			xMin, xMax = b2[0], b2[1]
			step = (xMax - xMin) / 10.0
		}
		exprs := splitTopLevel(body, ',')
		out, err := eng.RenderTable(exprs, xMin, xMax, step)
		if err != nil {
			fmt.Printf("Table error: %v\n\n", err)
		} else {
			fmt.Println(out)
		}
		return true
	}

	return false
}

func executeLine(eng *Engine, line string) {
	stmts := splitTopLevel(line, ';')
	var lastVal Value
	hasVal := false
	for _, stmt := range stmts {
		if stmt == "" {
			continue
		}
		if handleCommand(eng, stmt) {
			hasVal = false
			continue
		}
		val, err := eng.Eval(stmt)
		if err != nil {
			fmt.Printf("Error     = %v\n\n", err)
			return
		}
		lastVal = val
		hasVal = true
	}
	if hasVal {
		printFormattedResult(lastVal)
	}
}

func main() {
	if len(os.Args) > 1 {
		joined := strings.TrimSpace(strings.Join(os.Args[1:], " "))
		executeLine(defaultEngine, joined)
		return
	}

	scanner := bufio.NewScanner(os.Stdin)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	printBanner()

	for {
		if defaultEngine.Mode == ModeDeg {
			fmt.Print("Calculate [deg] > ")
		} else {
			fmt.Print("Calculate > ")
		}
		if !scanner.Scan() {
			break
		}

		input := strings.TrimSpace(scanner.Text())
		if input == "" {
			continue
		}

		lower := strings.ToLower(input)
		if lower == "exit" || lower == "quit" {
			break
		}

		executeLine(defaultEngine, input)
	}
}
