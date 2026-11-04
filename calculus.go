package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"strings"
)

func factorialReal(x float64) float64 {
	if math.IsNaN(x) || math.IsInf(x, -1) {
		return math.NaN()
	}
	if math.IsInf(x, 1) {
		return math.Inf(1)
	}
	if x < 0 {
		if x == math.Trunc(x) {
			return math.NaN()
		}
		return math.Gamma(x + 1)
	}
	n := math.Floor(x)
	if n != x {
		return math.Gamma(x + 1)
	}
	if n > 170 {
		return math.Inf(1)
	}
	res := 1.0
	for i := 2.0; i <= n; i++ {
		res *= i
	}
	return res
}

func doubleFactorial(x float64) float64 {
	if x < -1 || x != math.Trunc(x) {
		return math.NaN()
	}
	n := int64(x)
	if n <= 0 {
		return 1
	}
	res := 1.0
	for k := n; k >= 2; k -= 2 {
		res *= float64(k)
		if math.IsInf(res, 1) {
			return res
		}
	}
	return res
}

var lanczosCoeffs = []float64{
	0.99999999999980993,
	676.5203681218851,
	-1259.1392167224028,
	771.32342877765313,
	-176.61502916214059,
	12.507343278686905,
	-0.13857109526572012,
	9.9843695780195716e-6,
	1.5056327351493116e-7,
}

func gammaComplex(z complex128) complex128 {
	if imag(z) == 0 {
		x := real(z)
		if x <= 0 && x == math.Trunc(x) {
			return complex(math.NaN(), 0)
		}
		return complex(math.Gamma(x), 0)
	}
	if real(z) < 0.5 {
		return complex(math.Pi, 0) / (cmplx.Sin(complex(math.Pi, 0)*z) * gammaComplex(1-z))
	}
	z -= 1
	x := complex(lanczosCoeffs[0], 0)
	for i := 1; i < len(lanczosCoeffs); i++ {
		x += complex(lanczosCoeffs[i], 0) / (z + complex(float64(i), 0))
	}
	t := z + complex(7.5, 0)
	return complex(math.Sqrt(2*math.Pi), 0) * cmplx.Pow(t, z+0.5) * cmplx.Exp(-t) * x
}

func digammaReal(x float64) float64 {
	if x <= 0 && x == math.Trunc(x) {
		return math.NaN()
	}
	res := 0.0
	for x < 7 {
		res -= 1.0 / x
		x += 1.0
	}
	x -= 0.5
	xx := 1.0 / x
	xx2 := xx * xx
	xx4 := xx2 * xx2
	res += math.Log(x) + (1.0/24.0)*xx2 - (7.0/960.0)*xx4 + (31.0/8064.0)*xx4*xx2 - (127.0/30720.0)*xx4*xx4
	return res
}

func zetaComplex(s complex128) complex128 {
	if cmplx.Abs(s-1) < 1e-15 {
		return complex(math.Inf(1), 0)
	}
	if s == 0 {
		return complex(-0.5, 0)
	}
	if imag(s) == 0 && real(s) < 0 {
		rs := real(s)
		if math.Mod(-rs, 2) == 0 {
			return 0
		}
	}
	if real(s) < 0.5 {
		factor := cmplx.Pow(2, s) * cmplx.Pow(complex(math.Pi, 0), s-1) *
			cmplx.Sin(complex(math.Pi/2, 0)*s) * gammaComplex(1-s)
		return cleanComplex(factor * zetaComplex(1-s))
	}

	const n = 28
	d := make([]float64, n+1)
	sumTerm := 0.0
	for k := 0; k <= n; k++ {
		term := factorialReal(float64(n+k-1)) * math.Pow(4, float64(k)) /
			(factorialReal(float64(n-k)) * factorialReal(float64(2*k)))
		sumTerm += term
		d[k] = float64(n) * sumTerm
	}
	dn := d[n]
	var acc complex128
	for k := 0; k < n; k++ {
		sign := 1.0
		if k%2 == 1 {
			sign = -1.0
		}
		weight := complex(-sign*(d[k]-dn)/dn, 0)
		denom := cmplx.Pow(complex(float64(k+1), 0), s)
		acc += weight / denom
	}
	denomEta := 1 - cmplx.Pow(2, 1-s)
	res := acc / denomEta
	if cmplx.Abs(res) <= 1e-14 {
		return 0
	}
	return cleanComplex(res)
}

func lambertW(x float64, branch int) (float64, error) {
	invE := -1.0 / math.E
	if branch == 0 {
		if x < invE-1e-15 {
			return 0, fmt.Errorf("lambertw principal branch undefined for x < -1/e")
		}
		if math.Abs(x-invE) <= 1e-15 {
			return -1, nil
		}
		if x == 0 {
			return 0, nil
		}
		var w float64
		if x < 1 {
			p := math.Sqrt(2 * (math.E*x + 1))
			w = -1 + p - p*p/3.0
		} else {
			lx := math.Log(x)
			w = lx - math.Log(math.Max(1.0, lx))
		}
		return halleyLambertW(x, w), nil
	}
	if branch == -1 {
		if x < invE-1e-15 || x >= 0 {
			return 0, fmt.Errorf("lambertw(-1) branch requires -1/e <= x < 0")
		}
		if math.Abs(x-invE) <= 1e-15 {
			return -1, nil
		}
		lx := math.Log(-x)
		w := lx - math.Log(-lx)
		return halleyLambertW(x, w), nil
	}
	return 0, fmt.Errorf("unsupported lambertw branch %d (use 0 or -1)", branch)
}

func halleyLambertW(x, w float64) float64 {
	for iter := 0; iter < 50; iter++ {
		ew := math.Exp(w)
		wew := w * ew
		f := wew - x
		if math.Abs(f) <= 1e-15*(1+math.Abs(x)) {
			break
		}
		denom := ew*(w+1) - (w+2)*f/(2*w+2)
		if denom == 0 {
			break
		}
		dw := f / denom
		w -= dw
		if math.Abs(dw) <= 1e-15*(1+math.Abs(w)) {
			break
		}
	}
	return w
}

func erfInv(x float64) float64 {
	if x <= -1 {
		if x == -1 {
			return math.Inf(-1)
		}
		return math.NaN()
	}
	if x >= 1 {
		if x == 1 {
			return math.Inf(1)
		}
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	a := 0.147
	ln1x2 := math.Log(1 - x*x)
	t1 := 2/(math.Pi*a) + ln1x2/2
	t2 := ln1x2 / a
	sign := 1.0
	if x < 0 {
		sign = -1.0
	}
	y := sign * math.Sqrt(math.Sqrt(t1*t1-t2)-t1)
	for i := 0; i < 4; i++ {
		errVal := math.Erf(y) - x
		deriv := (2.0 / math.SqrtPi) * math.Exp(-y*y)
		y -= errVal / (deriv - y*errVal)
	}
	return y
}

func regGammaP(a, x float64) float64 {
	if x <= 0 || a <= 0 {
		return 0
	}
	if x < a+1.0 {
		ap := a
		sum := 1.0 / a
		del := sum
		for n := 1; n <= 200; n++ {
			ap += 1.0
			del *= x / ap
			sum += del
			if math.Abs(del) < math.Abs(sum)*1e-15 {
				break
			}
		}
		lg, _ := math.Lgamma(a)
		return sum * math.Exp(-x+a*math.Log(x)-lg)
	}
	return 1.0 - regGammaQCF(a, x)
}

func regGammaQCF(a, x float64) float64 {
	lg, _ := math.Lgamma(a)
	b := x + 1.0 - a
	c := 1.0 / 1e-30
	d := 1.0 / b
	h := d
	for i := 1; i <= 200; i++ {
		an := -float64(i) * (float64(i) - a)
		b += 2.0
		d = an*d + b
		if math.Abs(d) < 1e-30 {
			d = 1e-30
		}
		c = b + an/c
		if math.Abs(c) < 1e-30 {
			c = 1e-30
		}
		d = 1.0 / d
		del := d * c
		h *= del
		if math.Abs(del-1.0) < 1e-15 {
			break
		}
	}
	return math.Exp(-x+a*math.Log(x)-lg) * h
}

func regBetaInc(a, b, x float64) float64 {
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}
	lga, _ := math.Lgamma(a)
	lgb, _ := math.Lgamma(b)
	lgab, _ := math.Lgamma(a + b)
	bt := math.Exp(lgab - lga - lgb + a*math.Log(x) + b*math.Log(1-x))
	if x < (a+1.0)/(a+b+2.0) {
		return bt * betaCF(a, b, x) / a
	}
	return 1.0 - bt*betaCF(b, a, 1.0-x)/b
}

func betaCF(a, b, x float64) float64 {
	qab := a + b
	qap := a + 1.0
	qam := a - 1.0
	c := 1.0
	d := 1.0 - qab*x/qap
	if math.Abs(d) < 1e-30 {
		d = 1e-30
	}
	d = 1.0 / d
	h := d
	for m := 1; m <= 200; m++ {
		mf := float64(m)
		m2 := 2.0 * mf
		aa := mf * (b - mf) * x / ((qam + m2) * (a + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < 1e-30 {
			d = 1e-30
		}
		c = 1.0 + aa/c
		if math.Abs(c) < 1e-30 {
			c = 1e-30
		}
		d = 1.0 / d
		h *= d * c
		aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < 1e-30 {
			d = 1e-30
		}
		c = 1.0 + aa/c
		if math.Abs(c) < 1e-30 {
			c = 1e-30
		}
		d = 1.0 / d
		del := d * c
		h *= del
		if math.Abs(del-1.0) < 1e-15 {
			break
		}
	}
	return h
}

func evalSpecialMathFunc(fn string, args []Value) (bool, Value, error) {
	switch fn {
	case "gamma":
		if len(args) != 1 || args[0].Kind != KindScalar {
			return true, Value{}, fmt.Errorf("gamma(z) requires 1 scalar argument")
		}
		return true, ScalarComplex(gammaComplex(args[0].S)), nil

	case "lgamma":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("lgamma(x) requires 1 real argument")
		}
		r, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		lg, _ := math.Lgamma(r)
		return true, ScalarReal(lg), nil

	case "digamma", "psi":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("digamma(x) requires 1 real argument")
		}
		r, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(digammaReal(r)), nil

	case "beta":
		if len(args) != 2 || args[0].Kind != KindScalar || args[1].Kind != KindScalar {
			return true, Value{}, fmt.Errorf("beta(a, b) requires 2 scalar arguments")
		}
		a, b := args[0].S, args[1].S
		return true, ScalarComplex(gammaComplex(a) * gammaComplex(b) / gammaComplex(a+b)), nil

	case "zeta":
		if len(args) != 1 || args[0].Kind != KindScalar {
			return true, Value{}, fmt.Errorf("zeta(s) requires 1 scalar argument")
		}
		return true, ScalarComplex(zetaComplex(args[0].S)), nil

	case "lambertw", "productlog":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("lambertw(x [, branch]) takes 1 or 2 arguments")
		}
		x, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		branch := 0
		if len(args) == 2 {
			br, err := args[1].Real()
			if err != nil {
				return true, Value{}, err
			}
			branch = int(math.Round(br))
		}
		w, err := lambertW(x, branch)
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(w), nil

	case "erf":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("erf(x) takes 1 argument")
		}
		r, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(math.Erf(r)), nil

	case "erfc":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("erfc(x) takes 1 argument")
		}
		r, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(math.Erfc(r)), nil

	case "erfinv":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("erfinv(x) takes 1 argument")
		}
		r, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(erfInv(r)), nil

	case "besselj":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("besselj(n, x) takes 2 real arguments")
		}
		n, err1 := args[0].Real()
		x, err2 := args[1].Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("besselj requires real arguments")
		}
		return true, ScalarReal(math.Jn(int(math.Round(n)), x)), nil

	case "bessely":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("bessely(n, x) takes 2 real arguments")
		}
		n, err1 := args[0].Real()
		x, err2 := args[1].Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("bessely requires real arguments")
		}
		return true, ScalarReal(math.Yn(int(math.Round(n)), x)), nil

	case "sinc":
		if len(args) != 1 || args[0].Kind != KindScalar {
			return true, Value{}, fmt.Errorf("sinc(z) requires 1 scalar argument")
		}
		z := args[0].S
		if cmplx.Abs(z) < 1e-12 {
			return true, ScalarReal(1), nil
		}
		return true, ScalarComplex(cmplx.Sin(z) / z), nil

	case "gammainc":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("gammainc(a, x) requires 2 real arguments")
		}
		a, err1 := args[0].Real()
		x, err2 := args[1].Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("gammainc requires real arguments")
		}
		return true, ScalarReal(regGammaP(a, x)), nil

	case "betainc":
		if len(args) != 3 {
			return true, Value{}, fmt.Errorf("betainc(a, b, x) requires 3 real arguments")
		}
		a, err1 := args[0].Real()
		b, err2 := args[1].Real()
		x, err3 := args[2].Real()
		if err1 != nil || err2 != nil || err3 != nil {
			return true, Value{}, fmt.Errorf("betainc requires real arguments")
		}
		return true, ScalarReal(regBetaInc(a, b, x)), nil
	}

	return false, Value{}, nil
}

func (e *Engine) expandSymbolicRefs(node ASTNode) ASTNode {
	switch n := node.(type) {
	case IdentNode:
		if n.Name == "ans" || n.Name == "_" {
			if e.LastValue.Kind == KindSymbolic && e.LastValue.Node != nil {
				return e.LastValue.Node
			}
			if e.LastValue.Kind == KindScalar && e.LastValue.IsReal() {
				return NumNode{Val: real(e.LastValue.S)}
			}
		}
		if v, ok := e.Vars[n.Name]; ok {
			if v.Kind == KindSymbolic && v.Node != nil {
				return v.Node
			}
		}
		return n
	case UnaryNode:
		return UnaryNode{Op: n.Op, X: e.expandSymbolicRefs(n.X)}
	case BinaryNode:
		return BinaryNode{Op: n.Op, Left: e.expandSymbolicRefs(n.Left), Right: e.expandSymbolicRefs(n.Right)}
	case PostfixNode:
		return PostfixNode{Op: n.Op, X: e.expandSymbolicRefs(n.X)}
	case CallNode:
		args := make([]ASTNode, len(n.Args))
		for i, a := range n.Args {
			args[i] = e.expandSymbolicRefs(a)
		}
		if uf, ok := e.Funcs[n.Name]; ok && len(args) == len(uf.Params) && n.PrimeCount == 0 {
			body := uf.Body
			for i, p := range uf.Params {
				body = substituteAST(body, p, args[i])
			}
			return e.expandSymbolicRefs(body)
		}
		return CallNode{Name: n.Name, PrimeCount: n.PrimeCount, Args: args}
	default:
		return node
	}
}

func dependsOn(node ASTNode, varName string) bool {
	switch n := node.(type) {
	case IdentNode:
		return n.Name == varName
	case UnaryNode:
		return dependsOn(n.X, varName)
	case BinaryNode:
		return dependsOn(n.Left, varName) || dependsOn(n.Right, varName)
	case PostfixNode:
		return dependsOn(n.X, varName)
	case CallNode:
		if n.Name == varName {
			return true
		}
		for _, a := range n.Args {
			if dependsOn(a, varName) {
				return true
			}
		}
		return false
	default:
		return false
	}
}

func substituteAST(node ASTNode, varName string, rep ASTNode) ASTNode {
	switch n := node.(type) {
	case IdentNode:
		if n.Name == varName {
			return rep
		}
		return n
	case UnaryNode:
		return UnaryNode{Op: n.Op, X: substituteAST(n.X, varName, rep)}
	case BinaryNode:
		return BinaryNode{
			Op:    n.Op,
			Left:  substituteAST(n.Left, varName, rep),
			Right: substituteAST(n.Right, varName, rep),
		}
	case PostfixNode:
		return PostfixNode{Op: n.Op, X: substituteAST(n.X, varName, rep)}
	case CallNode:
		args := make([]ASTNode, len(n.Args))
		for i, a := range n.Args {
			args[i] = substituteAST(a, varName, rep)
		}
		return CallNode{Name: n.Name, PrimeCount: n.PrimeCount, Args: args}
	default:
		return node
	}
}

func symbolicDiff(node ASTNode, v string) (ASTNode, error) {
	if !dependsOn(node, v) {
		return NumNode{Val: 0}, nil
	}

	switch n := node.(type) {
	case IdentNode:
		if n.Name == v {
			return NumNode{Val: 1}, nil
		}
		return NumNode{Val: 0}, nil

	case UnaryNode:
		dx, err := symbolicDiff(n.X, v)
		if err != nil {
			return nil, err
		}
		if n.Op == "-" {
			return UnaryNode{Op: "-", X: dx}, nil
		}
		return dx, nil

	case BinaryNode:
		switch n.Op {
		case "+", "-":
			dl, err := symbolicDiff(n.Left, v)
			if err != nil {
				return nil, err
			}
			dr, err := symbolicDiff(n.Right, v)
			if err != nil {
				return nil, err
			}
			return BinaryNode{Op: n.Op, Left: dl, Right: dr}, nil

		case "*":
			dl, err := symbolicDiff(n.Left, v)
			if err != nil {
				return nil, err
			}
			dr, err := symbolicDiff(n.Right, v)
			if err != nil {
				return nil, err
			}
			return BinaryNode{
				Op:    "+",
				Left:  BinaryNode{Op: "*", Left: dl, Right: n.Right},
				Right: BinaryNode{Op: "*", Left: n.Left, Right: dr},
			}, nil

		case "/":
			dl, err := symbolicDiff(n.Left, v)
			if err != nil {
				return nil, err
			}
			if !dependsOn(n.Right, v) {
				return BinaryNode{Op: "/", Left: dl, Right: n.Right}, nil
			}
			dr, err := symbolicDiff(n.Right, v)
			if err != nil {
				return nil, err
			}
			num := BinaryNode{
				Op:    "-",
				Left:  BinaryNode{Op: "*", Left: dl, Right: n.Right},
				Right: BinaryNode{Op: "*", Left: n.Left, Right: dr},
			}
			den := BinaryNode{Op: "^", Left: n.Right, Right: NumNode{Val: 2}}
			return BinaryNode{Op: "/", Left: num, Right: den}, nil

		case "^":
			u, w := n.Left, n.Right
			uDep := dependsOn(u, v)
			wDep := dependsOn(w, v)
			if uDep && !wDep {
				du, err := symbolicDiff(u, v)
				if err != nil {
					return nil, err
				}
				expMinus1 := BinaryNode{Op: "-", Left: w, Right: NumNode{Val: 1}}
				return BinaryNode{
					Op: "*",
					Left: BinaryNode{
						Op:    "*",
						Left:  w,
						Right: BinaryNode{Op: "^", Left: u, Right: expMinus1},
					},
					Right: du,
				}, nil
			}
			if !uDep && wDep {
				dw, err := symbolicDiff(w, v)
				if err != nil {
					return nil, err
				}
				if id, ok := u.(IdentNode); ok && strings.ToLower(id.Name) == "e" {
					return BinaryNode{Op: "*", Left: n, Right: dw}, nil
				}
				lnU := CallNode{Name: "ln", Args: []ASTNode{u}}
				return BinaryNode{
					Op:    "*",
					Left:  BinaryNode{Op: "*", Left: n, Right: lnU},
					Right: dw,
				}, nil
			}
			du, err := symbolicDiff(u, v)
			if err != nil {
				return nil, err
			}
			dw, err := symbolicDiff(w, v)
			if err != nil {
				return nil, err
			}
			term1 := BinaryNode{Op: "*", Left: dw, Right: CallNode{Name: "ln", Args: []ASTNode{u}}}
			term2 := BinaryNode{Op: "/", Left: BinaryNode{Op: "*", Left: w, Right: du}, Right: u}
			return BinaryNode{Op: "*", Left: n, Right: BinaryNode{Op: "+", Left: term1, Right: term2}}, nil
		}

	case CallNode:
		fn := strings.ToLower(n.Name)
		if len(n.Args) == 1 {
			u := n.Args[0]
			du, err := symbolicDiff(u, v)
			if err != nil {
				return nil, err
			}
			var outer ASTNode
			switch fn {
			case "sin":
				outer = CallNode{Name: "cos", Args: []ASTNode{u}}
			case "cos":
				outer = UnaryNode{Op: "-", X: CallNode{Name: "sin", Args: []ASTNode{u}}}
			case "tan":
				outer = BinaryNode{Op: "^", Left: CallNode{Name: "sec", Args: []ASTNode{u}}, Right: NumNode{Val: 2}}
			case "sec":
				outer = BinaryNode{Op: "*", Left: CallNode{Name: "sec", Args: []ASTNode{u}}, Right: CallNode{Name: "tan", Args: []ASTNode{u}}}
			case "csc":
				outer = UnaryNode{Op: "-", X: BinaryNode{Op: "*", Left: CallNode{Name: "csc", Args: []ASTNode{u}}, Right: CallNode{Name: "cot", Args: []ASTNode{u}}}}
			case "cot":
				outer = UnaryNode{Op: "-", X: BinaryNode{Op: "^", Left: CallNode{Name: "csc", Args: []ASTNode{u}}, Right: NumNode{Val: 2}}}
			case "asin", "arcsin":
				denom := CallNode{Name: "sqrt", Args: []ASTNode{BinaryNode{Op: "-", Left: NumNode{Val: 1}, Right: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "acos", "arccos":
				denom := CallNode{Name: "sqrt", Args: []ASTNode{BinaryNode{Op: "-", Left: NumNode{Val: 1}, Right: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}}}}
				return UnaryNode{Op: "-", X: BinaryNode{Op: "/", Left: du, Right: denom}}, nil
			case "atan", "arctan":
				denom := BinaryNode{Op: "+", Left: NumNode{Val: 1}, Right: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "sinh":
				outer = CallNode{Name: "cosh", Args: []ASTNode{u}}
			case "cosh":
				outer = CallNode{Name: "sinh", Args: []ASTNode{u}}
			case "tanh":
				outer = BinaryNode{Op: "^", Left: CallNode{Name: "sech", Args: []ASTNode{u}}, Right: NumNode{Val: 2}}
			case "asinh", "arcsinh":
				denom := CallNode{Name: "sqrt", Args: []ASTNode{BinaryNode{Op: "+", Left: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}, Right: NumNode{Val: 1}}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "acosh", "arccosh":
				denom := CallNode{Name: "sqrt", Args: []ASTNode{BinaryNode{Op: "-", Left: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}, Right: NumNode{Val: 1}}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "atanh", "arctanh":
				denom := BinaryNode{Op: "-", Left: NumNode{Val: 1}, Right: BinaryNode{Op: "^", Left: u, Right: NumNode{Val: 2}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "exp":
				outer = CallNode{Name: "exp", Args: []ASTNode{u}}
			case "ln":
				return BinaryNode{Op: "/", Left: du, Right: u}, nil
			case "log", "log10":
				denom := BinaryNode{Op: "*", Left: u, Right: CallNode{Name: "ln", Args: []ASTNode{NumNode{Val: 10}}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "log2":
				denom := BinaryNode{Op: "*", Left: u, Right: CallNode{Name: "ln", Args: []ASTNode{NumNode{Val: 2}}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "sqrt":
				denom := BinaryNode{Op: "*", Left: NumNode{Val: 2}, Right: CallNode{Name: "sqrt", Args: []ASTNode{u}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "cbrt":
				denom := BinaryNode{Op: "*", Left: NumNode{Val: 3}, Right: BinaryNode{Op: "^", Left: CallNode{Name: "cbrt", Args: []ASTNode{u}}, Right: NumNode{Val: 2}}}
				return BinaryNode{Op: "/", Left: du, Right: denom}, nil
			case "abs":
				outer = BinaryNode{Op: "/", Left: u, Right: CallNode{Name: "abs", Args: []ASTNode{u}}}
			default:
				if n.Name == v {
					return BinaryNode{
						Op:    "+",
						Left:  u,
						Right: BinaryNode{Op: "*", Left: IdentNode{Name: v}, Right: du},
					}, nil
				}
				return nil, fmt.Errorf("cannot symbolically differentiate %s", n.Name)
			}
			return BinaryNode{Op: "*", Left: outer, Right: du}, nil
		}
		if fn == "pow" && len(n.Args) == 2 {
			return symbolicDiff(BinaryNode{Op: "^", Left: n.Args[0], Right: n.Args[1]}, v)
		}
	}

	return nil, fmt.Errorf("cannot symbolically differentiate %s", node.Format())
}

func isNumVal(n ASTNode, target float64) bool {
	if num, ok := n.(NumNode); ok {
		return math.Abs(num.Val-target) <= 1e-14
	}
	return false
}

func getNum(n ASTNode) (float64, bool) {
	if num, ok := n.(NumNode); ok {
		return num.Val, true
	}
	if un, ok := n.(UnaryNode); ok && un.Op == "-" {
		if num, ok := un.X.(NumNode); ok {
			return -num.Val, true
		}
	}
	return 0, false
}

func astEqual(a, b ASTNode) bool {
	return a.Format() == b.Format()
}

func gcdInt64(a, b int64) int64 {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

func simplifyAST(node ASTNode) ASTNode {
	for pass := 0; pass < 4; pass++ {
		next := simplifyPass(node)
		if next.Format() == node.Format() {
			return next
		}
		node = next
	}
	return node
}

func makeSignedNum(v float64) ASTNode {
	if v < 0 {
		return UnaryNode{Op: "-", X: NumNode{Val: -v}}
	}
	return NumNode{Val: v}
}

func extractLeadingNeg(n ASTNode) (ASTNode, bool) {
	if un, ok := n.(UnaryNode); ok && un.Op == "-" {
		return un.X, true
	}
	if bin, ok := n.(BinaryNode); ok && (bin.Op == "*" || bin.Op == "/") {
		if posL, neg := extractLeadingNeg(bin.Left); neg {
			return BinaryNode{Op: bin.Op, Left: posL, Right: bin.Right}, true
		}
	}
	return nil, false
}

func splitCoeffTerm(n ASTNode) (float64, ASTNode) {
	if val, ok := getNum(n); ok {
		return val, NumNode{Val: 1}
	}
	if un, ok := n.(UnaryNode); ok && un.Op == "-" {
		c, base := splitCoeffTerm(un.X)
		return -c, base
	}
	if bin, ok := n.(BinaryNode); ok && bin.Op == "*" {
		if c, ok := getNum(bin.Left); ok {
			return c, bin.Right
		}
		if c, ok := getNum(bin.Right); ok {
			return c, bin.Left
		}
	}
	return 1, n
}

func simplifyPass(node ASTNode) ASTNode {
	switch n := node.(type) {
	case UnaryNode:
		x := simplifyPass(n.X)
		if n.Op == "-" {
			if isNumVal(x, 0) {
				return NumNode{Val: 0}
			}
			if inner, ok := x.(UnaryNode); ok && inner.Op == "-" {
				return inner.X
			}
		}
		return UnaryNode{Op: n.Op, X: x}

	case BinaryNode:
		l := simplifyPass(n.Left)
		r := simplifyPass(n.Right)
		lv, lNum := getNum(l)
		rv, rNum := getNum(r)

		switch n.Op {
		case "+":
			if isNumVal(l, 0) {
				return r
			}
			if isNumVal(r, 0) {
				return l
			}
			if lNum && rNum {
				return makeSignedNum(lv + rv)
			}
			if posR, neg := extractLeadingNeg(r); neg {
				return simplifyPass(BinaryNode{Op: "-", Left: l, Right: posR})
			}
			c1, t1 := splitCoeffTerm(l)
			c2, t2 := splitCoeffTerm(r)
			if astEqual(t1, t2) && !isNumVal(t1, 1) {
				sumC := c1 + c2
				if math.Abs(sumC) <= 1e-14 {
					return NumNode{Val: 0}
				}
				if math.Abs(sumC-1) <= 1e-14 {
					return t1
				}
				if math.Abs(sumC+1) <= 1e-14 {
					return UnaryNode{Op: "-", X: t1}
				}
				return simplifyPass(BinaryNode{Op: "*", Left: makeSignedNum(sumC), Right: t1})
			}

		case "-":
			if isNumVal(r, 0) {
				return l
			}
			if isNumVal(l, 0) {
				return simplifyPass(UnaryNode{Op: "-", X: r})
			}
			if lNum && rNum {
				return makeSignedNum(lv - rv)
			}
			if astEqual(l, r) {
				return NumNode{Val: 0}
			}
			if posR, neg := extractLeadingNeg(r); neg {
				return simplifyPass(BinaryNode{Op: "+", Left: l, Right: posR})
			}
			c1, t1 := splitCoeffTerm(l)
			c2, t2 := splitCoeffTerm(r)
			if astEqual(t1, t2) && !isNumVal(t1, 1) {
				diffC := c1 - c2
				if math.Abs(diffC) <= 1e-14 {
					return NumNode{Val: 0}
				}
				return simplifyPass(BinaryNode{Op: "*", Left: makeSignedNum(diffC), Right: t1})
			}

		case "*":
			if isNumVal(l, 0) || isNumVal(r, 0) {
				return NumNode{Val: 0}
			}
			if isNumVal(l, 1) {
				return r
			}
			if isNumVal(r, 1) {
				return l
			}
			if lNum && math.Abs(lv+1) <= 1e-14 {
				return simplifyPass(UnaryNode{Op: "-", X: r})
			}
			if rNum && math.Abs(rv+1) <= 1e-14 {
				return simplifyPass(UnaryNode{Op: "-", X: l})
			}
			if lNum && rNum {
				return makeSignedNum(lv * rv)
			}
			if rNum && !lNum {
				l, r = r, l
				lv, lNum = rv, true
				rv, rNum = 0, false
			}
			if lNum {
				if rb, ok := r.(BinaryNode); ok {
					if rb.Op == "*" {
						if rc, ok := getNum(rb.Left); ok {
							return simplifyPass(BinaryNode{Op: "*", Left: makeSignedNum(lv * rc), Right: rb.Right})
						}
					}
					if rb.Op == "/" {
						return simplifyPass(BinaryNode{Op: "/", Left: BinaryNode{Op: "*", Left: l, Right: rb.Left}, Right: rb.Right})
					}
				}
			}
			if astEqual(l, r) {
				return BinaryNode{Op: "^", Left: l, Right: NumNode{Val: 2}}
			}
			if lb, ok := l.(BinaryNode); ok && lb.Op == "/" && astEqual(lb.Right, r) {
				return lb.Left
			}
			if rb, ok := r.(BinaryNode); ok && rb.Op == "/" && astEqual(rb.Right, l) {
				return rb.Left
			}

		case "/":
			if isNumVal(l, 0) {
				return NumNode{Val: 0}
			}
			if isNumVal(r, 1) {
				return l
			}
			if astEqual(l, r) {
				return NumNode{Val: 1}
			}
			if lNum && rNum && rv != 0 {
				if lv == math.Trunc(lv) && rv == math.Trunc(rv) && math.Abs(lv) < 1e12 && math.Abs(rv) < 1e12 {
					li, ri := int64(lv), int64(rv)
					if ri < 0 {
						li, ri = -li, -ri
					}
					g := gcdInt64(li, ri)
					li /= g
					ri /= g
					if ri == 1 {
						return makeSignedNum(float64(li))
					}
					return BinaryNode{Op: "/", Left: makeSignedNum(float64(li)), Right: NumNode{Val: float64(ri)}}
				}
				return makeSignedNum(lv / rv)
			}
			if lb, ok := l.(BinaryNode); ok && lb.Op == "*" {
				if astEqual(lb.Left, r) {
					return lb.Right
				}
				if astEqual(lb.Right, r) {
					return lb.Left
				}
				if lc, ok := getNum(lb.Left); ok && rNum && lc == math.Trunc(lc) && rv == math.Trunc(rv) {
					g := gcdInt64(int64(lc), int64(rv))
					if g > 1 {
						nl := lc / float64(g)
						nr := rv / float64(g)
						return simplifyPass(BinaryNode{
							Op:    "/",
							Left:  BinaryNode{Op: "*", Left: makeSignedNum(nl), Right: lb.Right},
							Right: makeSignedNum(nr),
						})
					}
				}
			}

		case "^":
			if isNumVal(r, 0) {
				return NumNode{Val: 1}
			}
			if isNumVal(r, 1) {
				return l
			}
			if isNumVal(l, 1) {
				return NumNode{Val: 1}
			}
			if lNum && rNum && (lv >= 0 || rv == math.Trunc(rv)) {
				p := math.Pow(lv, rv)
				if !math.IsNaN(p) && !math.IsInf(p, 0) && p == math.Trunc(p) && math.Abs(p) < 1e12 {
					return makeSignedNum(p)
				}
			}
			if lb, ok := l.(BinaryNode); ok && lb.Op == "^" {
				return simplifyPass(BinaryNode{Op: "^", Left: lb.Left, Right: BinaryNode{Op: "*", Left: lb.Right, Right: r}})
			}
		}
		return BinaryNode{Op: n.Op, Left: l, Right: r}

	case CallNode:
		args := make([]ASTNode, len(n.Args))
		for i, a := range n.Args {
			args[i] = simplifyPass(a)
		}
		fn := strings.ToLower(n.Name)
		if len(args) == 1 {
			if id, ok := args[0].(IdentNode); ok && strings.ToLower(id.Name) == "e" && fn == "ln" {
				return NumNode{Val: 1}
			}
			if v, ok := getNum(args[0]); ok {
				switch fn {
				case "ln", "log", "log10", "log2":
					if v == 1 {
						return NumNode{Val: 0}
					}
				case "exp":
					if v == 0 {
						return NumNode{Val: 1}
					}
				case "sin", "tan", "sinh", "tanh":
					if v == 0 {
						return NumNode{Val: 0}
					}
				case "cos", "cosh":
					if v == 0 {
						return NumNode{Val: 1}
					}
				case "sqrt":
					if v >= 0 {
						rt := math.Sqrt(v)
						if rt == math.Trunc(rt) {
							return NumNode{Val: rt}
						}
					}
				}
			}
		}
		return CallNode{Name: n.Name, PrimeCount: n.PrimeCount, Args: args}
	}

	return node
}

func extractLinear(node ASTNode, v string) (a, b ASTNode, ok bool) {
	if !dependsOn(node, v) {
		return NumNode{Val: 0}, node, true
	}
	if id, isId := node.(IdentNode); isId && id.Name == v {
		return NumNode{Val: 1}, NumNode{Val: 0}, true
	}
	if un, isUn := node.(UnaryNode); isUn && un.Op == "-" {
		sa, sb, sok := extractLinear(un.X, v)
		if !sok {
			return nil, nil, false
		}
		return simplifyAST(UnaryNode{Op: "-", X: sa}), simplifyAST(UnaryNode{Op: "-", X: sb}), true
	}
	if bin, isBin := node.(BinaryNode); isBin {
		switch bin.Op {
		case "+", "-":
			la, lb, lok := extractLinear(bin.Left, v)
			ra, rb, rok := extractLinear(bin.Right, v)
			if !lok || !rok {
				return nil, nil, false
			}
			return simplifyAST(BinaryNode{Op: bin.Op, Left: la, Right: ra}),
				simplifyAST(BinaryNode{Op: bin.Op, Left: lb, Right: rb}), true
		case "*":
			if !dependsOn(bin.Left, v) {
				ra, rb, rok := extractLinear(bin.Right, v)
				if !rok {
					return nil, nil, false
				}
				return simplifyAST(BinaryNode{Op: "*", Left: bin.Left, Right: ra}),
					simplifyAST(BinaryNode{Op: "*", Left: bin.Left, Right: rb}), true
			}
			if !dependsOn(bin.Right, v) {
				la, lb, lok := extractLinear(bin.Left, v)
				if !lok {
					return nil, nil, false
				}
				return simplifyAST(BinaryNode{Op: "*", Left: la, Right: bin.Right}),
					simplifyAST(BinaryNode{Op: "*", Left: lb, Right: bin.Right}), true
			}
		case "/":
			if !dependsOn(bin.Right, v) {
				la, lb, lok := extractLinear(bin.Left, v)
				if !lok {
					return nil, nil, false
				}
				return simplifyAST(BinaryNode{Op: "/", Left: la, Right: bin.Right}),
					simplifyAST(BinaryNode{Op: "/", Left: lb, Right: bin.Right}), true
			}
		}
	}
	return nil, nil, false
}

func symbolicAntideriv(node ASTNode, v string) (ASTNode, error) {
	node = simplifyAST(node)
	if !dependsOn(node, v) {
		return BinaryNode{Op: "*", Left: node, Right: IdentNode{Name: v}}, nil
	}

	switch n := node.(type) {
	case IdentNode:
		if n.Name == v {
			return BinaryNode{
				Op:    "/",
				Left:  BinaryNode{Op: "^", Left: IdentNode{Name: v}, Right: NumNode{Val: 2}},
				Right: NumNode{Val: 2},
			}, nil
		}

	case UnaryNode:
		if n.Op == "-" {
			sub, err := symbolicAntideriv(n.X, v)
			if err != nil {
				return nil, err
			}
			return UnaryNode{Op: "-", X: sub}, nil
		}

	case BinaryNode:
		switch n.Op {
		case "+", "-":
			li, err := symbolicAntideriv(n.Left, v)
			if err != nil {
				return nil, err
			}
			ri, err := symbolicAntideriv(n.Right, v)
			if err != nil {
				return nil, err
			}
			return BinaryNode{Op: n.Op, Left: li, Right: ri}, nil

		case "*":
			if !dependsOn(n.Left, v) {
				ri, err := symbolicAntideriv(n.Right, v)
				if err != nil {
					return nil, err
				}
				return BinaryNode{Op: "*", Left: n.Left, Right: ri}, nil
			}
			if !dependsOn(n.Right, v) {
				li, err := symbolicAntideriv(n.Left, v)
				if err != nil {
					return nil, err
				}
				return BinaryNode{Op: "*", Left: n.Right, Right: li}, nil
			}

		case "/":
			if !dependsOn(n.Right, v) {
				li, err := symbolicAntideriv(n.Left, v)
				if err != nil {
					return nil, err
				}
				return BinaryNode{Op: "/", Left: li, Right: n.Right}, nil
			}
			if !dependsOn(n.Left, v) {
				if a, _, ok := extractLinear(n.Right, v); ok && !isNumVal(a, 0) {
					lnPart := CallNode{Name: "ln", Args: []ASTNode{CallNode{Name: "abs", Args: []ASTNode{n.Right}}}}
					return BinaryNode{Op: "/", Left: BinaryNode{Op: "*", Left: n.Left, Right: lnPart}, Right: a}, nil
				}
			}

		case "^":
			if !dependsOn(n.Right, v) {
				if a, _, ok := extractLinear(n.Left, v); ok && !isNumVal(a, 0) {
					if rv, ok := getNum(n.Right); ok && math.Abs(rv+1) <= 1e-14 {
						lnPart := CallNode{Name: "ln", Args: []ASTNode{CallNode{Name: "abs", Args: []ASTNode{n.Left}}}}
						return BinaryNode{Op: "/", Left: lnPart, Right: a}, nil
					}
					newExp := simplifyAST(BinaryNode{Op: "+", Left: n.Right, Right: NumNode{Val: 1}})
					num := BinaryNode{Op: "^", Left: n.Left, Right: newExp}
					den := simplifyAST(BinaryNode{Op: "*", Left: a, Right: newExp})
					return BinaryNode{Op: "/", Left: num, Right: den}, nil
				}
			}
			if !dependsOn(n.Left, v) {
				if a, _, ok := extractLinear(n.Right, v); ok && !isNumVal(a, 0) {
					if id, isE := n.Left.(IdentNode); isE && strings.ToLower(id.Name) == "e" {
						return BinaryNode{Op: "/", Left: n, Right: a}, nil
					}
					den := BinaryNode{Op: "*", Left: a, Right: CallNode{Name: "ln", Args: []ASTNode{n.Left}}}
					return BinaryNode{Op: "/", Left: n, Right: den}, nil
				}
			}
		}

	case CallNode:
		if len(n.Args) == 1 {
			u := n.Args[0]
			if a, _, ok := extractLinear(u, v); ok && !isNumVal(a, 0) {
				fn := strings.ToLower(n.Name)
				var anti ASTNode
				switch fn {
				case "sin":
					anti = UnaryNode{Op: "-", X: CallNode{Name: "cos", Args: []ASTNode{u}}}
				case "cos":
					anti = CallNode{Name: "sin", Args: []ASTNode{u}}
				case "tan":
					anti = UnaryNode{Op: "-", X: CallNode{Name: "ln", Args: []ASTNode{CallNode{Name: "abs", Args: []ASTNode{CallNode{Name: "cos", Args: []ASTNode{u}}}}}}}
				case "exp":
					anti = CallNode{Name: "exp", Args: []ASTNode{u}}
				case "sinh":
					anti = CallNode{Name: "cosh", Args: []ASTNode{u}}
				case "cosh":
					anti = CallNode{Name: "sinh", Args: []ASTNode{u}}
				case "sqrt":
					num := BinaryNode{Op: "*", Left: NumNode{Val: 2}, Right: BinaryNode{Op: "^", Left: u, Right: BinaryNode{Op: "/", Left: NumNode{Val: 3}, Right: NumNode{Val: 2}}}}
					anti = BinaryNode{Op: "/", Left: num, Right: NumNode{Val: 3}}
				case "ln":
					t1 := BinaryNode{Op: "*", Left: u, Right: CallNode{Name: "ln", Args: []ASTNode{u}}}
					anti = BinaryNode{Op: "-", Left: t1, Right: u}
				}
				if anti != nil {
					return BinaryNode{Op: "/", Left: anti, Right: a}, nil
				}
			}
		}
	}

	if poly, ok := extractRealPoly(node, v); ok {
		var terms []ASTNode
		for deg, coeff := range poly {
			if math.Abs(coeff) <= 1e-14 {
				continue
			}
			newDeg := deg + 1
			var xPow ASTNode = IdentNode{Name: v}
			if newDeg > 1 {
				xPow = BinaryNode{Op: "^", Left: xPow, Right: NumNode{Val: float64(newDeg)}}
			}
			term := simplifyAST(BinaryNode{
				Op:    "/",
				Left:  BinaryNode{Op: "*", Left: makeSignedNum(coeff), Right: xPow},
				Right: NumNode{Val: float64(newDeg)},
			})
			terms = append(terms, term)
		}
		if len(terms) == 0 {
			return NumNode{Val: 0}, nil
		}
		acc := terms[len(terms)-1]
		for i := len(terms) - 2; i >= 0; i-- {
			acc = BinaryNode{Op: "+", Left: acc, Right: terms[i]}
		}
		return simplifyAST(acc), nil
	}

	return nil, fmt.Errorf("closed-form antiderivative not found; pass integration bounds int(expr, %s, a, b) for numerical quadrature", v)
}

var gkNodes = [8]float64{
	0.9914553711208126,
	0.9491079123427585,
	0.8648644233597691,
	0.7415311855993945,
	0.5860872354676911,
	0.4058451513773972,
	0.2077849550078985,
	0.0000000000000000,
}

var gkWeights = [8]float64{
	0.02293532201052922,
	0.06309209262997856,
	0.1047900103222502,
	0.1406532597155259,
	0.1690047266392679,
	0.1903505780647854,
	0.2044329400752989,
	0.2094821410847278,
}

var g7Weights = [4]float64{
	0.1294849661688697,
	0.2797053914892766,
	0.3818300505051189,
	0.4179591836734694,
}

func gk15Step(f func(float64) complex128, a, b float64) (complex128, float64) {
	center := 0.5 * (a + b)
	half := 0.5 * (b - a)

	fCenter := f(center)
	resK := fCenter * complex(gkWeights[7], 0)
	resG := fCenter * complex(g7Weights[3], 0)

	for j := 0; j < 7; j++ {
		absc := half * gkNodes[j]
		fSum := f(center-absc) + f(center+absc)
		resK += fSum * complex(gkWeights[j], 0)
		if j%2 == 1 {
			resG += fSum * complex(g7Weights[j/2], 0)
		}
	}

	resK *= complex(half, 0)
	resG *= complex(half, 0)
	errEst := cmplx.Abs(resK - resG)
	return resK, errEst
}

func adaptiveGK(f func(float64) complex128, a, b, tol float64, depth int) complex128 {
	val, errEst := gk15Step(f, a, b)
	if math.IsNaN(real(val)) || math.IsNaN(imag(val)) {
		mid := 0.5 * (a + b)
		if depth >= 16 {
			return 0
		}
		return adaptiveGK(f, a, mid, tol/2, depth+1) + adaptiveGK(f, mid, b, tol/2, depth+1)
	}
	if errEst <= tol || depth >= 18 {
		return val
	}
	mid := 0.5 * (a + b)
	return adaptiveGK(f, a, mid, tol/2, depth+1) + adaptiveGK(f, mid, b, tol/2, depth+1)
}

func integrateComplex(f func(float64) complex128, a, b float64) complex128 {
	if a == b {
		return 0
	}
	sign := complex(1, 0)
	if a > b {
		a, b = b, a
		sign = complex(-1, 0)
	}

	aNegInf := math.IsInf(a, -1)
	bPosInf := math.IsInf(b, 1)

	if aNegInf && bPosInf {
		g := func(t float64) complex128 {
			denom := 1 - t*t
			if math.Abs(denom) < 1e-15 {
				return 0
			}
			x := t / denom
			jac := (1 + t*t) / (denom * denom)
			v := f(x) * complex(jac, 0)
			if math.IsNaN(real(v)) || math.IsInf(real(v), 0) || math.IsNaN(imag(v)) || math.IsInf(imag(v), 0) {
				return 0
			}
			return v
		}
		return cleanComplex(sign * adaptiveGK(g, -1, 1, 1e-11, 0))
	}

	if !aNegInf && bPosInf {
		g := func(t float64) complex128 {
			denom := 1 - t
			if math.Abs(denom) < 1e-15 {
				return 0
			}
			x := a + t/denom
			jac := 1.0 / (denom * denom)
			v := f(x) * complex(jac, 0)
			if math.IsNaN(real(v)) || math.IsInf(real(v), 0) || math.IsNaN(imag(v)) || math.IsInf(imag(v), 0) {
				return 0
			}
			return v
		}
		return cleanComplex(sign * adaptiveGK(g, 0, 1, 1e-11, 0))
	}

	if aNegInf && !bPosInf {
		g := func(t float64) complex128 {
			denom := 1 - t
			if math.Abs(denom) < 1e-15 {
				return 0
			}
			x := b - t/denom
			jac := 1.0 / (denom * denom)
			v := f(x) * complex(jac, 0)
			if math.IsNaN(real(v)) || math.IsInf(real(v), 0) || math.IsNaN(imag(v)) || math.IsInf(imag(v), 0) {
				return 0
			}
			return v
		}
		return cleanComplex(sign * adaptiveGK(g, 0, 1, 1e-11, 0))
	}

	return cleanComplex(sign * adaptiveGK(f, a, b, 1e-12, 0))
}

func extractComplexPoly(node ASTNode, v string, e *Engine, scope map[string]Value) ([]complex128, bool) {
	if !dependsOn(node, v) {
		val, err := e.evalAST(node, scope, 0)
		if err != nil || val.Kind != KindScalar {
			return nil, false
		}
		return []complex128{val.S}, true
	}

	switch n := node.(type) {
	case IdentNode:
		if n.Name == v {
			return []complex128{0, 1}, true
		}
		return nil, false

	case UnaryNode:
		sub, ok := extractComplexPoly(n.X, v, e, scope)
		if !ok {
			return nil, false
		}
		if n.Op == "-" {
			out := make([]complex128, len(sub))
			for i, z := range sub {
				out[i] = -z
			}
			return out, true
		}
		return sub, true

	case BinaryNode:
		switch n.Op {
		case "+", "-":
			lp, ok1 := extractComplexPoly(n.Left, v, e, scope)
			rp, ok2 := extractComplexPoly(n.Right, v, e, scope)
			if !ok1 || !ok2 {
				return nil, false
			}
			maxLen := len(lp)
			if len(rp) > maxLen {
				maxLen = len(rp)
			}
			out := make([]complex128, maxLen)
			copy(out, lp)
			for i, z := range rp {
				if n.Op == "+" {
					out[i] += z
				} else {
					out[i] -= z
				}
			}
			return out, true

		case "*":
			lp, ok1 := extractComplexPoly(n.Left, v, e, scope)
			rp, ok2 := extractComplexPoly(n.Right, v, e, scope)
			if !ok1 || !ok2 {
				return nil, false
			}
			return polyConvolveComplex(lp, rp), true

		case "/":
			if dependsOn(n.Right, v) {
				return nil, false
			}
			denomVal, err := e.evalAST(n.Right, scope, 0)
			if err != nil || denomVal.Kind != KindScalar || denomVal.S == 0 {
				return nil, false
			}
			lp, ok1 := extractComplexPoly(n.Left, v, e, scope)
			if !ok1 {
				return nil, false
			}
			out := make([]complex128, len(lp))
			for i, z := range lp {
				out[i] = z / denomVal.S
			}
			return out, true

		case "^":
			if dependsOn(n.Right, v) {
				return nil, false
			}
			expVal, err := e.evalAST(n.Right, scope, 0)
			if err != nil || !expVal.IsReal() {
				return nil, false
			}
			er := real(expVal.S)
			if er < 0 || er != math.Trunc(er) || er > 64 {
				return nil, false
			}
			deg := int(er)
			bp, ok := extractComplexPoly(n.Left, v, e, scope)
			if !ok {
				return nil, false
			}
			res := []complex128{1}
			base := bp
			for deg > 0 {
				if deg%2 == 1 {
					res = polyConvolveComplex(res, base)
				}
				deg /= 2
				if deg > 0 {
					base = polyConvolveComplex(base, base)
				}
			}
			return res, true
		}
	}

	return nil, false
}

func extractRealPoly(node ASTNode, v string) ([]float64, bool) {
	cp, ok := extractComplexPoly(node, v, defaultEngine, nil)
	if !ok {
		return nil, false
	}
	rp := make([]float64, len(cp))
	for i, z := range cp {
		if math.Abs(imag(z)) > 1e-13 {
			return nil, false
		}
		rp[i] = real(z)
	}
	return rp, true
}

func polyConvolveComplex(a, b []complex128) []complex128 {
	if len(a) == 0 || len(b) == 0 {
		return []complex128{0}
	}
	out := make([]complex128, len(a)+len(b)-1)
	for i, za := range a {
		for j, zb := range b {
			out[i+j] += za * zb
		}
	}
	return out
}

func solvePolynomialAberth(coeffs []complex128) ([]complex128, error) {
	d := len(coeffs) - 1
	for d > 0 && cmplx.Abs(coeffs[d]) <= 1e-14 {
		d--
	}
	if d == 0 {
		return nil, fmt.Errorf("constant non-zero expression has no roots")
	}
	coeffs = coeffs[:d+1]

	var roots []complex128
	zeroShift := 0
	for zeroShift < d && cmplx.Abs(coeffs[zeroShift]) <= 1e-14 {
		roots = append(roots, 0)
		zeroShift++
	}
	coeffs = coeffs[zeroShift:]
	d = len(coeffs) - 1
	if d == 0 {
		return sortComplexSlice(roots), nil
	}

	lead := coeffs[d]
	p := make([]complex128, d+1)
	for i := 0; i <= d; i++ {
		p[i] = coeffs[i] / lead
	}

	if d == 1 {
		roots = append(roots, cleanComplex(-p[0]))
		return sortComplexSlice(roots), nil
	}

	if d == 2 {
		b, c := p[1], p[0]
		disc := cmplx.Sqrt(b*b - 4*c)
		roots = append(roots, cleanComplex((-b-disc)/2), cleanComplex((-b+disc)/2))
		return sortComplexSlice(roots), nil
	}

	radius := 0.0
	for i := 0; i < d; i++ {
		if m := cmplx.Abs(p[i]); m > radius {
			radius = m
		}
	}
	radius = 1.0 + radius

	z := make([]complex128, d)
	for k := 0; k < d; k++ {
		angle := 2*math.Pi*float64(k)/float64(d) + 0.371238898
		r := radius * (0.85 + 0.3*float64(k%3)/2.0)
		z[k] = cmplx.Rect(r, angle)
	}

	for iter := 0; iter < 300; iter++ {
		maxStep := 0.0
		offsets := make([]complex128, d)
		for k := 0; k < d; k++ {
			zk := z[k]
			val := p[d]
			deriv := complex(0, 0)
			for i := d - 1; i >= 0; i-- {
				deriv = deriv*zk + val
				val = val*zk + p[i]
			}
			if cmplx.Abs(val) <= 1e-15 {
				continue
			}
			if cmplx.Abs(deriv) <= 1e-30 {
				offsets[k] = complex(1e-5, 1e-5)
				continue
			}
			newton := val / deriv
			var sumInv complex128
			for j := 0; j < d; j++ {
				if j != k {
					diff := zk - z[j]
					if cmplx.Abs(diff) > 1e-30 {
						sumInv += 1 / diff
					}
				}
			}
			denom := 1 - newton*sumInv
			if cmplx.Abs(denom) <= 1e-30 {
				offsets[k] = newton
			} else {
				offsets[k] = newton / denom
			}
		}
		for k := 0; k < d; k++ {
			z[k] -= offsets[k]
			if s := cmplx.Abs(offsets[k]); s > maxStep {
				maxStep = s
			}
		}
		if maxStep <= 1e-14 {
			break
		}
	}

	for _, r := range z {
		roots = append(roots, cleanComplex(r))
	}
	return sortComplexSlice(roots), nil
}

func sortComplexSlice(vals []complex128) []complex128 {
	out := make([]complex128, len(vals))
	for i, z := range vals {
		out[i] = cleanComplex(z)
	}
	sort.Slice(out, func(i, j int) bool {
		ri, rj := real(out[i]), real(out[j])
		if math.Abs(ri-rj) > 1e-9 {
			return ri < rj
		}
		return imag(out[i]) < imag(out[j])
	})
	return out
}

func brentRoot(f func(float64) float64, a, b float64) (float64, bool) {
	fa, fb := f(a), f(b)
	if math.IsNaN(fa) || math.IsNaN(fb) || math.IsInf(fa, 0) || math.IsInf(fb, 0) {
		return 0, false
	}
	if math.Abs(fa) <= 1e-14 {
		return a, true
	}
	if math.Abs(fb) <= 1e-14 {
		return b, true
	}
	if (fa > 0 && fb > 0) || (fa < 0 && fb < 0) {
		return 0, false
	}
	if math.Abs(fa) < math.Abs(fb) {
		a, b = b, a
		fa, fb = fb, fa
	}
	c, fc := a, fa
	mflag := true
	var d float64

	for iter := 0; iter < 120; iter++ {
		if math.Abs(fb) <= 1e-14 || math.Abs(b-a) <= 1e-13*(1+math.Abs(b)) {
			if math.Abs(f(b)) > 1e-4 {
				return 0, false
			}
			return b, true
		}
		var s float64
		if fa != fc && fb != fc {
			s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb))
		} else {
			s = b - fb*(b-a)/(fb-fa)
		}
		cond1 := !((s >= (3*a+b)/4 && s <= b) || (s <= (3*a+b)/4 && s >= b))
		cond2 := mflag && math.Abs(s-b) >= math.Abs(b-c)/2
		cond3 := !mflag && math.Abs(s-b) >= math.Abs(c-d)/2
		cond4 := mflag && math.Abs(b-c) < 1e-14
		cond5 := !mflag && math.Abs(c-d) < 1e-14
		if cond1 || cond2 || cond3 || cond4 || cond5 {
			s = 0.5 * (a + b)
			mflag = true
		} else {
			mflag = false
		}
		fs := f(s)
		if math.IsNaN(fs) || math.IsInf(fs, 0) {
			return 0, false
		}
		d = c
		c, fc = b, fb
		if (fa > 0 && fs < 0) || (fa < 0 && fs > 0) {
			b, fb = s, fs
		} else {
			a, fa = s, fs
		}
		if math.Abs(fa) < math.Abs(fb) {
			a, b = b, a
			fa, fb = fb, fa
		}
	}
	if math.Abs(f(b)) <= 1e-6 {
		return b, true
	}
	return 0, false
}

func solveTranscendentalReal(f func(float64) float64, xMin, xMax float64) []complex128 {
	const samples = 800
	step := (xMax - xMin) / float64(samples)
	var found []float64

	addUnique := func(r float64) {
		if math.Abs(r) < 1e-13 {
			r = 0
		}
		for _, ex := range found {
			if math.Abs(ex-r) <= 1e-7*(1+math.Abs(r)) {
				return
			}
		}
		found = append(found, r)
	}

	prevX := xMin
	prevY := f(prevX)
	for i := 1; i <= samples; i++ {
		curX := xMin + float64(i)*step
		curY := f(curX)
		if !math.IsNaN(prevY) && math.Abs(prevY) <= 1e-13 {
			addUnique(prevX)
		} else if !math.IsNaN(prevY) && !math.IsNaN(curY) && !math.IsInf(prevY, 0) && !math.IsInf(curY, 0) {
			if (prevY > 0 && curY < 0) || (prevY < 0 && curY > 0) {
				if root, ok := brentRoot(f, prevX, curX); ok {
					addUnique(root)
				}
			} else if math.Abs(curY) < 1e-2 && i < samples {
				nextX := xMin + float64(i+1)*step
				nextY := f(nextX)
				if !math.IsNaN(nextY) && ((curY >= 0 && prevY > curY && nextY > curY) || (curY <= 0 && prevY < curY && nextY < curY)) {
					x := curX
					for k := 0; k < 30; k++ {
						fx := f(x)
						if math.Abs(fx) <= 1e-13 {
							addUnique(x)
							break
						}
						h := 1e-7 * (1 + math.Abs(x))
						df := (f(x+h) - f(x-h)) / (2 * h)
						if math.Abs(df) < 1e-14 {
							break
						}
						x -= fx / df
					}
				}
			}
		}
		prevX, prevY = curX, curY
	}
	if !math.IsNaN(prevY) && math.Abs(prevY) <= 1e-13 {
		addUnique(prevX)
	}

	sort.Float64s(found)
	out := make([]complex128, len(found))
	for i, r := range found {
		out[i] = complex(r, 0)
	}
	return out
}

func oneSidedLimit(f func(float64) complex128, target float64, dir float64) (complex128, bool) {
	const m = 9
	seq := make([]complex128, 0, m)
	for k := 0; k < m; k++ {
		var x float64
		if math.IsInf(target, 1) {
			x = math.Pow(4, float64(k+3))
		} else if math.IsInf(target, -1) {
			x = -math.Pow(4, float64(k+3))
		} else {
			h := 0.05 * math.Pow(0.25, float64(k))
			x = target + dir*h
		}
		y := f(x)
		if math.IsNaN(real(y)) || math.IsNaN(imag(y)) {
			return 0, false
		}
		seq = append(seq, y)
	}

	mags := make([]float64, len(seq))
	for i, z := range seq {
		mags[i] = cmplx.Abs(z)
	}
	if math.IsInf(mags[len(mags)-1], 0) || (mags[len(mags)-1] > 1e7 && mags[len(mags)-1] > 3*mags[len(mags)-2] && mags[len(mags)-2] > 3*mags[len(mags)-3]) {
		last := seq[len(seq)-1]
		if math.Abs(imag(last)) <= 1e-8*math.Abs(real(last)) {
			if real(last) > 0 {
				return complex(math.Inf(1), 0), true
			}
			return complex(math.Inf(-1), 0), true
		}
		return complex(math.Inf(1), 0), true
	}

	table := make([][]complex128, m)
	for i := 0; i < m; i++ {
		table[i] = make([]complex128, i+1)
		table[i][0] = seq[i]
		factor := 4.0
		for j := 1; j <= i; j++ {
			table[i][j] = (complex(factor, 0)*table[i][j-1] - table[i-1][j-1]) / complex(factor-1, 0)
			factor *= 4.0
		}
	}

	best := seq[m-1]
	bestErr := cmplx.Abs(seq[m-1] - seq[m-2])
	for i := 2; i < m; i++ {
		for j := 1; j <= i && j <= 4; j++ {
			cand := table[i][j]
			prev := table[i-1][j-1]
			if math.IsNaN(real(cand)) || math.IsInf(real(cand), 0) {
				continue
			}
			err := cmplx.Abs(cand - prev)
			if err < bestErr {
				bestErr = err
				best = cand
			}
		}
	}

	if bestErr > 1e-3*math.Max(1.0, cmplx.Abs(best)) {
		return 0, false
	}
	if cmplx.Abs(best) < 1e-11 && cmplx.Abs(seq[m-1]) < 1e-8 {
		best = 0
	}
	return cleanComplex(best), true
}

func extractIdentList(node ASTNode) ([]string, error) {
	if id, ok := node.(IdentNode); ok {
		return []string{id.Name}, nil
	}
	if vec, ok := node.(VectorNode); ok {
		names := make([]string, len(vec.Elems))
		for i, el := range vec.Elems {
			id, ok := el.(IdentNode)
			if !ok {
				return nil, fmt.Errorf("variable list element %d must be identifier", i)
			}
			names[i] = id.Name
		}
		return names, nil
	}
	return nil, fmt.Errorf("expected variable or [x, y, ...] variable list")
}

func (e *Engine) evalSpecialForm(n CallNode, scope map[string]Value, depth int) (bool, Value, error) {
	fn := strings.ToLower(n.Name)

	switch fn {
	case "simplify", "simp":
		if len(n.Args) != 1 {
			return true, Value{}, fmt.Errorf("simplify(expr) takes 1 argument")
		}
		expanded := e.expandSymbolicRefs(n.Args[0])
		simp := simplifyAST(expanded)
		return true, SymbolicVal(simp.Format(), simp), nil

	case "expand":
		if len(n.Args) < 1 || len(n.Args) > 2 {
			return true, Value{}, fmt.Errorf("expand(expr [, x]) takes 1 or 2 arguments")
		}
		varName := "x"
		if len(n.Args) == 2 {
			id, ok := n.Args[1].(IdentNode)
			if !ok {
				return true, Value{}, fmt.Errorf("second argument to expand must be variable name")
			}
			varName = id.Name
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		poly, ok := extractRealPoly(expr, varName)
		if !ok {
			simp := simplifyAST(expr)
			return true, SymbolicVal(simp.Format(), simp), nil
		}
		polyNode := buildPolyAST(poly, varName)
		return true, SymbolicVal(polyNode.Format(), polyNode), nil

	case "at", "subs":
		if len(n.Args) != 3 {
			return true, Value{}, fmt.Errorf("at(expr, x, value) takes 3 arguments")
		}
		id, ok := n.Args[1].(IdentNode)
		if !ok {
			return true, Value{}, fmt.Errorf("second argument to at must be variable name")
		}
		targetVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		local := copyScope(scope)
		local[id.Name] = targetVal
		expr := e.expandSymbolicRefs(n.Args[0])
		res, err := e.evalAST(expr, local, depth+1)
		return true, res, err

	case "diff", "deriv":
		if len(n.Args) < 1 || len(n.Args) > 4 {
			return true, Value{}, fmt.Errorf("diff(expr [, x [, order [, point]]]) takes 1 to 4 arguments")
		}
		varName := "x"
		if len(n.Args) >= 2 {
			id, ok := n.Args[1].(IdentNode)
			if !ok {
				return true, Value{}, fmt.Errorf("second argument to diff must be variable name")
			}
			varName = id.Name
		}
		order := 1
		if len(n.Args) >= 3 {
			ordVal, err := e.evalAST(n.Args[2], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			or, err := ordVal.Real()
			if err != nil || or < 0 || or != math.Trunc(or) {
				return true, Value{}, fmt.Errorf("derivative order must be non-negative integer")
			}
			order = int(or)
		}
		cur := e.expandSymbolicRefs(n.Args[0])
		var err error
		for i := 0; i < order; i++ {
			cur, err = symbolicDiff(cur, varName)
			if err != nil {
				return true, Value{}, err
			}
			cur = simplifyAST(cur)
		}
		cur = simplifyAST(cur)
		if len(n.Args) == 4 {
			ptVal, err := e.evalAST(n.Args[3], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			local := copyScope(scope)
			local[varName] = ptVal
			val, err := e.evalAST(cur, local, depth+1)
			return true, val, err
		}
		return true, SymbolicVal(cur.Format(), cur), nil

	case "ndiff":
		if len(n.Args) < 3 || len(n.Args) > 4 {
			return true, Value{}, fmt.Errorf("ndiff(expr, x, point [, order]) takes 3 or 4 arguments")
		}
		id, ok := n.Args[1].(IdentNode)
		if !ok {
			return true, Value{}, fmt.Errorf("second argument to ndiff must be variable name")
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		pt, err := ptVal.Real()
		if err != nil {
			return true, Value{}, err
		}
		order := 1
		if len(n.Args) == 4 {
			ordVal, err := e.evalAST(n.Args[3], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			or, err := ordVal.Real()
			if err != nil || or < 1 || or > 4 || or != math.Trunc(or) {
				return true, Value{}, fmt.Errorf("ndiff order must be 1, 2, 3, or 4")
			}
			order = int(or)
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		evalAt := func(x float64) complex128 {
			local := copyScope(scope)
			local[id.Name] = ScalarReal(x)
			v, err := e.evalAST(expr, local, depth+1)
			if err != nil || v.Kind != KindScalar {
				return complex(math.NaN(), 0)
			}
			return v.S
		}
		h := 1e-4 * math.Max(1.0, math.Abs(pt))
		var res complex128
		switch order {
		case 1:
			d1 := (-evalAt(pt+2*h) + 8*evalAt(pt+h) - 8*evalAt(pt-h) + evalAt(pt-2*h)) / complex(12*h, 0)
			h2 := h / 2
			d2 := (-evalAt(pt+2*h2) + 8*evalAt(pt+h2) - 8*evalAt(pt-h2) + evalAt(pt-2*h2)) / complex(12*h2, 0)
			res = (16*d2 - d1) / 15
		case 2:
			h = 1e-3 * math.Max(1.0, math.Abs(pt))
			res = (-evalAt(pt+2*h) + 16*evalAt(pt+h) - 30*evalAt(pt) + 16*evalAt(pt-h) - evalAt(pt-2*h)) / complex(12*h*h, 0)
		case 3:
			h = 5e-3 * math.Max(1.0, math.Abs(pt))
			res = (evalAt(pt+2*h) - 2*evalAt(pt+h) + 2*evalAt(pt-h) - evalAt(pt-2*h)) / complex(2*h*h*h, 0)
		case 4:
			h = 1e-2 * math.Max(1.0, math.Abs(pt))
			res = (evalAt(pt+2*h) - 4*evalAt(pt+h) + 6*evalAt(pt) - 4*evalAt(pt-h) + evalAt(pt-2*h)) / complex(h*h*h*h, 0)
		}
		return true, ScalarComplex(res), nil

	case "taylor", "maclaurin":
		varName := "x"
		center := 0.0
		order := 5
		if fn == "maclaurin" {
			if len(n.Args) < 1 || len(n.Args) > 3 {
				return true, Value{}, fmt.Errorf("maclaurin(expr [, x [, order]]) takes 1 to 3 arguments")
			}
			if len(n.Args) >= 2 {
				id, ok := n.Args[1].(IdentNode)
				if !ok {
					return true, Value{}, fmt.Errorf("second argument must be variable name")
				}
				varName = id.Name
			}
			if len(n.Args) == 3 {
				ov, err := e.evalAST(n.Args[2], scope, depth+1)
				if err != nil {
					return true, Value{}, err
				}
				or, err := ov.Real()
				if err != nil {
					return true, Value{}, err
				}
				order = int(math.Round(or))
			}
		} else {
			if len(n.Args) < 1 || len(n.Args) > 4 {
				return true, Value{}, fmt.Errorf("taylor(expr [, x [, center [, order]]]) takes 1 to 4 arguments")
			}
			if len(n.Args) >= 2 {
				id, ok := n.Args[1].(IdentNode)
				if !ok {
					return true, Value{}, fmt.Errorf("second argument must be variable name")
				}
				varName = id.Name
			}
			if len(n.Args) >= 3 {
				cv, err := e.evalAST(n.Args[2], scope, depth+1)
				if err != nil {
					return true, Value{}, err
				}
				cr, err := cv.Real()
				if err != nil {
					return true, Value{}, err
				}
				center = cr
			}
			if len(n.Args) == 4 {
				ov, err := e.evalAST(n.Args[3], scope, depth+1)
				if err != nil {
					return true, Value{}, err
				}
				or, err := ov.Real()
				if err != nil {
					return true, Value{}, err
				}
				order = int(math.Round(or))
			}
		}
		if order < 0 || order > 20 {
			return true, Value{}, fmt.Errorf("taylor order must be between 0 and 20")
		}

		cur := e.expandSymbolicRefs(n.Args[0])
		var terms []ASTNode
		for k := 0; k <= order; k++ {
			if k > 0 {
				var err error
				cur, err = symbolicDiff(cur, varName)
				if err != nil {
					return true, Value{}, err
				}
				cur = simplifyAST(cur)
			}
			local := copyScope(scope)
			local[varName] = ScalarReal(center)
			dkVal, err := e.evalAST(cur, local, depth+1)
			if err != nil {
				return true, Value{}, fmt.Errorf("cannot evaluate derivative %d at %v", k, center)
			}
			dk, err := dkVal.Real()
			if err != nil {
				return true, Value{}, err
			}
			if math.Abs(dk) <= 1e-13 {
				continue
			}
			factK := factorialReal(float64(k))
			var coeffNode ASTNode
			if math.Abs(dk-math.Round(dk)) <= 1e-11 && factK < 1e12 {
				numI := int64(math.Round(dk))
				denI := int64(factK)
				g := gcdInt64(numI, denI)
				numI /= g
				denI /= g
				if denI == 1 {
					coeffNode = makeSignedNum(float64(numI))
				} else {
					coeffNode = BinaryNode{Op: "/", Left: makeSignedNum(float64(numI)), Right: NumNode{Val: float64(denI)}}
				}
			} else {
				coeffNode = makeSignedNum(dk / factK)
			}

			if k == 0 {
				terms = append(terms, coeffNode)
				continue
			}
			var base ASTNode = IdentNode{Name: varName}
			if math.Abs(center) > 1e-14 {
				if center > 0 {
					base = BinaryNode{Op: "-", Left: IdentNode{Name: varName}, Right: NumNode{Val: center}}
				} else {
					base = BinaryNode{Op: "+", Left: IdentNode{Name: varName}, Right: NumNode{Val: -center}}
				}
			}
			var powNode ASTNode = base
			if k > 1 {
				powNode = BinaryNode{Op: "^", Left: base, Right: NumNode{Val: float64(k)}}
			}
			term := simplifyAST(BinaryNode{Op: "*", Left: coeffNode, Right: powNode})
			terms = append(terms, term)
		}

		if len(terms) == 0 {
			return true, SymbolicVal("0", NumNode{Val: 0}), nil
		}
		acc := terms[0]
		for i := 1; i < len(terms); i++ {
			acc = simplifyAST(BinaryNode{Op: "+", Left: acc, Right: terms[i]})
		}
		return true, SymbolicVal(acc.Format(), acc), nil

	case "int", "integrate", "antideriv":
		if len(n.Args) == 1 || len(n.Args) == 2 {
			varName := "x"
			if len(n.Args) == 2 {
				id, ok := n.Args[1].(IdentNode)
				if !ok {
					return true, Value{}, fmt.Errorf("second argument to int must be variable name")
				}
				varName = id.Name
			}
			expr := e.expandSymbolicRefs(n.Args[0])
			anti, err := symbolicAntideriv(expr, varName)
			if err != nil {
				return true, Value{}, err
			}
			anti = simplifyAST(anti)
			return true, SymbolicVal(anti.Format()+" + C", anti), nil
		}
		if len(n.Args) == 4 && fn != "antideriv" {
			id, ok := n.Args[1].(IdentNode)
			if !ok {
				return true, Value{}, fmt.Errorf("second argument to int must be variable name")
			}
			aVal, err := e.evalAST(n.Args[2], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			bVal, err := e.evalAST(n.Args[3], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			a, err1 := aVal.Real()
			b, err2 := bVal.Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("integration bounds must be real")
			}
			expr := e.expandSymbolicRefs(n.Args[0])
			res := integrateComplex(func(x float64) complex128 {
				local := copyScope(scope)
				local[id.Name] = ScalarReal(x)
				v, err := e.evalAST(expr, local, depth+1)
				if err != nil || v.Kind != KindScalar {
					return complex(math.NaN(), 0)
				}
				return v.S
			}, a, b)
			return true, ScalarComplex(res), nil
		}
		return true, Value{}, fmt.Errorf("usage: int(expr [, x]) for symbolic or int(expr, x, a, b) for definite integral")

	case "lim", "limit":
		if len(n.Args) < 3 || len(n.Args) > 4 {
			return true, Value{}, fmt.Errorf("lim(expr, x, target [, dir]) takes 3 or 4 arguments")
		}
		id, ok := n.Args[1].(IdentNode)
		if !ok {
			return true, Value{}, fmt.Errorf("second argument to lim must be variable name")
		}
		tVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		target, err := tVal.Real()
		if err != nil {
			return true, Value{}, fmt.Errorf("limit target must be real or +/-inf")
		}
		dir := 0.0
		if len(n.Args) == 4 {
			dv, err := e.evalAST(n.Args[3], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			dr, err := dv.Real()
			if err != nil {
				return true, Value{}, err
			}
			dir = dr
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		evalF := func(x float64) complex128 {
			local := copyScope(scope)
			local[id.Name] = ScalarReal(x)
			v, err := e.evalAST(expr, local, depth+1)
			if err != nil || v.Kind != KindScalar {
				return complex(math.NaN(), 0)
			}
			return v.S
		}

		if math.IsInf(target, 0) {
			res, ok := oneSidedLimit(evalF, target, 1)
			if !ok {
				return true, SymbolicVal("DNE", nil), nil
			}
			return true, ScalarComplex(res), nil
		}

		if dir > 0 {
			res, ok := oneSidedLimit(evalF, target, 1)
			if !ok {
				return true, SymbolicVal("DNE", nil), nil
			}
			return true, ScalarComplex(res), nil
		}
		if dir < 0 {
			res, ok := oneSidedLimit(evalF, target, -1)
			if !ok {
				return true, SymbolicVal("DNE", nil), nil
			}
			return true, ScalarComplex(res), nil
		}

		rightLim, okR := oneSidedLimit(evalF, target, 1)
		leftLim, okL := oneSidedLimit(evalF, target, -1)
		if !okR || !okL {
			return true, SymbolicVal("DNE", nil), nil
		}
		if math.IsInf(real(rightLim), 0) || math.IsInf(real(leftLim), 0) {
			if math.IsInf(real(rightLim), 1) && math.IsInf(real(leftLim), 1) {
				return true, ScalarReal(math.Inf(1)), nil
			}
			if math.IsInf(real(rightLim), -1) && math.IsInf(real(leftLim), -1) {
				return true, ScalarReal(math.Inf(-1)), nil
			}
			return true, SymbolicVal("DNE", nil), nil
		}
		if cmplx.Abs(rightLim-leftLim) > 1e-4*math.Max(1.0, cmplx.Abs(rightLim)) {
			return true, SymbolicVal("DNE", nil), nil
		}
		return true, ScalarComplex(0.5 * (rightLim + leftLim)), nil

	case "sum", "prod":
		if len(n.Args) == 1 || (len(n.Args) > 1 && n.Args[0] != nil) {
			isIterForm := false
			if len(n.Args) == 4 {
				if _, ok := n.Args[1].(IdentNode); ok {
					isIterForm = true
				}
			}
			if !isIterForm {
				return false, Value{}, nil
			}
		}
		id := n.Args[1].(IdentNode)
		sVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		eVal, err := e.evalAST(n.Args[3], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		startF, err1 := sVal.Real()
		endF, err2 := eVal.Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("bounds must be real")
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		start := int64(math.Round(startF))

		if fn == "sum" && math.IsInf(endF, 1) {
			res, err := evalInfiniteSum(e, expr, id.Name, start, scope, depth+1)
			return true, res, err
		}

		end := int64(math.Round(endF))
		if end-start > 1_000_000 {
			return true, Value{}, fmt.Errorf("iteration count exceeds 1,000,000")
		}
		local := copyScope(scope)
		if fn == "sum" {
			var acc complex128
			for k := start; k <= end; k++ {
				local[id.Name] = ScalarReal(float64(k))
				v, err := e.evalAST(expr, local, depth+1)
				if err != nil {
					return true, Value{}, err
				}
				if v.Kind != KindScalar {
					return true, Value{}, fmt.Errorf("summand must be scalar")
				}
				acc += v.S
			}
			return true, ScalarComplex(acc), nil
		}
		acc := complex(1, 0)
		for k := start; k <= end; k++ {
			local[id.Name] = ScalarReal(float64(k))
			v, err := e.evalAST(expr, local, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			if v.Kind != KindScalar {
				return true, Value{}, fmt.Errorf("product term must be scalar")
			}
			acc *= v.S
		}
		return true, ScalarComplex(acc), nil

	case "solve", "roots":
		if len(n.Args) < 1 || len(n.Args) > 4 {
			return true, Value{}, fmt.Errorf("solve(expr_or_eq [, x [, min, max]]) takes 1 to 4 arguments")
		}
		varName := "x"
		if len(n.Args) >= 2 {
			id, ok := n.Args[1].(IdentNode)
			if !ok {
				return true, Value{}, fmt.Errorf("second argument to solve must be variable name")
			}
			varName = id.Name
		}

		rawExpr := e.expandSymbolicRefs(n.Args[0])
		if bin, ok := rawExpr.(BinaryNode); ok && (bin.Op == "=" || bin.Op == "==") {
			rawExpr = BinaryNode{Op: "-", Left: bin.Left, Right: bin.Right}
		}

		if len(n.Args) <= 2 {
			if polyCoeffs, ok := extractComplexPoly(rawExpr, varName, e, scope); ok {
				roots, err := solvePolynomialAberth(polyCoeffs)
				if err != nil {
					return true, Value{}, err
				}
				if len(roots) == 1 {
					return true, ScalarComplex(roots[0]), nil
				}
				return true, VectorVal(roots), nil
			}
		}

		xMin, xMax := -50.0, 50.0
		if len(n.Args) == 3 {
			gv, err := e.evalAST(n.Args[2], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			guess, err := gv.Real()
			if err != nil {
				return true, Value{}, err
			}
			xMin, xMax = guess-10, guess+10
		} else if len(n.Args) == 4 {
			v1, err1 := e.evalAST(n.Args[2], scope, depth+1)
			v2, err2 := e.evalAST(n.Args[3], scope, depth+1)
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("invalid range bounds")
			}
			r1, err1 := v1.Real()
			r2, err2 := v2.Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("range bounds must be real")
			}
			xMin, xMax = r1, r2
		}

		evalF := func(x float64) float64 {
			local := copyScope(scope)
			local[varName] = ScalarReal(x)
			v, err := e.evalAST(rawExpr, local, depth+1)
			if err != nil || !v.IsReal() {
				return math.NaN()
			}
			return real(v.S)
		}

		roots := solveTranscendentalReal(evalF, xMin, xMax)
		if len(roots) == 0 {
			return true, VectorVal(nil), nil
		}
		if len(roots) == 1 {
			return true, ScalarComplex(roots[0]), nil
		}
		return true, VectorVal(roots), nil

	case "rk4":
		if len(n.Args) < 6 || len(n.Args) > 7 {
			return true, Value{}, fmt.Errorf("rk4(dydx, x, y, x0, y0, xEnd [, steps]) takes 6 or 7 arguments")
		}
		xId, ok1 := n.Args[1].(IdentNode)
		yId, ok2 := n.Args[2].(IdentNode)
		if !ok1 || !ok2 {
			return true, Value{}, fmt.Errorf("arguments 2 and 3 to rk4 must be variable names x, y")
		}
		x0Val, err := e.evalAST(n.Args[3], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		y0Val, err := e.evalAST(n.Args[4], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		xEndVal, err := e.evalAST(n.Args[5], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		x0, err1 := x0Val.Real()
		y0, err2 := y0Val.Real()
		xEnd, err3 := xEndVal.Real()
		if err1 != nil || err2 != nil || err3 != nil {
			return true, Value{}, fmt.Errorf("rk4 initial values and endpoint must be real")
		}
		steps := 1000
		if len(n.Args) == 7 {
			sv, err := e.evalAST(n.Args[6], scope, depth+1)
			if err != nil {
				return true, Value{}, err
			}
			sr, err := sv.Real()
			if err != nil || sr < 1 || sr > 500000 {
				return true, Value{}, fmt.Errorf("steps must be between 1 and 500000")
			}
			steps = int(math.Round(sr))
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		local := copyScope(scope)
		f := func(xv, yv float64) (float64, error) {
			local[xId.Name] = ScalarReal(xv)
			local[yId.Name] = ScalarReal(yv)
			val, err := e.evalAST(expr, local, depth+1)
			if err != nil {
				return 0, err
			}
			return val.Real()
		}
		h := (xEnd - x0) / float64(steps)
		x, y := x0, y0
		for i := 0; i < steps; i++ {
			k1, err := f(x, y)
			if err != nil {
				return true, Value{}, err
			}
			k2, err := f(x+0.5*h, y+0.5*h*k1)
			if err != nil {
				return true, Value{}, err
			}
			k3, err := f(x+0.5*h, y+0.5*h*k2)
			if err != nil {
				return true, Value{}, err
			}
			k4, err := f(x+h, y+h*k3)
			if err != nil {
				return true, Value{}, err
			}
			y += (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
			x += h
		}
		return true, ScalarReal(y), nil

	case "grad", "gradient":
		if len(n.Args) < 2 || len(n.Args) > 3 {
			return true, Value{}, fmt.Errorf("grad(expr, [x, y, ...] [, [p1, p2, ...]]) takes 2 or 3 arguments")
		}
		vars, err := extractIdentList(n.Args[1])
		if err != nil {
			return true, Value{}, err
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		partials := make([]ASTNode, len(vars))
		for i, v := range vars {
			d, err := symbolicDiff(expr, v)
			if err != nil {
				return true, Value{}, err
			}
			partials[i] = simplifyAST(d)
		}
		if len(n.Args) == 2 {
			parts := make([]string, len(partials))
			for i, p := range partials {
				parts[i] = p.Format()
			}
			return true, SymbolicVal("["+strings.Join(parts, ", ")+"]", VectorNode{Elems: partials}), nil
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		if ptVal.Kind != KindVector || len(ptVal.V) != len(vars) {
			return true, Value{}, fmt.Errorf("point must be vector of length %d", len(vars))
		}
		local := copyScope(scope)
		for i, v := range vars {
			local[v] = ScalarComplex(ptVal.V[i])
		}
		out := make([]complex128, len(vars))
		for i, p := range partials {
			val, err := e.evalAST(p, local, depth+1)
			if err != nil || val.Kind != KindScalar {
				return true, Value{}, fmt.Errorf("failed evaluating partial derivative %d", i)
			}
			out[i] = val.S
		}
		return true, VectorVal(out), nil

	case "hessian":
		if len(n.Args) != 3 {
			return true, Value{}, fmt.Errorf("hessian(expr, [x, y, ...], [p1, p2, ...]) takes 3 arguments")
		}
		vars, err := extractIdentList(n.Args[1])
		if err != nil {
			return true, Value{}, err
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		if ptVal.Kind != KindVector || len(ptVal.V) != len(vars) {
			return true, Value{}, fmt.Errorf("point must be vector of length %d", len(vars))
		}
		local := copyScope(scope)
		for i, v := range vars {
			local[v] = ScalarComplex(ptVal.V[i])
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		k := len(vars)
		mat := make([][]complex128, k)
		for i := 0; i < k; i++ {
			mat[i] = make([]complex128, k)
			di, err := symbolicDiff(expr, vars[i])
			if err != nil {
				return true, Value{}, err
			}
			di = simplifyAST(di)
			for j := 0; j < k; j++ {
				dij, err := symbolicDiff(di, vars[j])
				if err != nil {
					return true, Value{}, err
				}
				val, err := e.evalAST(simplifyAST(dij), local, depth+1)
				if err != nil || val.Kind != KindScalar {
					return true, Value{}, fmt.Errorf("failed evaluating Hessian [%d,%d]", i, j)
				}
				mat[i][j] = val.S
			}
		}
		return true, MatrixVal(mat), nil

	case "divergence", "div":
		if len(n.Args) < 2 || len(n.Args) > 3 {
			return true, Value{}, fmt.Errorf("divergence([Fx, Fy, ...], [x, y, ...] [, point]) takes 2 or 3 arguments")
		}
		fieldVec, ok := n.Args[0].(VectorNode)
		if !ok {
			return true, Value{}, fmt.Errorf("first argument to divergence must be vector field [Fx, Fy, ...]")
		}
		vars, err := extractIdentList(n.Args[1])
		if err != nil {
			return true, Value{}, err
		}
		if len(fieldVec.Elems) != len(vars) {
			return true, Value{}, fmt.Errorf("field dimension (%d) must match variables (%d)", len(fieldVec.Elems), len(vars))
		}
		var sumAST ASTNode = NumNode{Val: 0}
		for i, v := range vars {
			d, err := symbolicDiff(e.expandSymbolicRefs(fieldVec.Elems[i]), v)
			if err != nil {
				return true, Value{}, err
			}
			sumAST = simplifyAST(BinaryNode{Op: "+", Left: sumAST, Right: d})
		}
		if len(n.Args) == 2 {
			return true, SymbolicVal(sumAST.Format(), sumAST), nil
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		if ptVal.Kind != KindVector || len(ptVal.V) != len(vars) {
			return true, Value{}, fmt.Errorf("point must be vector of length %d", len(vars))
		}
		local := copyScope(scope)
		for i, v := range vars {
			local[v] = ScalarComplex(ptVal.V[i])
		}
		res, err := e.evalAST(sumAST, local, depth+1)
		return true, res, err

	case "curl":
		if len(n.Args) < 2 || len(n.Args) > 3 {
			return true, Value{}, fmt.Errorf("curl([Fx, Fy, Fz], [x, y, z] [, point]) takes 2 or 3 arguments")
		}
		fieldVec, ok := n.Args[0].(VectorNode)
		if !ok || len(fieldVec.Elems) != 3 {
			return true, Value{}, fmt.Errorf("curl requires 3D vector field [Fx, Fy, Fz]")
		}
		vars, err := extractIdentList(n.Args[1])
		if err != nil || len(vars) != 3 {
			return true, Value{}, fmt.Errorf("curl requires 3 variables [x, y, z]")
		}
		fx := e.expandSymbolicRefs(fieldVec.Elems[0])
		fy := e.expandSymbolicRefs(fieldVec.Elems[1])
		fz := e.expandSymbolicRefs(fieldVec.Elems[2])

		dfzDy, err1 := symbolicDiff(fz, vars[1])
		dfyDz, err2 := symbolicDiff(fy, vars[2])
		dfxDz, err3 := symbolicDiff(fx, vars[2])
		dfzDx, err4 := symbolicDiff(fz, vars[0])
		dfyDx, err5 := symbolicDiff(fy, vars[0])
		dfxDy, err6 := symbolicDiff(fx, vars[1])
		if err1 != nil || err2 != nil || err3 != nil || err4 != nil || err5 != nil || err6 != nil {
			return true, Value{}, fmt.Errorf("cannot differentiate vector field")
		}
		cx := simplifyAST(BinaryNode{Op: "-", Left: dfzDy, Right: dfyDz})
		cy := simplifyAST(BinaryNode{Op: "-", Left: dfxDz, Right: dfzDx})
		cz := simplifyAST(BinaryNode{Op: "-", Left: dfyDx, Right: dfxDy})

		if len(n.Args) == 2 {
			sym := fmt.Sprintf("[%s, %s, %s]", cx.Format(), cy.Format(), cz.Format())
			return true, SymbolicVal(sym, VectorNode{Elems: []ASTNode{cx, cy, cz}}), nil
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		if ptVal.Kind != KindVector || len(ptVal.V) != 3 {
			return true, Value{}, fmt.Errorf("point must be 3D vector")
		}
		local := copyScope(scope)
		for i, v := range vars {
			local[v] = ScalarComplex(ptVal.V[i])
		}
		vx, err1 := e.evalAST(cx, local, depth+1)
		vy, err2 := e.evalAST(cy, local, depth+1)
		vz, err3 := e.evalAST(cz, local, depth+1)
		if err1 != nil || err2 != nil || err3 != nil {
			return true, Value{}, fmt.Errorf("failed evaluating curl at point")
		}
		return true, VectorVal([]complex128{vx.S, vy.S, vz.S}), nil

	case "laplacian":
		if len(n.Args) < 2 || len(n.Args) > 3 {
			return true, Value{}, fmt.Errorf("laplacian(expr, [x, y, ...] [, point]) takes 2 or 3 arguments")
		}
		vars, err := extractIdentList(n.Args[1])
		if err != nil {
			return true, Value{}, err
		}
		expr := e.expandSymbolicRefs(n.Args[0])
		var sumAST ASTNode = NumNode{Val: 0}
		for _, v := range vars {
			d1, err := symbolicDiff(expr, v)
			if err != nil {
				return true, Value{}, err
			}
			d2, err := symbolicDiff(simplifyAST(d1), v)
			if err != nil {
				return true, Value{}, err
			}
			sumAST = simplifyAST(BinaryNode{Op: "+", Left: sumAST, Right: d2})
		}
		if len(n.Args) == 2 {
			return true, SymbolicVal(sumAST.Format(), sumAST), nil
		}
		ptVal, err := e.evalAST(n.Args[2], scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		if ptVal.Kind != KindVector || len(ptVal.V) != len(vars) {
			return true, Value{}, fmt.Errorf("point must be vector of length %d", len(vars))
		}
		local := copyScope(scope)
		for i, v := range vars {
			local[v] = ScalarComplex(ptVal.V[i])
		}
		res, err := e.evalAST(sumAST, local, depth+1)
		return true, res, err
	}

	return false, Value{}, nil
}

func evalInfiniteSum(e *Engine, expr ASTNode, varName string, start int64, scope map[string]Value, depth int) (Value, error) {
	local := copyScope(scope)
	termAt := func(n float64) (complex128, error) {
		local[varName] = ScalarReal(n)
		v, err := e.evalAST(expr, local, depth+1)
		if err != nil || v.Kind != KindScalar {
			return 0, fmt.Errorf("invalid series term at %v", n)
		}
		return v.S, nil
	}

	const directN = 500
	var partial complex128
	for k := int64(0); k < directN; k++ {
		t, err := termAt(float64(start + k))
		if err != nil {
			return Value{}, err
		}
		partial += t
		if k > 20 && cmplx.Abs(t) <= 1e-16*cmplx.Abs(partial) {
			return ScalarComplex(partial), nil
		}
	}

	tTail, err := termAt(float64(start + directN))
	if err == nil && cmplx.Abs(tTail) > 1e-13 {
		tPrev, _ := termAt(float64(start + directN - 1))
		ratio := cmplx.Abs(tTail) / cmplx.Abs(tPrev)
		if ratio > 1.0001 {
			return ScalarReal(math.Inf(1)), nil
		}
		tailIntegral := integrateComplex(func(x float64) complex128 {
			z, err := termAt(x)
			if err != nil {
				return 0
			}
			return z
		}, float64(start+directN), math.Inf(1))

		if !math.IsNaN(real(tailIntegral)) && !math.IsInf(real(tailIntegral), 0) && imag(tTail) == 0 && real(tTail)*real(tPrev) > 0 {
			N := float64(start + directN)
			fN, _ := termAt(N)
			h := 1e-4 * N
			fP, _ := termAt(N + h)
			fM, _ := termAt(N - h)
			d1 := (fP - fM) / complex(2*h, 0)
			h3 := 1e-2 * N
			fp2, _ := termAt(N + 2*h3)
			fp1, _ := termAt(N + h3)
			fm1, _ := termAt(N - h3)
			fm2, _ := termAt(N - 2*h3)
			d3 := (fp2 - 2*fp1 + 2*fm1 - fm2) / complex(2*h3*h3*h3, 0)
			emSum := partial + tailIntegral + 0.5*fN - d1/12.0 + d3/720.0
			return ScalarComplex(emSum), nil
		}
	}

	const wynnTerms = 60
	partials := make([]complex128, wynnTerms)
	var s complex128
	for k := 0; k < wynnTerms; k++ {
		t, err := termAt(float64(start) + float64(k))
		if err != nil {
			return Value{}, err
		}
		s += t
		partials[k] = s
	}
	eps := make([][]complex128, wynnTerms+1)
	for i := range eps {
		eps[i] = make([]complex128, wynnTerms+1)
	}
	for i := 0; i < wynnTerms; i++ {
		eps[i][1] = partials[i]
	}
	best := partials[wynnTerms-1]
	for j := 2; j <= wynnTerms; j++ {
		for i := 0; i <= wynnTerms-j; i++ {
			diff := eps[i+1][j-1] - eps[i][j-1]
			if cmplx.Abs(diff) < 1e-30 {
				eps[i][j] = eps[i+1][j-1]
			} else {
				eps[i][j] = eps[i+1][j-2] + 1/diff
			}
			if j%2 == 1 && !math.IsNaN(real(eps[i][j])) && !math.IsInf(real(eps[i][j]), 0) {
				best = eps[i][j]
			}
		}
	}
	return ScalarComplex(best), nil
}

func buildPolyAST(coeffs []float64, varName string) ASTNode {
	var terms []ASTNode
	for deg := len(coeffs) - 1; deg >= 0; deg-- {
		c := coeffs[deg]
		if math.Abs(c) <= 1e-13 {
			continue
		}
		cClean := c
		if math.Abs(c-math.Round(c)) <= 1e-11 {
			cClean = math.Round(c)
		}
		if deg == 0 {
			terms = append(terms, makeSignedNum(cClean))
			continue
		}
		var xPart ASTNode = IdentNode{Name: varName}
		if deg > 1 {
			xPart = BinaryNode{Op: "^", Left: xPart, Right: NumNode{Val: float64(deg)}}
		}
		term := simplifyAST(BinaryNode{Op: "*", Left: makeSignedNum(cClean), Right: xPart})
		terms = append(terms, term)
	}
	if len(terms) == 0 {
		return NumNode{Val: 0}
	}
	acc := terms[0]
	for i := 1; i < len(terms); i++ {
		acc = simplifyAST(BinaryNode{Op: "+", Left: acc, Right: terms[i]})
	}
	return acc
}

func copyScope(scope map[string]Value) map[string]Value {
	out := make(map[string]Value, len(scope)+2)
	for k, v := range scope {
		out[k] = v
	}
	return out
}
