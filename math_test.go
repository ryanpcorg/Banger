package main

import (
	"math"
	"math/cmplx"
	"strings"
	"testing"
)

func TestPrecedenceAndImplicitMul(t *testing.T) {
	eng := NewEngine()
	tests := []struct {
		expr     string
		expected float64
	}{
		{"2 + 2", 4},
		{"2 * 3 ^ 2", 18},
		{"-2^2", -4},
		{"(-2)^2", 4},
		{"2^3^2", 512},
		{"2^-3", 0.125},
		{"5!", 120},
		{"7!!", 105},
		{"2pi", 2 * math.Pi},
		{"3(4 + 5)", 27},
		{"(2 + 3)(4 - 1)", 15},
		{"2sin(pi/6)", 1},
		{"0xFF + 0b1010 + 0o10", 255 + 10 + 8},
	}

	for _, tt := range tests {
		v, err := eng.Eval(tt.expr)
		if err != nil {
			t.Fatalf("expr %q failed: %v", tt.expr, err)
		}
		got, err := v.Real()
		if err != nil {
			t.Fatalf("expr %q not real: %v", tt.expr, err)
		}
		if math.Abs(got-tt.expected) > 1e-9 {
			t.Errorf("expr %q: got %v, want %v", tt.expr, got, tt.expected)
		}
	}
}

func TestComplexAnalysis(t *testing.T) {
	eng := NewEngine()

	v1, err := eng.Eval("e^(i*pi) + 1")
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(v1.S) > 1e-13 {
		t.Errorf("Euler identity got %v, want 0", v1.S)
	}

	v2, err := eng.Eval("sqrt(-4)")
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(v2.S-complex(0, 2)) > 1e-13 {
		t.Errorf("sqrt(-4) got %v, want 2i", v2.S)
	}

	v3, err := eng.Eval("i^i")
	if err != nil {
		t.Fatal(err)
	}
	wantIi := math.Exp(-math.Pi / 2)
	if cmplx.Abs(v3.S-complex(wantIi, 0)) > 1e-12 {
		t.Errorf("i^i got %v, want %v", v3.S, wantIi)
	}
}

func TestCalculusAndDifferentialEquations(t *testing.T) {
	eng := NewEngine()

	d1, err := eng.Eval("diff(x^3 + 2x^2 - 5x + 7, x)")
	if err != nil {
		t.Fatal(err)
	}
	if d1.Sym != "3*x^2 + 4*x - 5" {
		t.Errorf("diff polynomial got %q, want %q", d1.Sym, "3*x^2 + 4*x - 5")
	}

	_, err = eng.Eval("f(x) = x^4 - 2x^2")
	if err != nil {
		t.Fatal(err)
	}
	fp, err := eng.Eval("f'(2)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := fp.Real(); math.Abs(r-24) > 1e-9 {
		t.Errorf("f'(2) got %v, want 24", r)
	}

	tay, err := eng.Eval("taylor(sin(x), x, 0, 5)")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(tay.Sym, "x^3") || !strings.Contains(tay.Sym, "x^5") {
		t.Errorf("unexpected taylor(sin(x)): %q", tay.Sym)
	}

	gauss, err := eng.Eval("int(e^(-x^2), x, -inf, inf)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := gauss.Real(); math.Abs(r-math.SqrtPi) > 1e-9 {
		t.Errorf("Gaussian integral got %v, want %v", r, math.SqrtPi)
	}

	lim1, err := eng.Eval("lim(sin(x)/x, x, 0)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := lim1.Real(); math.Abs(r-1) > 1e-7 {
		t.Errorf("lim sin(x)/x got %v, want 1", r)
	}

	limDNE, err := eng.Eval("lim(1/x, x, 0)")
	if err != nil {
		t.Fatal(err)
	}
	if limDNE.Sym != "DNE" {
		t.Errorf("lim(1/x, x, 0) got %v, want DNE", limDNE.String())
	}

	basel, err := eng.Eval("sum(1/n^2, n, 1, inf)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := basel.Real(); math.Abs(r-(math.Pi*math.Pi/6)) > 1e-7 {
		t.Errorf("Basel series got %v, want %v", r, math.Pi*math.Pi/6)
	}

	ode, err := eng.Eval("rk4(y, x, y, 0, 1, 1, 1000)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := ode.Real(); math.Abs(r-math.E) > 1e-9 {
		t.Errorf("rk4 y'=y got %v, want %v", r, math.E)
	}

	curlVal, err := eng.Eval("curl([-y, x, 0], [x, y, z], [1, 2, 3])")
	if err != nil {
		t.Fatal(err)
	}
	if len(curlVal.V) != 3 || math.Abs(real(curlVal.V[2])-2) > 1e-9 {
		t.Errorf("curl got %v, want [0, 0, 2]", curlVal)
	}
}

func TestEquationSolver(t *testing.T) {
	eng := NewEngine()

	cubic, err := eng.Eval("solve(x^3 - 6x^2 + 11x - 6 = 0, x)")
	if err != nil {
		t.Fatal(err)
	}
	if cubic.Kind != KindVector || len(cubic.V) != 3 {
		t.Fatalf("expected 3 cubic roots, got %v", cubic)
	}
	for idx, want := range []float64{1, 2, 3} {
		if cmplx.Abs(cubic.V[idx]-complex(want, 0)) > 1e-9 {
			t.Errorf("root %d: got %v, want %v", idx, cubic.V[idx], want)
		}
	}

	unitImag, err := eng.Eval("solve(x^2 + 1 = 0, x)")
	if err != nil {
		t.Fatal(err)
	}
	if len(unitImag.V) != 2 || cmplx.Abs(unitImag.V[0]-complex(0, -1)) > 1e-9 || cmplx.Abs(unitImag.V[1]-complex(0, 1)) > 1e-9 {
		t.Errorf("solve(x^2+1=0) got %v", unitImag)
	}

	dottie, err := eng.Eval("solve(cos(x) = x, x, 0, 2)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := dottie.Real(); math.Abs(r-0.739085133215) > 1e-8 {
		t.Errorf("Dottie number got %v", r)
	}
}

func TestLinearAlgebra(t *testing.T) {
	eng := NewEngine()

	_, err := eng.Eval("A = [1, 2; 3, 4]")
	if err != nil {
		t.Fatal(err)
	}
	detVal, err := eng.Eval("det(A)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := detVal.Real(); math.Abs(r - -2) > 1e-9 {
		t.Errorf("det(A) got %v, want -2", r)
	}

	prodIdent, err := eng.Eval("A * A^-1")
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(prodIdent.M[0][0]-1) > 1e-12 || cmplx.Abs(prodIdent.M[0][1]) > 1e-12 || cmplx.Abs(prodIdent.M[1][1]-1) > 1e-12 {
		t.Errorf("A * A^-1 not identity: %v", prodIdent)
	}

	rotEig, err := eng.Eval("eig([0, -1; 1, 0])")
	if err != nil {
		t.Fatal(err)
	}
	if len(rotEig.V) != 2 || cmplx.Abs(rotEig.V[0]-complex(0, -1)) > 1e-9 || cmplx.Abs(rotEig.V[1]-complex(0, 1)) > 1e-9 {
		t.Errorf("rotation eigenvalues got %v, want [-i, i]", rotEig)
	}

	permEig, err := eng.Eval("eig([0, 0, 1; 1, 0, 0; 0, 1, 0])")
	if err != nil {
		t.Fatal(err)
	}
	if len(permEig.V) != 3 {
		t.Fatalf("expected 3 eigenvalues for 3x3 permutation matrix")
	}

	sol, err := eng.Eval("linsolve([2, 1; 5, 3], [4, 11])")
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(sol.V[0]-1) > 1e-9 || cmplx.Abs(sol.V[1]-2) > 1e-9 {
		t.Errorf("linsolve got %v, want [1, 2]", sol)
	}

	adjCheck, err := eng.Eval("A * adj(A)")
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(adjCheck.M[0][0] - -2) > 1e-9 || cmplx.Abs(adjCheck.M[1][1] - -2) > 1e-9 || cmplx.Abs(adjCheck.M[0][1]) > 1e-9 {
		t.Errorf("A * adj(A) != det(A)*I: %v", adjCheck)
	}
}

func TestNumberTheorySpecialAndStats(t *testing.T) {
	eng := NewEngine()

	crtVal, err := eng.Eval("crt([2, 3, 2], [3, 5, 7])")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := crtVal.Real(); r != 23 {
		t.Errorf("crt got %v, want 23", r)
	}

	facVal, err := eng.Eval("factor(360)")
	if err != nil || facVal.Sym != "2^3 * 3^2 * 5" {
		t.Errorf("factor(360) got %v, err=%v", facVal, err)
	}

	z2, err := eng.Eval("zeta(2)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := z2.Real(); math.Abs(r-(math.Pi*math.Pi/6)) > 1e-11 {
		t.Errorf("zeta(2) got %v, want %v", r, math.Pi*math.Pi/6)
	}

	zm1, err := eng.Eval("zeta(-1)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := zm1.Real(); math.Abs(r - -1.0/12.0) > 1e-11 {
		t.Errorf("zeta(-1) got %v, want -1/12", r)
	}

	lw, err := eng.Eval("lambertw(e)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := lw.Real(); math.Abs(r-1) > 1e-12 {
		t.Errorf("lambertw(e) got %v, want 1", r)
	}

	p100, err := eng.Eval("partition(100)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := p100.Real(); r != 190569292 {
		t.Errorf("partition(100) got %v, want 190569292", r)
	}

	z975, err := eng.Eval("norminv(0.975)")
	if err != nil {
		t.Fatal(err)
	}
	if r, _ := z975.Real(); math.Abs(r-1.95996398454) > 1e-6 {
		t.Errorf("norminv(0.975) got %v, want ~1.959964", r)
	}

	reg, err := eng.Eval("linreg([1, 2, 3, 4], [3, 5, 7, 9])")
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(real(reg.V[0])-2) > 1e-9 || math.Abs(real(reg.V[1])-1) > 1e-9 || math.Abs(real(reg.V[2])-1) > 1e-9 {
		t.Errorf("linreg got %v, want [2, 1, 1]", reg)
	}
}

func TestBrailleGraphing(t *testing.T) {
	output := RenderGraph("x^2 - 4", -5, 5, 48, 14)
	if !strings.Contains(output, "f(x) = x^2 - 4") {
		t.Fatalf("expected graph title in output, got %s", output)
	}
	if !strings.Contains(output, "Zeros in window: x = -2, 2") {
		t.Fatalf("expected detected zeros in output, got %s", output)
	}
	hasBraille := false
	for _, r := range output {
		if r >= 0x2801 && r <= 0x28FF {
			hasBraille = true
			break
		}
	}
	if !hasBraille {
		t.Fatalf("expected Unicode Braille curve characters in graph output")
	}
}
