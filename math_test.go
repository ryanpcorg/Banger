package main

import (
	"math"
	"strings"
	"testing"
)

func TestEvaluateExpression(t *testing.T) {
	tests := []struct {
		expr     string
		expected float64
	}{
		{"2 + 2", 4},
		{"2 * 3 + 4", 10},
		{"2 ^ 3", 8},
		{"sqrt(16)", 4},
		{"sin(pi / 2)", 1},
		{"cos(0)", 1},
		{"fact(5)", 120},
		{"pow(2, 4)", 16},
		{"log10(100)", 2},
		{"ln(e)", 1},
		{"abs(-42)", 42},
		{"max(10, 20)", 20},
		{"min(10, 20)", 10},
		{"10 % 3", 1},
		{"gcd(48, 18)", 6},
		{"lcm(12, 18)", 36},
		{"npr(5, 2)", 20},
		{"ncr(5, 2)", 10},
		{"root(3, 27)", 3},
		{"gamma(5)", 24},
	}

	for _, tt := range tests {
		val := EvaluateExpression(tt.expr)
		if math.Abs(val-tt.expected) > 1e-9 {
			t.Errorf("expr %q: got %v, want %v", tt.expr, val, tt.expected)
		}
	}
}

func TestRenderGraph(t *testing.T) {
	output := RenderGraph("x^2", -5, 5, 40, 15)
	if !strings.Contains(output, "f(x) = x^2") {
		t.Fatalf("expected graph title in output, got %s", output)
	}
	if !strings.Contains(output, "*") {
		t.Fatalf("expected plotted stars in graph output, got %s", output)
	}
}
