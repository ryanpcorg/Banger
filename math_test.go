package main

import (
	"math"
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
	}

	for _, tt := range tests {
		val := EvaluateExpression(tt.expr)
		if math.Abs(val-tt.expected) > 1e-9 {
			t.Errorf("expr %q: got %v, want %v", tt.expr, val, tt.expected)
		}
	}
}
