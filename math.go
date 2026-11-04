package main

import (
	"fmt"
	"go/ast"
	"go/parser"
	"go/token"
	"math"
	"strconv"
	"strings"
)

var constants = map[string]float64{
	"pi":      math.Pi,
	"PI":      math.Pi,
	"e":       math.E,
	"E":       math.E,
	"phi":     math.Phi,
	"tau":     2 * math.Pi,
	"inf":     math.Inf(1),
	"infinity": math.Inf(1),
}

var singleArgFuncs = map[string]func(float64) float64{
	"sin":     math.Sin,
	"cos":     math.Cos,
	"tan":     math.Tan,
	"asin":    math.Asin,
	"acos":    math.Acos,
	"atan":    math.Atan,
	"sinh":    math.Sinh,
	"cosh":    math.Cosh,
	"tanh":    math.Tanh,
	"asinh":   math.Asinh,
	"acosh":   math.Acosh,
	"atanh":   math.Atanh,
	"sqrt":    math.Sqrt,
	"cbrt":    math.Cbrt,
	"abs":     math.Abs,
	"floor":   math.Floor,
	"ceil":    math.Ceil,
	"round":   math.Round,
	"ln":      math.Log,
	"log":     math.Log10,
	"log10":   math.Log10,
	"log2":    math.Log2,
	"exp":     math.Exp,
	"deg2rad": func(x float64) float64 { return x * math.Pi / 180 },
	"rad2deg": func(x float64) float64 { return x * 180 / math.Pi },
	"fact":    factorial,
}

var twoArgFuncs = map[string]func(float64, float64) float64{
	"pow":   math.Pow,
	"mod":   math.Mod,
	"min":   math.Min,
	"max":   math.Max,
	"atan2": math.Atan2,
	"hypot": math.Hypot,
	"logb":  func(base, x float64) float64 { return math.Log(x) / math.Log(base) },
}

func factorial(x float64) float64 {
	if x < 0 || math.IsNaN(x) || math.IsInf(x, 0) {
		return math.NaN()
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

func EvaluateExpression(expr string) float64 {
	val, err := evalString(expr)
	if err != nil {
		return 0
	}
	return val
}

func evalString(expr string) (float64, error) {
	sanitized := strings.ReplaceAll(expr, "^", "**")
	sanitized = strings.ReplaceAll(sanitized, "**", "^")

	node, err := parser.ParseExpr(sanitized)
	if err != nil {
		return 0, err
	}
	return evalNode(node)
}

func evalNode(node ast.Expr) (float64, error) {
	switch n := node.(type) {
	case *ast.BasicLit:
		switch n.Kind {
		case token.INT, token.FLOAT:
			val, err := strconv.ParseFloat(n.Value, 64)
			if err != nil {
				return 0, err
			}
			return val, nil
		default:
			return 0, fmt.Errorf("unsupported literal %v", n.Value)
		}

	case *ast.Ident:
		if val, ok := constants[n.Name]; ok {
			return val, nil
		}
		if val, ok := constants[strings.ToLower(n.Name)]; ok {
			return val, nil
		}
		return 0, fmt.Errorf("unknown identifier %s", n.Name)

	case *ast.ParenExpr:
		return evalNode(n.X)

	case *ast.UnaryExpr:
		x, err := evalNode(n.X)
		if err != nil {
			return 0, err
		}
		switch n.Op {
		case token.SUB:
			return -x, nil
		case token.ADD:
			return x, nil
		default:
			return 0, fmt.Errorf("unsupported unary operator %s", n.Op)
		}

	case *ast.BinaryExpr:
		left, err := evalNode(n.X)
		if err != nil {
			return 0, err
		}
		right, err := evalNode(n.Y)
		if err != nil {
			return 0, err
		}

		switch n.Op {
		case token.ADD:
			return left + right, nil
		case token.SUB:
			return left - right, nil
		case token.MUL:
			return left * right, nil
		case token.QUO:
			return left / right, nil
		case token.REM:
			return math.Mod(left, right), nil
		case token.XOR:
			return math.Pow(left, right), nil
		default:
			return 0, fmt.Errorf("unsupported binary operator %s", n.Op)
		}

	case *ast.CallExpr:
		ident, ok := n.Fun.(*ast.Ident)
		if !ok {
			return 0, fmt.Errorf("unsupported function call format")
		}

		fnName := strings.ToLower(ident.Name)
		if fn, exists := singleArgFuncs[fnName]; exists {
			if len(n.Args) != 1 {
				return 0, fmt.Errorf("%s takes 1 argument", fnName)
			}
			arg, err := evalNode(n.Args[0])
			if err != nil {
				return 0, err
			}
			return fn(arg), nil
		}

		if fn, exists := twoArgFuncs[fnName]; exists {
			if len(n.Args) != 2 {
				return 0, fmt.Errorf("%s takes 2 arguments", fnName)
			}
			arg1, err := evalNode(n.Args[0])
			if err != nil {
				return 0, err
			}
			arg2, err := evalNode(n.Args[1])
			if err != nil {
				return 0, err
			}
			return fn(arg1, arg2), nil
		}

		return 0, fmt.Errorf("unknown function %s", ident.Name)

	default:
		return 0, fmt.Errorf("unsupported expression")
	}
}
