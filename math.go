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
	"pi":       math.Pi,
	"PI":       math.Pi,
	"e":        math.E,
	"E":        math.E,
	"phi":      math.Phi,
	"tau":      2 * math.Pi,
	"inf":      math.Inf(1),
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
	"gamma":   math.Gamma,
	"erf":     math.Erf,
	"erfc":    math.Erfc,
}

var twoArgFuncs = map[string]func(float64, float64) float64{
	"pow":   math.Pow,
	"mod":   math.Mod,
	"min":   math.Min,
	"max":   math.Max,
	"atan2": math.Atan2,
	"hypot": math.Hypot,
	"logb":  func(base, x float64) float64 { return math.Log(x) / math.Log(base) },
	"root":  func(n, x float64) float64 { return math.Pow(x, 1.0/n) },
	"npr":   nPr,
	"ncr":   nCr,
	"gcd":   gcd,
	"lcm":   lcm,
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

func gcd(a, b float64) float64 {
	x := int64(math.Abs(math.Round(a)))
	y := int64(math.Abs(math.Round(b)))
	for y != 0 {
		x, y = y, x%y
	}
	return float64(x)
}

func lcm(a, b float64) float64 {
	x := int64(math.Abs(math.Round(a)))
	y := int64(math.Abs(math.Round(b)))
	if x == 0 || y == 0 {
		return 0
	}
	g := gcd(float64(x), float64(y))
	return float64((x / int64(g)) * y)
}

func nPr(n, r float64) float64 {
	if n < 0 || r < 0 || r > n {
		return 0
	}
	return math.Round(factorial(n) / factorial(n-r))
}

func nCr(n, r float64) float64 {
	if n < 0 || r < 0 || r > n {
		return 0
	}
	return math.Round(factorial(n) / (factorial(r) * factorial(n-r)))
}

func EvaluateExpression(expr string) float64 {
	val, err := EvaluateWithScope(expr, nil)
	if err != nil {
		return 0
	}
	return val
}

func EvaluateWithScope(expr string, scope map[string]float64) (float64, error) {
	sanitized := strings.ReplaceAll(expr, "^", "**")
	sanitized = strings.ReplaceAll(sanitized, "**", "^")

	node, err := parser.ParseExpr(sanitized)
	if err != nil {
		return 0, err
	}
	return evalNode(node, scope)
}

func evalNode(node ast.Expr, scope map[string]float64) (float64, error) {
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
		if scope != nil {
			if val, ok := scope[n.Name]; ok {
				return val, nil
			}
			if val, ok := scope[strings.ToLower(n.Name)]; ok {
				return val, nil
			}
		}
		if val, ok := constants[n.Name]; ok {
			return val, nil
		}
		if val, ok := constants[strings.ToLower(n.Name)]; ok {
			return val, nil
		}
		return 0, fmt.Errorf("unknown identifier %s", n.Name)

	case *ast.ParenExpr:
		return evalNode(n.X, scope)

	case *ast.UnaryExpr:
		x, err := evalNode(n.X, scope)
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
		left, err := evalNode(n.X, scope)
		if err != nil {
			return 0, err
		}
		right, err := evalNode(n.Y, scope)
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
			arg, err := evalNode(n.Args[0], scope)
			if err != nil {
				return 0, err
			}
			return fn(arg), nil
		}

		if fn, exists := twoArgFuncs[fnName]; exists {
			if len(n.Args) != 2 {
				return 0, fmt.Errorf("%s takes 2 arguments", fnName)
			}
			arg1, err := evalNode(n.Args[0], scope)
			if err != nil {
				return 0, err
			}
			arg2, err := evalNode(n.Args[1], scope)
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

func RenderGraph(expr string, xMin, xMax float64, width, height int) string {
	if width <= 0 {
		width = 60
	}
	if height <= 0 {
		height = 20
	}
	if xMax <= xMin {
		xMin = -10
		xMax = 10
	}

	xVals := make([]float64, width)
	yVals := make([]float64, width)
	valid := make([]bool, width)

	stepX := (xMax - xMin) / float64(width-1)
	yMin := math.Inf(1)
	yMax := math.Inf(-1)

	scope := make(map[string]float64)

	for i := 0; i < width; i++ {
		x := xMin + float64(i)*stepX
		xVals[i] = x
		scope["x"] = x
		scope["X"] = x

		y, err := EvaluateWithScope(expr, scope)
		if err == nil && !math.IsNaN(y) && !math.IsInf(y, 0) {
			yVals[i] = y
			valid[i] = true
			if y < yMin {
				yMin = y
			}
			if y > yMax {
				yMax = y
			}
		}
	}

	if math.IsInf(yMin, 0) || math.IsInf(yMax, 0) {
		return "Unable to plot: no valid real points in domain.\n"
	}

	if yMin == yMax {
		yMin -= 1
		yMax += 1
	}

	grid := make([][]rune, height)
	for r := 0; r < height; r++ {
		grid[r] = make([]rune, width)
		for c := 0; c < width; c++ {
			grid[r][c] = ' '
		}
	}

	yZeroRow := -1
	if yMin <= 0 && yMax >= 0 {
		ratio := (0 - yMin) / (yMax - yMin)
		yZeroRow = height - 1 - int(math.Round(ratio*float64(height-1)))
		for c := 0; c < width; c++ {
			grid[yZeroRow][c] = '-'
		}
	}

	xZeroCol := -1
	if xMin <= 0 && xMax >= 0 {
		ratio := (0 - xMin) / (xMax - xMin)
		xZeroCol = int(math.Round(ratio * float64(width-1)))
		for r := 0; r < height; r++ {
			if grid[r][xZeroCol] == '-' {
				grid[r][xZeroCol] = '+'
			} else {
				grid[r][xZeroCol] = '|'
			}
		}
	}

	for c := 0; c < width; c++ {
		if !valid[c] {
			continue
		}
		ratio := (yVals[c] - yMin) / (yMax - yMin)
		r := height - 1 - int(math.Round(ratio*float64(height-1)))
		if r >= 0 && r < height {
			grid[r][c] = '*'
		}
	}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("f(x) = %s   [x: %.2f to %.2f, y: %.2f to %.2f]\n", expr, xMin, xMax, yMin, yMax))
	for r := 0; r < height; r++ {
		sb.WriteString(string(grid[r]))
		sb.WriteRune('\n')
	}

	return sb.String()
}
