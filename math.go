package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"strconv"
	"strings"
	"unicode"
)

type ValueKind int

const (
	KindScalar ValueKind = iota
	KindVector
	KindMatrix
	KindSymbolic
)

type Value struct {
	Kind ValueKind
	S    complex128
	V    []complex128
	M    [][]complex128
	Sym  string
	Node ASTNode
}

type UserFunc struct {
	Params []string
	Body   ASTNode
	Raw    string
}

type HistoryEntry struct {
	Input string
	Val   Value
}

type AngleMode int

const (
	ModeRad AngleMode = iota
	ModeDeg
)

type Engine struct {
	Vars      map[string]Value
	Funcs     map[string]UserFunc
	History   []HistoryEntry
	Mode      AngleMode
	LastValue Value
}

var defaultEngine = NewEngine()

func NewEngine() *Engine {
	return &Engine{
		Vars:      make(map[string]Value),
		Funcs:     make(map[string]UserFunc),
		History:   make([]HistoryEntry, 0),
		Mode:      ModeRad,
		LastValue: ScalarReal(0),
	}
}

func (e *Engine) Reset() {
	e.Vars = make(map[string]Value)
	e.Funcs = make(map[string]UserFunc)
	e.History = make([]HistoryEntry, 0)
	e.LastValue = ScalarReal(0)
}

func ScalarReal(x float64) Value {
	if x == 0 {
		x = 0
	}
	return Value{Kind: KindScalar, S: complex(x, 0)}
}

func ScalarComplex(z complex128) Value {
	return Value{Kind: KindScalar, S: cleanComplex(z)}
}

func VectorVal(v []complex128) Value {
	out := make([]complex128, len(v))
	for i, z := range v {
		out[i] = cleanComplex(z)
	}
	return Value{Kind: KindVector, V: out}
}

func MatrixVal(m [][]complex128) Value {
	out := make([][]complex128, len(m))
	for i := range m {
		out[i] = make([]complex128, len(m[i]))
		for j, z := range m[i] {
			out[i][j] = cleanComplex(z)
		}
	}
	return Value{Kind: KindMatrix, M: out}
}

func SymbolicVal(sym string, node ASTNode) Value {
	return Value{Kind: KindSymbolic, Sym: sym, Node: node}
}

func cleanComplex(z complex128) complex128 {
	r, im := real(z), imag(z)
	if math.IsNaN(r) || math.IsNaN(im) {
		return z
	}
	if r != 0 && !math.IsInf(r, 0) && math.Abs(im) <= 1e-14*math.Abs(r) {
		im = 0
	}
	if im != 0 && !math.IsInf(im, 0) && math.Abs(r) <= 1e-14*math.Abs(im) {
		r = 0
	}
	if r == 0 {
		r = 0
	}
	if im == 0 {
		im = 0
	}
	return complex(r, im)
}

func (v Value) IsReal() bool {
	return v.Kind == KindScalar && math.Abs(imag(v.S)) <= 1e-13*math.Max(1.0, math.Abs(real(v.S)))
}

func (v Value) Real() (float64, error) {
	if v.Kind != KindScalar {
		return 0, fmt.Errorf("expected real scalar")
	}
	if !v.IsReal() {
		return 0, fmt.Errorf("expected real value, got complex %s", formatComplex(v.S, false))
	}
	return real(v.S), nil
}

func formatFloatClean(x float64) string {
	if math.IsNaN(x) {
		return "NaN"
	}
	if math.IsInf(x, 1) {
		return "inf"
	}
	if math.IsInf(x, -1) {
		return "-inf"
	}
	if x == 0 {
		return "0"
	}
	rounded := math.Round(x)
	if rounded != 0 && math.Abs(x) < 1e15 && math.Abs(x-rounded) <= 1e-12*math.Abs(x) {
		return strconv.FormatFloat(rounded, 'f', -1, 64)
	}
	return strconv.FormatFloat(x, 'g', 12, 64)
}

func exactHint(x float64) string {
	if math.IsNaN(x) || math.IsInf(x, 0) || math.Abs(x) < 1e-9 || math.Abs(x) > 1e7 {
		return ""
	}
	if math.Abs(x-math.Round(x)) <= 1e-11*math.Max(1.0, math.Abs(x)) {
		return ""
	}

	if num, den, ok := bestRational(x, 1000, 1e-11); ok && den > 1 {
		decStr := formatFloatClean(x)
		if len(decStr) > 5 || (den != 2 && den != 4 && den != 5 && den != 8 && den != 10) {
			return fmt.Sprintf("%d/%d", num, den)
		}
	}

	if num, den, ok := bestRational(x/math.Pi, 72, 1e-11); ok && num != 0 {
		return formatSymbolicMultiple(num, den, "pi")
	}

	if num, den, ok := bestRational(x/(math.Pi*math.Pi), 96, 1e-11); ok && num != 0 && absInt64(num) <= 24 {
		return formatSymbolicMultiple(num, den, "pi^2")
	}

	if math.Abs(math.Abs(x)-math.Sqrt(math.Pi)) <= 1e-11 {
		if x < 0 {
			return "-sqrt(pi)"
		}
		return "sqrt(pi)"
	}
	if num, den, ok := bestRational(x/math.Sqrt(math.Pi), 24, 1e-11); ok && num != 0 && absInt64(num) <= 12 {
		return formatSymbolicMultiple(num, den, "sqrt(pi)")
	}

	radicals := []int{2, 3, 5, 6, 7, 10}
	for _, d := range radicals {
		rt := math.Sqrt(float64(d))
		if num, den, ok := bestRational(x/rt, 24, 1e-11); ok && num != 0 && absInt64(num) <= 24 {
			return formatSymbolicMultiple(num, den, fmt.Sprintf("sqrt(%d)", d))
		}
	}

	return ""
}

func absInt64(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

func formatSymbolicMultiple(num, den int64, sym string) string {
	sign := ""
	if num < 0 {
		sign = "-"
		num = -num
	}
	var top string
	if num == 1 {
		top = sym
	} else {
		top = fmt.Sprintf("%d*%s", num, sym)
	}
	if den == 1 {
		return sign + top
	}
	return fmt.Sprintf("%s%s/%d", sign, top, den)
}

func bestRational(x float64, maxDen int64, tol float64) (int64, int64, bool) {
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return 0, 1, false
	}
	sign := int64(1)
	if x < 0 {
		sign = -1
		x = -x
	}
	h0, h1 := int64(0), int64(1)
	k0, k1 := int64(1), int64(0)
	b := x
	for i := 0; i < 40; i++ {
		a := math.Floor(b)
		if a > 1e12 {
			break
		}
		ai := int64(a)
		h2 := ai*h1 + h0
		k2 := ai*k1 + k0
		if k2 > maxDen || k2 <= 0 {
			break
		}
		h0, h1 = h1, h2
		k0, k1 = k1, k2
		approx := float64(h1) / float64(k1)
		if math.Abs(x-approx) <= tol*math.Max(1.0, x) {
			return sign * h1, k1, true
		}
		rem := b - a
		if rem < 1e-14 {
			break
		}
		b = 1.0 / rem
	}
	if k1 > 0 && k1 <= maxDen {
		approx := float64(h1) / float64(k1)
		if math.Abs(x-approx) <= tol*math.Max(1.0, x) {
			return sign * h1, k1, true
		}
	}
	return 0, 1, false
}

func formatComplex(z complex128, withHint bool) string {
	z = cleanComplex(z)
	r, im := real(z), imag(z)
	if im == 0 {
		s := formatFloatClean(r)
		if withHint {
			if h := exactHint(r); h != "" && h != s {
				return fmt.Sprintf("%s  [= %s]", s, h)
			}
		}
		return s
	}
	if r == 0 {
		if im == 1 {
			return "i"
		}
		if im == -1 {
			return "-i"
		}
		return formatFloatClean(im) + "i"
	}
	rStr := formatFloatClean(r)
	absIm := math.Abs(im)
	var imStr string
	if absIm == 1 {
		imStr = "i"
	} else {
		imStr = formatFloatClean(absIm) + "i"
	}
	if im < 0 {
		return rStr + " - " + imStr
	}
	return rStr + " + " + imStr
}

func (v Value) String() string {
	switch v.Kind {
	case KindScalar:
		return formatComplex(v.S, true)
	case KindSymbolic:
		return v.Sym
	case KindVector:
		parts := make([]string, len(v.V))
		for i, z := range v.V {
			parts[i] = formatComplex(z, false)
		}
		return "[" + strings.Join(parts, ", ") + "]"
	case KindMatrix:
		return formatMatrix(v.M)
	default:
		return "0"
	}
}

func formatMatrix(m [][]complex128) string {
	rows := len(m)
	if rows == 0 {
		return "[]"
	}
	cols := len(m[0])
	cells := make([][]string, rows)
	widths := make([]int, cols)
	for r := 0; r < rows; r++ {
		cells[r] = make([]string, cols)
		for c := 0; c < cols; c++ {
			s := formatComplex(m[r][c], false)
			cells[r][c] = s
			if len(s) > widths[c] {
				widths[c] = len(s)
			}
		}
	}
	var sb strings.Builder
	for r := 0; r < rows; r++ {
		left, right := "│", "│"
		if rows == 1 {
			left, right = "[", "]"
		} else if r == 0 {
			left, right = "┌", "┐"
		} else if r == rows-1 {
			left, right = "└", "┘"
		}
		sb.WriteString(left + " ")
		for c := 0; c < cols; c++ {
			pad := widths[c] - len(cells[r][c])
			sb.WriteString(strings.Repeat(" ", pad) + cells[r][c])
			if c < cols-1 {
				sb.WriteString("   ")
			}
		}
		sb.WriteString(" " + right)
		if r < rows-1 {
			sb.WriteRune('\n')
		}
	}
	return sb.String()
}

var builtinConstants = map[string]complex128{
	"pi":       complex(math.Pi, 0),
	"tau":      complex(2*math.Pi, 0),
	"e":        complex(math.E, 0),
	"phi":      complex(math.Phi, 0),
	"i":        complex(0, 1),
	"j":        complex(0, 1),
	"euler":    complex(0.57721566490153286, 0),
	"catalan":  complex(0.91596559417721901, 0),
	"inf":      complex(math.Inf(1), 0),
	"infinity": complex(math.Inf(1), 0),
	"c_light":  complex(299792458, 0),
	"g_grav":   complex(6.67430e-11, 0),
	"planck":   complex(6.62607015e-34, 0),
	"hbar":     complex(1.054571817e-34, 0),
	"boltz":    complex(1.380649e-23, 0),
	"avogadro": complex(6.02214076e23, 0),
}

type TokenKind int

const (
	TokEOF TokenKind = iota
	TokNumber
	TokImag
	TokIdent
	TokHistoryRef
	TokPlus
	TokMinus
	TokStar
	TokSlash
	TokPercent
	TokCaret
	TokBang
	TokDoubleBang
	TokPrime
	TokLParen
	TokRParen
	TokLBracket
	TokRBracket
	TokComma
	TokSemi
	TokAssign
	TokEq
	TokNeq
	TokLt
	TokLte
	TokGt
	TokGte
)

type Token struct {
	Kind TokenKind
	Text string
	Num  float64
}

func tokenize(input string) ([]Token, error) {
	var tokens []Token
	n := len(input)
	idx := 0

	for idx < n {
		ch := input[idx]
		if ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n' {
			idx++
			continue
		}

		if ch == '$' && idx+1 < n && isDigit(input[idx+1]) {
			start := idx + 1
			idx++
			for idx < n && isDigit(input[idx]) {
				idx++
			}
			tokens = append(tokens, Token{Kind: TokHistoryRef, Text: input[start:idx]})
			continue
		}

		if ch == '0' && idx+2 <= n {
			if idx+2 < n && (input[idx+1] == 'x' || input[idx+1] == 'X') && isHexDigit(input[idx+2]) {
				start := idx + 2
				idx += 2
				for idx < n && (isHexDigit(input[idx]) || input[idx] == '_') {
					idx++
				}
				clean := strings.ReplaceAll(input[start:idx], "_", "")
				val, err := strconv.ParseInt(clean, 16, 64)
				if err != nil {
					return nil, err
				}
				tokens = append(tokens, Token{Kind: TokNumber, Text: input[start-2 : idx], Num: float64(val)})
				continue
			}
			if idx+2 < n && (input[idx+1] == 'b' || input[idx+1] == 'B') && (input[idx+2] == '0' || input[idx+2] == '1') {
				start := idx + 2
				idx += 2
				for idx < n && (input[idx] == '0' || input[idx] == '1' || input[idx] == '_') {
					idx++
				}
				clean := strings.ReplaceAll(input[start:idx], "_", "")
				val, err := strconv.ParseInt(clean, 2, 64)
				if err != nil {
					return nil, err
				}
				tokens = append(tokens, Token{Kind: TokNumber, Text: input[start-2 : idx], Num: float64(val)})
				continue
			}
			if idx+2 < n && (input[idx+1] == 'o' || input[idx+1] == 'O') && input[idx+2] >= '0' && input[idx+2] <= '7' {
				start := idx + 2
				idx += 2
				for idx < n && ((input[idx] >= '0' && input[idx] <= '7') || input[idx] == '_') {
					idx++
				}
				clean := strings.ReplaceAll(input[start:idx], "_", "")
				val, err := strconv.ParseInt(clean, 8, 64)
				if err != nil {
					return nil, err
				}
				tokens = append(tokens, Token{Kind: TokNumber, Text: input[start-2 : idx], Num: float64(val)})
				continue
			}
		}

		if isDigit(ch) || (ch == '.' && idx+1 < n && isDigit(input[idx+1])) {
			start := idx
			hasDot := false
			for idx < n && (isDigit(input[idx]) || input[idx] == '_' || (!hasDot && input[idx] == '.')) {
				if input[idx] == '.' {
					hasDot = true
				}
				idx++
			}
			if idx < n && (input[idx] == 'e' || input[idx] == 'E') {
				expIdx := idx + 1
				if expIdx < n && (input[expIdx] == '+' || input[expIdx] == '-') {
					expIdx++
				}
				if expIdx < n && isDigit(input[expIdx]) {
					idx = expIdx
					for idx < n && isDigit(input[idx]) {
						idx++
					}
				}
			}
			raw := strings.ReplaceAll(input[start:idx], "_", "")
			val, err := strconv.ParseFloat(raw, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid number %s", raw)
			}
			if idx < n && input[idx] == 'i' && (idx+1 == n || !isIdentPart(rune(input[idx+1]))) {
				idx++
				tokens = append(tokens, Token{Kind: TokImag, Text: raw + "i", Num: val})
			} else {
				tokens = append(tokens, Token{Kind: TokNumber, Text: raw, Num: val})
			}
			continue
		}

		if isIdentStart(rune(ch)) {
			start := idx
			idx++
			for idx < n && isIdentPart(rune(input[idx])) {
				idx++
			}
			tokens = append(tokens, Token{Kind: TokIdent, Text: input[start:idx]})
			continue
		}

		switch ch {
		case '+':
			tokens = append(tokens, Token{Kind: TokPlus, Text: "+"})
			idx++
		case '-':
			tokens = append(tokens, Token{Kind: TokMinus, Text: "-"})
			idx++
		case '*':
			if idx+1 < n && input[idx+1] == '*' {
				tokens = append(tokens, Token{Kind: TokCaret, Text: "^"})
				idx += 2
			} else {
				tokens = append(tokens, Token{Kind: TokStar, Text: "*"})
				idx++
			}
		case '/':
			tokens = append(tokens, Token{Kind: TokSlash, Text: "/"})
			idx++
		case '%':
			tokens = append(tokens, Token{Kind: TokPercent, Text: "%"})
			idx++
		case '^':
			tokens = append(tokens, Token{Kind: TokCaret, Text: "^"})
			idx++
		case '!':
			if idx+1 < n && input[idx+1] == '=' {
				tokens = append(tokens, Token{Kind: TokNeq, Text: "!="})
				idx += 2
			} else if idx+1 < n && input[idx+1] == '!' {
				tokens = append(tokens, Token{Kind: TokDoubleBang, Text: "!!"})
				idx += 2
			} else {
				tokens = append(tokens, Token{Kind: TokBang, Text: "!"})
				idx++
			}
		case '\'':
			tokens = append(tokens, Token{Kind: TokPrime, Text: "'"})
			idx++
		case '(':
			tokens = append(tokens, Token{Kind: TokLParen, Text: "("})
			idx++
		case ')':
			tokens = append(tokens, Token{Kind: TokRParen, Text: ")"})
			idx++
		case '[':
			tokens = append(tokens, Token{Kind: TokLBracket, Text: "["})
			idx++
		case ']':
			tokens = append(tokens, Token{Kind: TokRBracket, Text: "]"})
			idx++
		case ',':
			tokens = append(tokens, Token{Kind: TokComma, Text: ","})
			idx++
		case ';':
			tokens = append(tokens, Token{Kind: TokSemi, Text: ";"})
			idx++
		case '=':
			if idx+1 < n && input[idx+1] == '=' {
				tokens = append(tokens, Token{Kind: TokEq, Text: "=="})
				idx += 2
			} else {
				tokens = append(tokens, Token{Kind: TokAssign, Text: "="})
				idx++
			}
		case '<':
			if idx+1 < n && input[idx+1] == '=' {
				tokens = append(tokens, Token{Kind: TokLte, Text: "<="})
				idx += 2
			} else {
				tokens = append(tokens, Token{Kind: TokLt, Text: "<"})
				idx++
			}
		case '>':
			if idx+1 < n && input[idx+1] == '=' {
				tokens = append(tokens, Token{Kind: TokGte, Text: ">="})
				idx += 2
			} else {
				tokens = append(tokens, Token{Kind: TokGt, Text: ">"})
				idx++
			}
		default:
			return nil, fmt.Errorf("unexpected character %q", ch)
		}
	}

	tokens = append(tokens, Token{Kind: TokEOF})
	return tokens, nil
}

func isDigit(b byte) bool {
	return b >= '0' && b <= '9'
}

func isHexDigit(b byte) bool {
	return (b >= '0' && b <= '9') || (b >= 'a' && b <= 'f') || (b >= 'A' && b <= 'F')
}

func isIdentStart(r rune) bool {
	return unicode.IsLetter(r) || r == '_'
}

func isIdentPart(r rune) bool {
	return unicode.IsLetter(r) || unicode.IsDigit(r) || r == '_'
}

type ASTNode interface {
	Format() string
}

type NumNode struct {
	Val float64
	Raw string
}

func (n NumNode) Format() string {
	return formatFloatClean(n.Val)
}

type ImagNode struct {
	Val float64
}

func (n ImagNode) Format() string {
	if n.Val == 1 {
		return "i"
	}
	if n.Val == -1 {
		return "-i"
	}
	return formatFloatClean(n.Val) + "i"
}

type IdentNode struct {
	Name string
}

func (n IdentNode) Format() string {
	return n.Name
}

type HistoryNode struct {
	Index int
}

func (n HistoryNode) Format() string {
	return fmt.Sprintf("$%d", n.Index)
}

type UnaryNode struct {
	Op string
	X  ASTNode
}

func (n UnaryNode) Format() string {
	switch inner := n.X.(type) {
	case BinaryNode:
		if inner.Op == "+" || inner.Op == "-" {
			return fmt.Sprintf("%s(%s)", n.Op, inner.Format())
		}
	}
	return n.Op + n.X.Format()
}

type BinaryNode struct {
	Op    string
	Left  ASTNode
	Right ASTNode
}

func binaryPrec(op string) int {
	switch op {
	case "==", "!=", "<", "<=", ">", ">=", "=":
		return 1
	case "+", "-":
		return 2
	case "*", "/", "%":
		return 3
	case "^":
		return 4
	default:
		return 0
	}
}

func (n BinaryNode) Format() string {
	myP := binaryPrec(n.Op)
	lStr := n.Left.Format()
	rStr := n.Right.Format()

	if lb, ok := n.Left.(BinaryNode); ok {
		lp := binaryPrec(lb.Op)
		if lp < myP || (n.Op == "^" && lp == myP) {
			lStr = "(" + lStr + ")"
		}
	} else if _, ok := n.Left.(UnaryNode); ok && n.Op == "^" {
		lStr = "(" + lStr + ")"
	}

	if rb, ok := n.Right.(BinaryNode); ok {
		rp := binaryPrec(rb.Op)
		if rp < myP || ((n.Op == "-" || n.Op == "/" || n.Op == "%") && rp == myP) {
			rStr = "(" + rStr + ")"
		}
	}

	if n.Op == "+" || n.Op == "-" {
		return fmt.Sprintf("%s %s %s", lStr, n.Op, rStr)
	}
	return fmt.Sprintf("%s%s%s", lStr, n.Op, rStr)
}

type PostfixNode struct {
	Op string
	X  ASTNode
}

func (n PostfixNode) Format() string {
	switch n.X.(type) {
	case BinaryNode, UnaryNode:
		return fmt.Sprintf("(%s)%s", n.X.Format(), n.Op)
	default:
		return n.X.Format() + n.Op
	}
}

type CallNode struct {
	Name       string
	PrimeCount int
	Args       []ASTNode
}

func (n CallNode) Format() string {
	primes := strings.Repeat("'", n.PrimeCount)
	args := make([]string, len(n.Args))
	for i, a := range n.Args {
		args[i] = a.Format()
	}
	return fmt.Sprintf("%s%s(%s)", n.Name, primes, strings.Join(args, ", "))
}

type VectorNode struct {
	Elems []ASTNode
}

func (n VectorNode) Format() string {
	parts := make([]string, len(n.Elems))
	for i, e := range n.Elems {
		parts[i] = e.Format()
	}
	return "[" + strings.Join(parts, ", ") + "]"
}

type MatrixNode struct {
	Rows [][]ASTNode
}

func (n MatrixNode) Format() string {
	rows := make([]string, len(n.Rows))
	for r, row := range n.Rows {
		cols := make([]string, len(row))
		for c, cell := range row {
			cols[c] = cell.Format()
		}
		rows[r] = strings.Join(cols, ", ")
	}
	return "[" + strings.Join(rows, "; ") + "]"
}

type AssignVarNode struct {
	Name string
	Expr ASTNode
}

func (n AssignVarNode) Format() string {
	return fmt.Sprintf("%s = %s", n.Name, n.Expr.Format())
}

type AssignFuncNode struct {
	Name   string
	Params []string
	Body   ASTNode
}

func (n AssignFuncNode) Format() string {
	return fmt.Sprintf("%s(%s) = %s", n.Name, strings.Join(n.Params, ", "), n.Body.Format())
}

type mathParser struct {
	tokens []Token
	pos    int
}

func parseMath(input string) (ASTNode, error) {
	tokens, err := tokenize(input)
	if err != nil {
		return nil, err
	}
	p := &mathParser{tokens: tokens, pos: 0}
	node, err := p.parseTopLevel()
	if err != nil {
		return nil, err
	}
	if p.peek().Kind != TokEOF {
		return nil, fmt.Errorf("unexpected token %q", p.peek().Text)
	}
	return node, nil
}

func (p *mathParser) peek() Token {
	if p.pos >= len(p.tokens) {
		return Token{Kind: TokEOF}
	}
	return p.tokens[p.pos]
}

func (p *mathParser) next() Token {
	t := p.peek()
	if p.pos < len(p.tokens) {
		p.pos++
	}
	return t
}

func (p *mathParser) parseTopLevel() (ASTNode, error) {
	if p.peek().Kind == TokIdent {
		saved := p.pos
		nameTok := p.next()
		if p.peek().Kind == TokAssign {
			p.next()
			expr, err := p.parseCompare()
			if err != nil {
				return nil, err
			}
			return AssignVarNode{Name: nameTok.Text, Expr: expr}, nil
		}
		if p.peek().Kind == TokLParen {
			p.next()
			var params []string
			validParams := true
			if p.peek().Kind != TokRParen {
				for {
					if p.peek().Kind != TokIdent {
						validParams = false
						break
					}
					params = append(params, p.next().Text)
					if p.peek().Kind == TokComma {
						p.next()
					} else {
						break
					}
				}
			}
			if validParams && len(params) > 0 && p.peek().Kind == TokRParen {
				p.next()
				if p.peek().Kind == TokAssign {
					p.next()
					body, err := p.parseCompare()
					if err != nil {
						return nil, err
					}
					return AssignFuncNode{Name: nameTok.Text, Params: params, Body: body}, nil
				}
			}
		}
		p.pos = saved
	}
	return p.parseCompare()
}

func (p *mathParser) parseCompare() (ASTNode, error) {
	left, err := p.parseAdd()
	if err != nil {
		return nil, err
	}
	for {
		k := p.peek().Kind
		if k == TokEq || k == TokNeq || k == TokLt || k == TokLte || k == TokGt || k == TokGte || k == TokAssign {
			op := p.next().Text
			right, err := p.parseAdd()
			if err != nil {
				return nil, err
			}
			left = BinaryNode{Op: op, Left: left, Right: right}
		} else {
			break
		}
	}
	return left, nil
}

func (p *mathParser) parseAdd() (ASTNode, error) {
	left, err := p.parseMul()
	if err != nil {
		return nil, err
	}
	for {
		k := p.peek().Kind
		if k == TokPlus || k == TokMinus {
			op := p.next().Text
			right, err := p.parseMul()
			if err != nil {
				return nil, err
			}
			left = BinaryNode{Op: op, Left: left, Right: right}
		} else {
			break
		}
	}
	return left, nil
}

func (p *mathParser) canImplicitMul() bool {
	if p.pos == 0 {
		return false
	}
	prev := p.tokens[p.pos-1].Kind
	prevValid := prev == TokNumber || prev == TokImag || prev == TokIdent || prev == TokHistoryRef ||
		prev == TokRParen || prev == TokRBracket || prev == TokBang || prev == TokDoubleBang || prev == TokPrime
	if !prevValid {
		return false
	}
	nxt := p.peek().Kind
	return nxt == TokNumber || nxt == TokImag || nxt == TokIdent || nxt == TokHistoryRef ||
		nxt == TokLParen || nxt == TokLBracket
}

func (p *mathParser) parseMul() (ASTNode, error) {
	left, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	for {
		k := p.peek().Kind
		if k == TokStar || k == TokSlash || k == TokPercent {
			op := p.next().Text
			right, err := p.parseUnary()
			if err != nil {
				return nil, err
			}
			left = BinaryNode{Op: op, Left: left, Right: right}
		} else if p.canImplicitMul() {
			right, err := p.parseUnary()
			if err != nil {
				return nil, err
			}
			left = BinaryNode{Op: "*", Left: left, Right: right}
		} else {
			break
		}
	}
	return left, nil
}

func (p *mathParser) parseUnary() (ASTNode, error) {
	k := p.peek().Kind
	if k == TokMinus {
		p.next()
		sub, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return UnaryNode{Op: "-", X: sub}, nil
	}
	if k == TokPlus {
		p.next()
		return p.parseUnary()
	}
	return p.parsePow()
}

func (p *mathParser) parsePow() (ASTNode, error) {
	base, err := p.parsePostfix()
	if err != nil {
		return nil, err
	}
	if p.peek().Kind == TokCaret {
		p.next()
		exp, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return BinaryNode{Op: "^", Left: base, Right: exp}, nil
	}
	return base, nil
}

func (p *mathParser) parsePostfix() (ASTNode, error) {
	node, err := p.parsePrimary()
	if err != nil {
		return nil, err
	}
	for {
		k := p.peek().Kind
		if k == TokBang {
			p.next()
			node = PostfixNode{Op: "!", X: node}
		} else if k == TokDoubleBang {
			p.next()
			node = PostfixNode{Op: "!!", X: node}
		} else if k == TokPrime {
			p.next()
			node = PostfixNode{Op: "'", X: node}
		} else {
			break
		}
	}
	return node, nil
}

func (p *mathParser) parsePrimary() (ASTNode, error) {
	tok := p.peek()
	switch tok.Kind {
	case TokNumber:
		p.next()
		return NumNode{Val: tok.Num, Raw: tok.Text}, nil

	case TokImag:
		p.next()
		return ImagNode{Val: tok.Num}, nil

	case TokHistoryRef:
		p.next()
		idx, err := strconv.Atoi(tok.Text)
		if err != nil {
			return nil, err
		}
		return HistoryNode{Index: idx}, nil

	case TokIdent:
		p.next()
		primes := 0
		for p.peek().Kind == TokPrime {
			p.next()
			primes++
		}
		if p.peek().Kind == TokLParen {
			p.next()
			var args []ASTNode
			if p.peek().Kind != TokRParen {
				for {
					arg, err := p.parseCompare()
					if err != nil {
						return nil, err
					}
					args = append(args, arg)
					if p.peek().Kind == TokComma {
						p.next()
					} else {
						break
					}
				}
			}
			if p.peek().Kind != TokRParen {
				return nil, fmt.Errorf("missing ')' after arguments to %s", tok.Text)
			}
			p.next()
			return CallNode{Name: tok.Text, PrimeCount: primes, Args: args}, nil
		}
		var out ASTNode = IdentNode{Name: tok.Text}
		for i := 0; i < primes; i++ {
			out = PostfixNode{Op: "'", X: out}
		}
		return out, nil

	case TokLParen:
		p.next()
		inner, err := p.parseCompare()
		if err != nil {
			return nil, err
		}
		if p.peek().Kind != TokRParen {
			return nil, fmt.Errorf("missing closing ')'")
		}
		p.next()
		return inner, nil

	case TokLBracket:
		return p.parseBracketLiteral()

	default:
		return nil, fmt.Errorf("unexpected token %q", tok.Text)
	}
}

func (p *mathParser) parseBracketLiteral() (ASTNode, error) {
	p.next()
	if p.peek().Kind == TokRBracket {
		p.next()
		return VectorNode{Elems: nil}, nil
	}

	if p.peek().Kind == TokLBracket {
		var rows [][]ASTNode
		for {
			if p.peek().Kind != TokLBracket {
				return nil, fmt.Errorf("expected '[' for matrix row")
			}
			p.next()
			var row []ASTNode
			for p.peek().Kind != TokRBracket {
				el, err := p.parseCompare()
				if err != nil {
					return nil, err
				}
				row = append(row, el)
				if p.peek().Kind == TokComma {
					p.next()
				} else {
					break
				}
			}
			if p.peek().Kind != TokRBracket {
				return nil, fmt.Errorf("missing ']' at end of matrix row")
			}
			p.next()
			rows = append(rows, row)
			if p.peek().Kind == TokComma || p.peek().Kind == TokSemi {
				p.next()
			} else {
				break
			}
		}
		if p.peek().Kind != TokRBracket {
			return nil, fmt.Errorf("missing closing ']' for matrix")
		}
		p.next()
		return MatrixNode{Rows: rows}, nil
	}

	var rows [][]ASTNode
	var curRow []ASTNode
	hasSemi := false
	for {
		el, err := p.parseCompare()
		if err != nil {
			return nil, err
		}
		curRow = append(curRow, el)
		if p.peek().Kind == TokComma {
			p.next()
		} else if p.peek().Kind == TokSemi {
			hasSemi = true
			p.next()
			rows = append(rows, curRow)
			curRow = nil
			if p.peek().Kind == TokRBracket {
				break
			}
		} else {
			break
		}
	}
	if p.peek().Kind != TokRBracket {
		return nil, fmt.Errorf("missing closing ']'")
	}
	p.next()
	if hasSemi {
		if len(curRow) > 0 {
			rows = append(rows, curRow)
		}
		return MatrixNode{Rows: rows}, nil
	}
	return VectorNode{Elems: curRow}, nil
}

func (e *Engine) Eval(input string) (Value, error) {
	node, err := parseMath(input)
	if err != nil {
		return Value{}, err
	}
	val, err := e.evalAST(node, nil, 0)
	if err != nil {
		return Value{}, err
	}
	e.LastValue = val
	e.History = append(e.History, HistoryEntry{Input: input, Val: val})
	return val, nil
}

func (e *Engine) resolveIdent(name string, scope map[string]Value) (Value, bool) {
	if scope != nil {
		if v, ok := scope[name]; ok {
			return v, true
		}
	}
	if name == "ans" || name == "_" {
		return e.LastValue, true
	}
	if v, ok := e.Vars[name]; ok {
		return v, true
	}
	lower := strings.ToLower(name)
	if z, ok := builtinConstants[lower]; ok {
		return ScalarComplex(z), true
	}
	return Value{}, false
}

var reservedConstantNames = map[string]bool{
	"pi": true, "tau": true, "e": true, "i": true, "inf": true, "infinity": true, "ans": true,
}

func (e *Engine) evalAST(node ASTNode, scope map[string]Value, depth int) (Value, error) {
	if depth > 256 {
		return Value{}, fmt.Errorf("maximum recursion depth exceeded")
	}

	switch n := node.(type) {
	case AssignVarNode:
		if reservedConstantNames[strings.ToLower(n.Name)] {
			return Value{}, fmt.Errorf("cannot overwrite constant %s", n.Name)
		}
		val, err := e.evalAST(n.Expr, scope, depth+1)
		if err != nil {
			return Value{}, err
		}
		e.Vars[n.Name] = val
		return val, nil

	case AssignFuncNode:
		e.Funcs[n.Name] = UserFunc{
			Params: n.Params,
			Body:   n.Body,
			Raw:    n.Format(),
		}
		return SymbolicVal(n.Format(), n.Body), nil

	case NumNode:
		return ScalarReal(n.Val), nil

	case ImagNode:
		return ScalarComplex(complex(0, n.Val)), nil

	case HistoryNode:
		if n.Index < 1 || n.Index > len(e.History) {
			return Value{}, fmt.Errorf("history entry $%d out of range", n.Index)
		}
		return e.History[n.Index-1].Val, nil

	case IdentNode:
		if val, ok := e.resolveIdent(n.Name, scope); ok {
			if val.Kind == KindSymbolic && val.Node != nil && scope != nil && len(scope) > 0 {
				return e.evalAST(val.Node, scope, depth+1)
			}
			return val, nil
		}
		return Value{}, fmt.Errorf("unknown identifier %s", n.Name)

	case VectorNode:
		elems := make([]complex128, len(n.Elems))
		for idx, el := range n.Elems {
			v, err := e.evalAST(el, scope, depth+1)
			if err != nil {
				return Value{}, err
			}
			if v.Kind != KindScalar {
				return Value{}, fmt.Errorf("vector element %d must be scalar", idx)
			}
			elems[idx] = v.S
		}
		return VectorVal(elems), nil

	case MatrixNode:
		if len(n.Rows) == 0 {
			return MatrixVal(nil), nil
		}
		cols := len(n.Rows[0])
		mat := make([][]complex128, len(n.Rows))
		for r, row := range n.Rows {
			if len(row) != cols {
				return Value{}, fmt.Errorf("matrix row %d has length %d, expected %d", r, len(row), cols)
			}
			mat[r] = make([]complex128, cols)
			for c, cell := range row {
				v, err := e.evalAST(cell, scope, depth+1)
				if err != nil {
					return Value{}, err
				}
				if v.Kind != KindScalar {
					return Value{}, fmt.Errorf("matrix entry [%d,%d] must be scalar", r, c)
				}
				mat[r][c] = v.S
			}
		}
		return MatrixVal(mat), nil

	case UnaryNode:
		x, err := e.evalAST(n.X, scope, depth+1)
		if err != nil {
			return Value{}, err
		}
		if n.Op == "-" {
			return negateValue(x)
		}
		return x, nil

	case PostfixNode:
		x, err := e.evalAST(n.X, scope, depth+1)
		if err != nil {
			return Value{}, err
		}
		switch n.Op {
		case "!":
			if x.Kind != KindScalar {
				return Value{}, fmt.Errorf("factorial requires scalar")
			}
			if x.IsReal() {
				return ScalarReal(factorialReal(real(x.S))), nil
			}
			return ScalarComplex(gammaComplex(x.S + 1)), nil
		case "!!":
			r, err := x.Real()
			if err != nil {
				return Value{}, err
			}
			return ScalarReal(doubleFactorial(r)), nil
		case "'":
			switch x.Kind {
			case KindMatrix:
				return MatrixVal(matConjTranspose(x.M)), nil
			case KindVector:
				col := make([][]complex128, len(x.V))
				for i, z := range x.V {
					col[i] = []complex128{cmplx.Conj(z)}
				}
				return MatrixVal(col), nil
			case KindSymbolic:
				if x.Node != nil {
					dNode, err := symbolicDiff(x.Node, "x")
					if err != nil {
						return Value{}, err
					}
					simp := simplifyAST(dNode)
					return SymbolicVal(simp.Format(), simp), nil
				}
			}
			return Value{}, fmt.Errorf("transpose/prime not supported on scalar")
		}
		return Value{}, fmt.Errorf("unknown postfix %s", n.Op)

	case BinaryNode:
		left, err := e.evalAST(n.Left, scope, depth+1)
		if err != nil {
			return Value{}, err
		}
		right, err := e.evalAST(n.Right, scope, depth+1)
		if err != nil {
			return Value{}, err
		}
		return evalBinaryValues(n.Op, left, right)

	case CallNode:
		if n.PrimeCount > 0 {
			return e.evalPrimedCall(n, scope, depth+1)
		}
		if handled, val, err := e.evalSpecialForm(n, scope, depth+1); handled {
			return val, err
		}
		if uf, ok := e.Funcs[n.Name]; ok {
			if len(n.Args) != len(uf.Params) {
				return Value{}, fmt.Errorf("%s expects %d arguments, got %d", n.Name, len(uf.Params), len(n.Args))
			}
			localScope := make(map[string]Value, len(uf.Params))
			if scope != nil {
				for k, v := range scope {
					localScope[k] = v
				}
			}
			for i, pName := range uf.Params {
				argVal, err := e.evalAST(n.Args[i], scope, depth+1)
				if err != nil {
					return Value{}, err
				}
				localScope[pName] = argVal
			}
			return e.evalAST(uf.Body, localScope, depth+1)
		}
		if handled, val, err := e.evalBuiltinFunc(n, scope, depth+1); handled {
			return val, err
		}
		if len(n.Args) == 1 {
			if identVal, ok := e.resolveIdent(n.Name, scope); ok {
				argVal, err := e.evalAST(n.Args[0], scope, depth+1)
				if err != nil {
					return Value{}, err
				}
				return evalBinaryValues("*", identVal, argVal)
			}
		}
		return Value{}, fmt.Errorf("unknown function %s", n.Name)

	default:
		return Value{}, fmt.Errorf("unsupported expression")
	}
}

func (e *Engine) evalPrimedCall(n CallNode, scope map[string]Value, depth int) (Value, error) {
	if len(n.Args) != 1 {
		return Value{}, fmt.Errorf("prime derivative notation requires 1 argument")
	}
	var baseExpr ASTNode
	paramName := "x"
	if uf, ok := e.Funcs[n.Name]; ok {
		if len(uf.Params) != 1 {
			return Value{}, fmt.Errorf("prime notation requires single-variable function")
		}
		paramName = uf.Params[0]
		baseExpr = uf.Body
	} else {
		baseExpr = CallNode{Name: n.Name, Args: []ASTNode{IdentNode{Name: paramName}}}
	}

	cur := baseExpr
	var err error
	for i := 0; i < n.PrimeCount; i++ {
		cur, err = symbolicDiff(cur, paramName)
		if err != nil {
			return Value{}, err
		}
		cur = simplifyAST(cur)
	}

	if id, ok := n.Args[0].(IdentNode); ok {
		if _, bound := e.resolveIdent(id.Name, scope); !bound {
			sub := substituteAST(cur, paramName, id)
			simp := simplifyAST(sub)
			return SymbolicVal(simp.Format(), simp), nil
		}
	}

	argVal, err := e.evalAST(n.Args[0], scope, depth+1)
	if err != nil {
		return Value{}, err
	}
	localScope := make(map[string]Value)
	if scope != nil {
		for k, v := range scope {
			localScope[k] = v
		}
	}
	localScope[paramName] = argVal
	return e.evalAST(cur, localScope, depth+1)
}

func negateValue(v Value) (Value, error) {
	switch v.Kind {
	case KindScalar:
		return ScalarComplex(-v.S), nil
	case KindVector:
		out := make([]complex128, len(v.V))
		for i, z := range v.V {
			out[i] = -z
		}
		return VectorVal(out), nil
	case KindMatrix:
		out := make([][]complex128, len(v.M))
		for r := range v.M {
			out[r] = make([]complex128, len(v.M[r]))
			for c, z := range v.M[r] {
				out[r][c] = -z
			}
		}
		return MatrixVal(out), nil
	case KindSymbolic:
		if v.Node != nil {
			simp := simplifyAST(UnaryNode{Op: "-", X: v.Node})
			return SymbolicVal(simp.Format(), simp), nil
		}
	}
	return Value{}, fmt.Errorf("cannot negate %v", v.Kind)
}

func evalBinaryValues(op string, a, b Value) (Value, error) {
	if a.Kind == KindScalar && b.Kind == KindScalar {
		return evalScalarBinary(op, a, b)
	}

	switch op {
	case "+":
		if a.Kind == KindVector && b.Kind == KindVector {
			if len(a.V) != len(b.V) {
				return Value{}, fmt.Errorf("vector length mismatch (%d vs %d)", len(a.V), len(b.V))
			}
			out := make([]complex128, len(a.V))
			for i := range a.V {
				out[i] = a.V[i] + b.V[i]
			}
			return VectorVal(out), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindMatrix {
			m, err := matAdd(a.M, b.M, 1)
			if err != nil {
				return Value{}, err
			}
			return MatrixVal(m), nil
		}

	case "-":
		if a.Kind == KindVector && b.Kind == KindVector {
			if len(a.V) != len(b.V) {
				return Value{}, fmt.Errorf("vector length mismatch (%d vs %d)", len(a.V), len(b.V))
			}
			out := make([]complex128, len(a.V))
			for i := range a.V {
				out[i] = a.V[i] - b.V[i]
			}
			return VectorVal(out), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindMatrix {
			m, err := matAdd(a.M, b.M, -1)
			if err != nil {
				return Value{}, err
			}
			return MatrixVal(m), nil
		}

	case "*":
		if a.Kind == KindScalar && b.Kind == KindVector {
			out := make([]complex128, len(b.V))
			for i, z := range b.V {
				out[i] = a.S * z
			}
			return VectorVal(out), nil
		}
		if a.Kind == KindVector && b.Kind == KindScalar {
			out := make([]complex128, len(a.V))
			for i, z := range a.V {
				out[i] = z * b.S
			}
			return VectorVal(out), nil
		}
		if a.Kind == KindScalar && b.Kind == KindMatrix {
			return MatrixVal(matScale(b.M, a.S)), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindScalar {
			return MatrixVal(matScale(a.M, b.S)), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindMatrix {
			m, err := matMul(a.M, b.M)
			if err != nil {
				return Value{}, err
			}
			return MatrixVal(m), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindVector {
			v, err := matVecMul(a.M, b.V)
			if err != nil {
				return Value{}, err
			}
			return VectorVal(v), nil
		}
		if a.Kind == KindVector && b.Kind == KindMatrix {
			v, err := vecMatMul(a.V, b.M)
			if err != nil {
				return Value{}, err
			}
			return VectorVal(v), nil
		}

	case "/":
		if a.Kind == KindVector && b.Kind == KindScalar {
			out := make([]complex128, len(a.V))
			for i, z := range a.V {
				out[i] = z / b.S
			}
			return VectorVal(out), nil
		}
		if a.Kind == KindMatrix && b.Kind == KindScalar {
			return MatrixVal(matScale(a.M, 1/b.S)), nil
		}

	case "^":
		if a.Kind == KindMatrix && b.Kind == KindScalar {
			expReal, err := b.Real()
			if err != nil || expReal != math.Trunc(expReal) {
				return Value{}, fmt.Errorf("matrix power requires integer exponent")
			}
			m, err := matPow(a.M, int(expReal))
			if err != nil {
				return Value{}, err
			}
			return MatrixVal(m), nil
		}
	}

	return Value{}, fmt.Errorf("unsupported operation %s between types", op)
}

func evalScalarBinary(op string, a, b Value) (Value, error) {
	switch op {
	case "+":
		return ScalarComplex(a.S + b.S), nil
	case "-":
		return ScalarComplex(a.S - b.S), nil
	case "*":
		return ScalarComplex(a.S * b.S), nil
	case "/":
		if b.S == 0 {
			if a.IsReal() {
				if real(a.S) > 0 {
					return ScalarReal(math.Inf(1)), nil
				} else if real(a.S) < 0 {
					return ScalarReal(math.Inf(-1)), nil
				}
				return ScalarReal(math.NaN()), nil
			}
			return ScalarReal(math.NaN()), nil
		}
		return ScalarComplex(a.S / b.S), nil
	case "%":
		ar, err1 := a.Real()
		br, err2 := b.Real()
		if err1 != nil || err2 != nil {
			return Value{}, fmt.Errorf("modulo requires real operands")
		}
		return ScalarReal(math.Mod(ar, br)), nil
	case "^":
		if a.S == 0 {
			if b.IsReal() {
				br := real(b.S)
				if br > 0 {
					return ScalarReal(0), nil
				}
				if br == 0 {
					return ScalarReal(1), nil
				}
				return ScalarReal(math.Inf(1)), nil
			}
		}
		if a.IsReal() && b.IsReal() {
			ar, br := real(a.S), real(b.S)
			if ar >= 0 || br == math.Trunc(br) {
				return ScalarReal(math.Pow(ar, br)), nil
			}
		}
		return ScalarComplex(cmplx.Pow(a.S, b.S)), nil
	case "==", "=":
		if cmplx.Abs(a.S-b.S) <= 1e-12*math.Max(1.0, math.Max(cmplx.Abs(a.S), cmplx.Abs(b.S))) {
			return ScalarReal(1), nil
		}
		return ScalarReal(0), nil
	case "!=":
		if cmplx.Abs(a.S-b.S) > 1e-12*math.Max(1.0, math.Max(cmplx.Abs(a.S), cmplx.Abs(b.S))) {
			return ScalarReal(1), nil
		}
		return ScalarReal(0), nil
	case "<", "<=", ">", ">=":
		ar, err1 := a.Real()
		br, err2 := b.Real()
		if err1 != nil || err2 != nil {
			return Value{}, fmt.Errorf("ordering comparison requires real numbers")
		}
		var ok bool
		switch op {
		case "<":
			ok = ar < br
		case "<=":
			ok = ar <= br
		case ">":
			ok = ar > br
		case ">=":
			ok = ar >= br
		}
		if ok {
			return ScalarReal(1), nil
		}
		return ScalarReal(0), nil
	default:
		return Value{}, fmt.Errorf("unknown scalar operator %s", op)
	}
}

func (e *Engine) toRad(z complex128) complex128 {
	if e.Mode == ModeDeg && math.Abs(imag(z)) == 0 {
		return complex(real(z)*math.Pi/180.0, 0)
	}
	return z
}

func (e *Engine) fromRad(z complex128) complex128 {
	z = cleanComplex(z)
	if e.Mode == ModeDeg && imag(z) == 0 {
		return complex(real(z)*180.0/math.Pi, 0)
	}
	return z
}

func cleanSin(z complex128) complex128 {
	if imag(z) == 0 {
		x := real(z)
		if !math.IsInf(x, 0) && !math.IsNaN(x) && x != 0 {
			q := x / math.Pi
			if math.Abs(q-math.Round(q)) < 1e-15 {
				return 0
			}
		}
		return complex(math.Sin(x), 0)
	}
	return cmplx.Sin(z)
}

func cleanCos(z complex128) complex128 {
	if imag(z) == 0 {
		x := real(z)
		if !math.IsInf(x, 0) && !math.IsNaN(x) {
			q := x/math.Pi - 0.5
			if math.Abs(q-math.Round(q)) < 1e-15 {
				return 0
			}
		}
		return complex(math.Cos(x), 0)
	}
	return cmplx.Cos(z)
}

func cleanTan(z complex128) complex128 {
	if imag(z) == 0 {
		x := real(z)
		if !math.IsInf(x, 0) && !math.IsNaN(x) {
			q := x/math.Pi - 0.5
			if math.Abs(q-math.Round(q)) < 1e-14 {
				return complex(math.NaN(), 0)
			}
			qZero := x / math.Pi
			if x != 0 && math.Abs(qZero-math.Round(qZero)) < 1e-15 {
				return 0
			}
		}
		return complex(math.Tan(x), 0)
	}
	return cmplx.Tan(z)
}

func (e *Engine) evalBuiltinFunc(n CallNode, scope map[string]Value, depth int) (bool, Value, error) {
	fn := strings.ToLower(n.Name)

	args := make([]Value, len(n.Args))
	for i, a := range n.Args {
		v, err := e.evalAST(a, scope, depth+1)
		if err != nil {
			return true, Value{}, err
		}
		args[i] = v
	}

	if handled, val, err := evalLinAlgFunc(fn, args); handled {
		return true, val, err
	}
	if handled, val, err := evalNumStatsFunc(fn, args); handled {
		return true, val, err
	}
	if handled, val, err := evalSpecialMathFunc(fn, args); handled {
		return true, val, err
	}

	if len(args) == 1 && args[0].Kind == KindScalar {
		z := args[0].S
		switch fn {
		case "sin":
			return true, ScalarComplex(cleanSin(e.toRad(z))), nil
		case "cos":
			return true, ScalarComplex(cleanCos(e.toRad(z))), nil
		case "tan":
			return true, ScalarComplex(cleanTan(e.toRad(z))), nil
		case "sec":
			return true, ScalarComplex(1 / cleanCos(e.toRad(z))), nil
		case "csc":
			return true, ScalarComplex(1 / cleanSin(e.toRad(z))), nil
		case "cot":
			return true, ScalarComplex(cleanCos(e.toRad(z)) / cleanSin(e.toRad(z))), nil
		case "asin", "arcsin":
			if args[0].IsReal() && math.Abs(real(z)) <= 1 {
				return true, ScalarComplex(e.fromRad(complex(math.Asin(real(z)), 0))), nil
			}
			return true, ScalarComplex(e.fromRad(cmplx.Asin(z))), nil
		case "acos", "arccos":
			if args[0].IsReal() && math.Abs(real(z)) <= 1 {
				return true, ScalarComplex(e.fromRad(complex(math.Acos(real(z)), 0))), nil
			}
			return true, ScalarComplex(e.fromRad(cmplx.Acos(z))), nil
		case "atan", "arctan":
			if args[0].IsReal() {
				return true, ScalarComplex(e.fromRad(complex(math.Atan(real(z)), 0))), nil
			}
			return true, ScalarComplex(e.fromRad(cmplx.Atan(z))), nil
		case "asec", "arcsec":
			return true, ScalarComplex(e.fromRad(cmplx.Acos(1 / z))), nil
		case "acsc", "arccsc":
			return true, ScalarComplex(e.fromRad(cmplx.Asin(1 / z))), nil
		case "acot", "arccot":
			return true, ScalarComplex(e.fromRad(cmplx.Atan(1 / z))), nil
		case "sinh":
			return true, ScalarComplex(cmplx.Sinh(z)), nil
		case "cosh":
			return true, ScalarComplex(cmplx.Cosh(z)), nil
		case "tanh":
			return true, ScalarComplex(cmplx.Tanh(z)), nil
		case "sech":
			return true, ScalarComplex(1 / cmplx.Cosh(z)), nil
		case "csch":
			return true, ScalarComplex(1 / cmplx.Sinh(z)), nil
		case "coth":
			return true, ScalarComplex(1 / cmplx.Tanh(z)), nil
		case "asinh", "arcsinh":
			return true, ScalarComplex(cmplx.Asinh(z)), nil
		case "acosh", "arccosh":
			return true, ScalarComplex(cmplx.Acosh(z)), nil
		case "atanh", "arctanh":
			return true, ScalarComplex(cmplx.Atanh(z)), nil
		case "sqrt":
			if args[0].IsReal() && real(z) >= 0 {
				return true, ScalarReal(math.Sqrt(real(z))), nil
			}
			return true, ScalarComplex(cmplx.Sqrt(z)), nil
		case "cbrt":
			if args[0].IsReal() {
				return true, ScalarReal(math.Cbrt(real(z))), nil
			}
			return true, ScalarComplex(cmplx.Pow(z, complex(1.0/3.0, 0))), nil
		case "exp":
			if args[0].IsReal() {
				return true, ScalarReal(math.Exp(real(z))), nil
			}
			return true, ScalarComplex(cmplx.Exp(z)), nil
		case "ln":
			if args[0].IsReal() && real(z) > 0 {
				return true, ScalarReal(math.Log(real(z))), nil
			}
			if z == 0 {
				return true, ScalarReal(math.Inf(-1)), nil
			}
			return true, ScalarComplex(cmplx.Log(z)), nil
		case "log", "log10":
			if args[0].IsReal() && real(z) > 0 {
				return true, ScalarReal(math.Log10(real(z))), nil
			}
			if z == 0 {
				return true, ScalarReal(math.Inf(-1)), nil
			}
			return true, ScalarComplex(cmplx.Log10(z)), nil
		case "log2":
			if args[0].IsReal() && real(z) > 0 {
				return true, ScalarReal(math.Log2(real(z))), nil
			}
			if z == 0 {
				return true, ScalarReal(math.Inf(-1)), nil
			}
			return true, ScalarComplex(cmplx.Log(z) / complex(math.Ln2, 0)), nil
		case "abs":
			return true, ScalarReal(cmplx.Abs(z)), nil
		case "re", "real":
			return true, ScalarReal(real(z)), nil
		case "im", "imag":
			return true, ScalarReal(imag(z)), nil
		case "conj":
			return true, ScalarComplex(cmplx.Conj(z)), nil
		case "arg", "phase":
			p := cmplx.Phase(z)
			if e.Mode == ModeDeg {
				p = p * 180.0 / math.Pi
			}
			return true, ScalarReal(p), nil
		case "polar":
			r, theta := cmplx.Polar(z)
			if e.Mode == ModeDeg {
				theta = theta * 180.0 / math.Pi
			}
			return true, VectorVal([]complex128{complex(r, 0), complex(theta, 0)}), nil
		case "floor":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(math.Floor(r)), nil
		case "ceil":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(math.Ceil(r)), nil
		case "round":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(math.Round(r)), nil
		case "trunc":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(math.Trunc(r)), nil
		case "sgn", "sign":
			r, err := args[0].Real()
			if err != nil {
				if z == 0 {
					return true, ScalarReal(0), nil
				}
				return true, ScalarComplex(z / complex(cmplx.Abs(z), 0)), nil
			}
			if r > 0 {
				return true, ScalarReal(1), nil
			} else if r < 0 {
				return true, ScalarReal(-1), nil
			}
			return true, ScalarReal(0), nil
		case "heaviside", "step":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			if r < 0 {
				return true, ScalarReal(0), nil
			} else if r == 0 {
				return true, ScalarReal(0.5), nil
			}
			return true, ScalarReal(1), nil
		case "deg2rad", "rad":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(r * math.Pi / 180.0), nil
		case "rad2deg", "deg":
			r, err := args[0].Real()
			if err != nil {
				return true, Value{}, err
			}
			return true, ScalarReal(r * 180.0 / math.Pi), nil
		case "fact", "factorial":
			if args[0].IsReal() {
				return true, ScalarReal(factorialReal(real(z))), nil
			}
			return true, ScalarComplex(gammaComplex(z + 1)), nil
		}
	}

	if len(args) == 2 && args[0].Kind == KindScalar && args[1].Kind == KindScalar {
		z1, z2 := args[0].S, args[1].S
		switch fn {
		case "pow":
			v, err := evalScalarBinary("^", args[0], args[1])
			return true, v, err
		case "root":
			if args[0].IsReal() && args[1].IsReal() {
				nDeg := real(z1)
				val := real(z2)
				if val >= 0 {
					return true, ScalarReal(math.Pow(val, 1.0/nDeg)), nil
				}
				if nDeg == math.Trunc(nDeg) && int64(math.Abs(nDeg))%2 == 1 {
					return true, ScalarReal(-math.Pow(-val, 1.0/nDeg)), nil
				}
			}
			return true, ScalarComplex(cmplx.Pow(z2, 1/z1)), nil
		case "logb":
			return true, ScalarComplex(cmplx.Log(z2) / cmplx.Log(z1)), nil
		case "atan2":
			y, err1 := args[0].Real()
			x, err2 := args[1].Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("atan2 requires real arguments")
			}
			ang := math.Atan2(y, x)
			if e.Mode == ModeDeg {
				ang = ang * 180.0 / math.Pi
			}
			return true, ScalarReal(ang), nil
		case "hypot":
			a, err1 := args[0].Real()
			b, err2 := args[1].Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("hypot requires real arguments")
			}
			return true, ScalarReal(math.Hypot(a, b)), nil
		case "rect":
			r, err1 := args[0].Real()
			theta, err2 := args[1].Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("rect requires real r and theta")
			}
			if e.Mode == ModeDeg {
				theta = theta * math.Pi / 180.0
			}
			return true, ScalarComplex(cmplx.Rect(r, theta)), nil
		case "mod":
			a, err1 := args[0].Real()
			b, err2 := args[1].Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("mod requires real arguments")
			}
			m := math.Mod(a, b)
			if m < 0 && b > 0 {
				m += b
			}
			return true, ScalarReal(m), nil
		case "round":
			r, err1 := args[0].Real()
			places, err2 := args[1].Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("round requires real arguments")
			}
			scale := math.Pow(10, math.Round(places))
			return true, ScalarReal(math.Round(r*scale) / scale), nil
		}
	}

	return false, Value{}, nil
}

func (e *Engine) SortedVarNames() []string {
	keys := make([]string, 0, len(e.Vars))
	for k := range e.Vars {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

func (e *Engine) SortedFuncNames() []string {
	keys := make([]string, 0, len(e.Funcs))
	for k := range e.Funcs {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

func EvaluateExpression(expr string) float64 {
	val, err := defaultEngine.Eval(expr)
	if err != nil {
		return 0
	}
	if r, err := val.Real(); err == nil {
		return r
	}
	return 0
}

func EvaluateWithScope(expr string, scope map[string]float64) (float64, error) {
	node, err := parseMath(expr)
	if err != nil {
		return 0, err
	}
	valScope := make(map[string]Value, len(scope))
	for k, v := range scope {
		valScope[k] = ScalarReal(v)
	}
	val, err := defaultEngine.evalAST(node, valScope, 0)
	if err != nil {
		return 0, err
	}
	return val.Real()
}
