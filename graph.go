package main

import (
	"fmt"
	"math"
	"sort"
	"strings"
)

var brailleBits = [4][2]uint8{
	{0x01, 0x08},
	{0x02, 0x10},
	{0x04, 0x20},
	{0x40, 0x80},
}

var seriesColors = []string{
	"\033[36m",
	"\033[33m",
	"\033[35m",
	"\033[32m",
	"\033[31m",
}

const ansiReset = "\033[0m"
const ansiDim = "\033[90m"

type brailleCanvas struct {
	charW int
	charH int
	pixW  int
	pixH  int
	dots  [][]uint8
	owner [][]int
}

func newBrailleCanvas(charW, charH int) *brailleCanvas {
	dots := make([][]uint8, charH)
	owner := make([][]int, charH)
	for r := 0; r < charH; r++ {
		dots[r] = make([]uint8, charW)
		owner[r] = make([]int, charW)
		for c := 0; c < charW; c++ {
			owner[r][c] = -1
		}
	}
	return &brailleCanvas{
		charW: charW,
		charH: charH,
		pixW:  charW * 2,
		pixH:  charH * 4,
		dots:  dots,
		owner: owner,
	}
}

func (bc *brailleCanvas) setPixel(px, py, series int) {
	if px < 0 || px >= bc.pixW || py < 0 || py >= bc.pixH {
		return
	}
	cx, cy := px/2, py/4
	lx, ly := px%2, py%4
	bc.dots[cy][cx] |= brailleBits[ly][lx]
	if series >= 0 || bc.owner[cy][cx] < 0 {
		bc.owner[cy][cx] = series
	}
}

func (bc *brailleCanvas) drawLine(x0, y0, x1, y1, series int) {
	dx := absInt(x1 - x0)
	sx := -1
	if x0 < x1 {
		sx = 1
	}
	dy := -absInt(y1 - y0)
	sy := -1
	if y0 < y1 {
		sy = 1
	}
	err := dx + dy
	steps := 0
	maxSteps := bc.pixW + bc.pixH + 10
	for steps < maxSteps {
		bc.setPixel(x0, y0, series)
		if x0 == x1 && y0 == y1 {
			break
		}
		e2 := 2 * err
		if e2 >= dy {
			err += dy
			x0 += sx
		}
		if e2 <= dx {
			err += dx
			y0 += sy
		}
		steps++
	}
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func (bc *brailleCanvas) renderWithMargin(yMax, yMin float64, useColor bool) string {
	const labelW = 10
	var sb strings.Builder
	midRow := bc.charH / 2
	yMid := 0.5 * (yMax + yMin)

	topBorder := strings.Repeat(" ", labelW) + " +" + strings.Repeat("-", bc.charW) + "+\n"
	sb.WriteString(topBorder)

	for r := 0; r < bc.charH; r++ {
		lbl := ""
		if r == 0 {
			lbl = fmt.Sprintf("%*.4g", labelW, yMax)
		} else if r == midRow {
			lbl = fmt.Sprintf("%*.4g", labelW, yMid)
		} else if r == bc.charH-1 {
			lbl = fmt.Sprintf("%*.4g", labelW, yMin)
		} else {
			lbl = strings.Repeat(" ", labelW)
		}
		if len(lbl) > labelW {
			lbl = lbl[:labelW]
		}
		sb.WriteString(lbl + " |")
		for c := 0; c < bc.charW; c++ {
			mask := bc.dots[r][c]
			if mask == 0 {
				sb.WriteRune(' ')
				continue
			}
			ch := rune(0x2800 + int(mask))
			own := bc.owner[r][c]
			if useColor {
				if own == -2 {
					sb.WriteString(ansiDim + string(ch) + ansiReset)
				} else if own >= 0 {
					col := seriesColors[own%len(seriesColors)]
					sb.WriteString(col + string(ch) + ansiReset)
				} else {
					sb.WriteRune(ch)
				}
			} else {
				sb.WriteRune(ch)
			}
		}
		sb.WriteString("|\n")
	}

	botBorder := strings.Repeat(" ", labelW) + " +" + strings.Repeat("-", bc.charW) + "+\n"
	sb.WriteString(botBorder)
	return sb.String()
}

func formatXAxisFooter(xMin, xMax float64, charW int) string {
	const labelW = 12
	xMid := 0.5 * (xMin + xMax)
	sMin := fmt.Sprintf("%.4g", xMin)
	sMid := fmt.Sprintf("%.4g", xMid)
	sMax := fmt.Sprintf("%.4g", xMax)

	line := make([]byte, charW)
	for i := range line {
		line[i] = ' '
	}
	for i := 0; i < len(sMin) && i < charW; i++ {
		line[i] = sMin[i]
	}
	midStart := (charW - len(sMid)) / 2
	if midStart > len(sMin)+1 {
		for i := 0; i < len(sMid) && midStart+i < charW; i++ {
			line[midStart+i] = sMid[i]
		}
	}
	maxStart := charW - len(sMax)
	if maxStart > midStart+len(sMid)+1 {
		for i := 0; i < len(sMax) && maxStart+i < charW; i++ {
			line[maxStart+i] = sMax[i]
		}
	}
	return strings.Repeat(" ", labelW) + string(line) + "\n"
}

func autoRangeY(allY []float64) (float64, float64, bool) {
	if len(allY) == 0 {
		return 0, 0, false
	}
	sorted := append([]float64(nil), allY...)
	sort.Float64s(sorted)

	yMin := sorted[0]
	yMax := sorted[len(sorted)-1]

	if len(sorted) >= 20 {
		p05 := quantileSorted(sorted, 0.05)
		p95 := quantileSorted(sorted, 0.95)
		coreSpan := p95 - p05
		fullSpan := yMax - yMin
		if coreSpan > 1e-9 && fullSpan > 15*coreSpan {
			yMin = p05 - 0.6*coreSpan
			yMax = p95 + 0.6*coreSpan
		}
	}

	if math.Abs(yMax-yMin) < 1e-12 {
		pad := math.Max(1.0, 0.1*math.Abs(yMax))
		yMin -= pad
		yMax += pad
	} else {
		if yMin < 0 && yMax > 0 && math.Abs(yMax+yMin) <= 0.02*(yMax-yMin) {
			m := math.Max(-yMin, yMax)
			yMin, yMax = -m, m
		}
		pad := 0.06 * (yMax - yMin)
		yMin -= pad
		yMax += pad
	}
	return yMin, yMax, true
}

func drawZeroAxes(bc *brailleCanvas, xMin, xMax, yMin, yMax float64) {
	if yMin <= 0 && yMax >= 0 {
		py := int(math.Round((yMax - 0) / (yMax - yMin) * float64(bc.pixH-1)))
		for px := 0; px < bc.pixW; px += 2 {
			bc.setPixel(px, py, -2)
		}
	}
	if xMin <= 0 && xMax >= 0 {
		px := int(math.Round((0 - xMin) / (xMax - xMin) * float64(bc.pixW-1)))
		for py := 0; py < bc.pixH; py += 2 {
			bc.setPixel(px, py, -2)
		}
	}
}

func (e *Engine) PlotFunctions(exprs []string, xMin, xMax, yMinOpt, yMaxOpt float64, charW, charH int, useColor bool) (string, error) {
	if len(exprs) == 0 {
		return "", fmt.Errorf("at least one expression required to plot")
	}
	if charW <= 0 {
		charW = 64
	}
	if charH <= 0 {
		charH = 18
	}
	if xMax <= xMin {
		xMin, xMax = -10, 10
	}

	nodes := make([]ASTNode, len(exprs))
	for i, ex := range exprs {
		nd, err := parseMath(strings.TrimSpace(ex))
		if err != nil {
			return "", fmt.Errorf("error in expression %q: %v", ex, err)
		}
		if bin, ok := nd.(BinaryNode); ok && (bin.Op == "=" || bin.Op == "==") {
			if id, isY := bin.Left.(IdentNode); isY && (id.Name == "y" || id.Name == "Y") {
				nd = bin.Right
			} else {
				nd = BinaryNode{Op: "-", Left: bin.Left, Right: bin.Right}
			}
		}
		nodes[i] = e.expandSymbolicRefs(nd)
	}

	bc := newBrailleCanvas(charW, charH)
	samples := bc.pixW * 2
	xs := make([]float64, samples)
	stepX := (xMax - xMin) / float64(samples-1)
	for i := 0; i < samples; i++ {
		xs[i] = xMin + float64(i)*stepX
	}

	evalAt := func(nd ASTNode, x float64) (float64, bool) {
		scope := map[string]Value{
			"x": ScalarReal(x),
			"X": ScalarReal(x),
		}
		v, err := e.evalAST(nd, scope, 0)
		if err != nil || !v.IsReal() {
			return 0, false
		}
		r := real(v.S)
		if math.IsNaN(r) || math.IsInf(r, 0) {
			return 0, false
		}
		return r, true
	}

	ysAll := make([][]float64, len(nodes))
	validAll := make([][]bool, len(nodes))
	var finitePool []float64

	for sIdx, nd := range nodes {
		ysAll[sIdx] = make([]float64, samples)
		validAll[sIdx] = make([]bool, samples)
		for i, x := range xs {
			if y, ok := evalAt(nd, x); ok {
				ysAll[sIdx][i] = y
				validAll[sIdx][i] = true
				finitePool = append(finitePool, y)
			}
		}
	}

	if len(finitePool) == 0 {
		return "Unable to plot: no valid real points in domain.\n", nil
	}

	var yMin, yMax float64
	if !math.IsNaN(yMinOpt) && !math.IsNaN(yMaxOpt) && yMaxOpt > yMinOpt {
		yMin, yMax = yMinOpt, yMaxOpt
	} else {
		var ok bool
		yMin, yMax, ok = autoRangeY(finitePool)
		if !ok {
			return "Unable to plot: no valid real points in domain.\n", nil
		}
	}

	drawZeroAxes(bc, xMin, xMax, yMin, yMax)
	ySpan := yMax - yMin

	mapX := func(x float64) int {
		return int(math.Round((x - xMin) / (xMax - xMin) * float64(bc.pixW-1)))
	}
	mapY := func(y float64) int {
		return int(math.Round((yMax - y) / ySpan * float64(bc.pixH-1)))
	}

	for sIdx, nd := range nodes {
		ys := ysAll[sIdx]
		val := validAll[sIdx]
		for i := 0; i < samples; i++ {
			if !val[i] {
				continue
			}
			px1, py1 := mapX(xs[i]), mapY(ys[i])
			if i > 0 && val[i-1] {
				dy := math.Abs(ys[i] - ys[i-1])
				isDiscontinuous := false
				if dy > 0.65*ySpan {
					midY, okMid := evalAt(nd, 0.5*(xs[i-1]+xs[i]))
					if !okMid {
						isDiscontinuous = true
					} else {
						lo, hi := math.Min(ys[i-1], ys[i]), math.Max(ys[i-1], ys[i])
						if midY < lo || midY > hi {
							isDiscontinuous = true
						}
					}
				}
				if !isDiscontinuous {
					px0, py0 := mapX(xs[i-1]), mapY(ys[i-1])
					if (py0 >= 0 && py0 < bc.pixH) || (py1 >= 0 && py1 < bc.pixH) {
						bc.drawLine(px0, py0, px1, py1, sIdx)
						continue
					}
				}
			}
			bc.setPixel(px1, py1, sIdx)
		}
	}

	var header strings.Builder
	if len(exprs) == 1 {
		header.WriteString(fmt.Sprintf("f(x) = %s   [x: %.4g to %.4g, y: %.4g to %.4g]\n",
			strings.TrimSpace(exprs[0]), xMin, xMax, yMin, yMax))
	} else {
		header.WriteString("Legend: ")
		for i, ex := range exprs {
			tag := fmt.Sprintf("f%d(x) = %s", i+1, strings.TrimSpace(ex))
			if useColor {
				tag = seriesColors[i%len(seriesColors)] + tag + ansiReset
			}
			if i > 0 {
				header.WriteString("  |  ")
			}
			header.WriteString(tag)
		}
		header.WriteString(fmt.Sprintf("\nWindow: [x: %.4g to %.4g, y: %.4g to %.4g]\n", xMin, xMax, yMin, yMax))
	}

	header.WriteString(bc.renderWithMargin(yMax, yMin, useColor))
	header.WriteString(formatXAxisFooter(xMin, xMax, charW))

	if len(nodes) == 1 {
		roots := solveTranscendentalReal(func(x float64) float64 {
			y, ok := evalAt(nodes[0], x)
			if !ok {
				return math.NaN()
			}
			return y
		}, xMin, xMax)
		if len(roots) > 0 && len(roots) <= 8 {
			rStrs := make([]string, len(roots))
			for i, z := range roots {
				rStrs[i] = formatFloatClean(real(z))
			}
			header.WriteString("Zeros in window: x = " + strings.Join(rStrs, ", ") + "\n")
		}
	}

	return header.String(), nil
}

func (e *Engine) PlotParametric(xExpr, yExpr string, tMin, tMax float64, charW, charH int, useColor bool) (string, error) {
	if charW <= 0 {
		charW = 64
	}
	if charH <= 0 {
		charH = 18
	}
	if tMax <= tMin {
		tMin, tMax = 0, 2*math.Pi
	}

	xNode, err := parseMath(strings.TrimSpace(xExpr))
	if err != nil {
		return "", err
	}
	yNode, err := parseMath(strings.TrimSpace(yExpr))
	if err != nil {
		return "", err
	}
	xNode = e.expandSymbolicRefs(xNode)
	yNode = e.expandSymbolicRefs(yNode)

	const samples = 600
	xs := make([]float64, samples)
	ys := make([]float64, samples)
	valid := make([]bool, samples)
	var xPool, yPool []float64

	for i := 0; i < samples; i++ {
		t := tMin + (tMax-tMin)*float64(i)/float64(samples-1)
		scope := map[string]Value{"t": ScalarReal(t), "T": ScalarReal(t)}
		xv, err1 := e.evalAST(xNode, scope, 0)
		yv, err2 := e.evalAST(yNode, scope, 0)
		if err1 == nil && err2 == nil && xv.IsReal() && yv.IsReal() {
			xr, yr := real(xv.S), real(yv.S)
			if !math.IsNaN(xr) && !math.IsInf(xr, 0) && !math.IsNaN(yr) && !math.IsInf(yr, 0) {
				xs[i], ys[i] = xr, yr
				valid[i] = true
				xPool = append(xPool, xr)
				yPool = append(yPool, yr)
			}
		}
	}

	xMin, xMax, ok1 := autoRangeY(xPool)
	yMin, yMax, ok2 := autoRangeY(yPool)
	if !ok1 || !ok2 {
		return "Unable to plot parametric curve: no valid points.\n", nil
	}

	bc := newBrailleCanvas(charW, charH)
	drawZeroAxes(bc, xMin, xMax, yMin, yMax)

	mapX := func(x float64) int {
		return int(math.Round((x - xMin) / (xMax - xMin) * float64(bc.pixW-1)))
	}
	mapY := func(y float64) int {
		return int(math.Round((yMax - y) / (yMax - yMin) * float64(bc.pixH-1)))
	}

	for i := 0; i < samples; i++ {
		if !valid[i] {
			continue
		}
		px1, py1 := mapX(xs[i]), mapY(ys[i])
		if i > 0 && valid[i-1] {
			px0, py0 := mapX(xs[i-1]), mapY(ys[i-1])
			if absInt(px1-px0) < bc.pixW/2 && absInt(py1-py0) < bc.pixH/2 {
				bc.drawLine(px0, py0, px1, py1, 0)
				continue
			}
		}
		bc.setPixel(px1, py1, 0)
	}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("Parametric: x(t) = %s, y(t) = %s   [t: %.4g to %.4g]\n",
		strings.TrimSpace(xExpr), strings.TrimSpace(yExpr), tMin, tMax))
	sb.WriteString(bc.renderWithMargin(yMax, yMin, useColor))
	sb.WriteString(formatXAxisFooter(xMin, xMax, charW))
	return sb.String(), nil
}

func (e *Engine) PlotPolar(rExpr string, tMin, tMax float64, charW, charH int, useColor bool) (string, error) {
	if tMax <= tMin {
		tMin, tMax = 0, 2*math.Pi
	}
	rNode, err := parseMath(strings.TrimSpace(rExpr))
	if err != nil {
		return "", err
	}
	rNode = e.expandSymbolicRefs(rNode)

	if charW <= 0 {
		charW = 64
	}
	if charH <= 0 {
		charH = 18
	}

	const samples = 720
	xs := make([]float64, samples)
	ys := make([]float64, samples)
	valid := make([]bool, samples)
	var xPool, yPool []float64

	for i := 0; i < samples; i++ {
		theta := tMin + (tMax-tMin)*float64(i)/float64(samples-1)
		scope := map[string]Value{
			"t":     ScalarReal(theta),
			"theta": ScalarReal(theta),
			"x":     ScalarReal(theta),
		}
		rv, err := e.evalAST(rNode, scope, 0)
		if err == nil && rv.IsReal() {
			r := real(rv.S)
			if !math.IsNaN(r) && !math.IsInf(r, 0) {
				xr := r * math.Cos(theta)
				yr := r * math.Sin(theta)
				xs[i], ys[i] = xr, yr
				valid[i] = true
				xPool = append(xPool, xr)
				yPool = append(yPool, yr)
			}
		}
	}

	xMin, xMax, ok1 := autoRangeY(xPool)
	yMin, yMax, ok2 := autoRangeY(yPool)
	if !ok1 || !ok2 {
		return "Unable to plot polar curve: no valid points.\n", nil
	}

	bc := newBrailleCanvas(charW, charH)
	drawZeroAxes(bc, xMin, xMax, yMin, yMax)

	mapX := func(x float64) int {
		return int(math.Round((x - xMin) / (xMax - xMin) * float64(bc.pixW-1)))
	}
	mapY := func(y float64) int {
		return int(math.Round((yMax - y) / (yMax - yMin) * float64(bc.pixH-1)))
	}

	for i := 0; i < samples; i++ {
		if !valid[i] {
			continue
		}
		px1, py1 := mapX(xs[i]), mapY(ys[i])
		if i > 0 && valid[i-1] {
			px0, py0 := mapX(xs[i-1]), mapY(ys[i-1])
			if absInt(px1-px0) < bc.pixW/2 && absInt(py1-py0) < bc.pixH/2 {
				bc.drawLine(px0, py0, px1, py1, 1)
				continue
			}
		}
		bc.setPixel(px1, py1, 1)
	}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("Polar: r(t) = %s   [t: %.4g to %.4g]\n", strings.TrimSpace(rExpr), tMin, tMax))
	sb.WriteString(bc.renderWithMargin(yMax, yMin, useColor))
	sb.WriteString(formatXAxisFooter(xMin, xMax, charW))
	return sb.String(), nil
}

func (e *Engine) PlotSlopeField(dydxExpr string, xMin, xMax, yMin, yMax float64, charW, charH int, useColor bool) (string, error) {
	if charW <= 0 {
		charW = 64
	}
	if charH <= 0 {
		charH = 18
	}
	if xMax <= xMin {
		xMin, xMax = -5, 5
	}
	if yMax <= yMin {
		yMin, yMax = -5, 5
	}

	nd, err := parseMath(strings.TrimSpace(dydxExpr))
	if err != nil {
		return "", err
	}
	nd = e.expandSymbolicRefs(nd)

	bc := newBrailleCanvas(charW, charH)
	drawZeroAxes(bc, xMin, xMax, yMin, yMax)

	gridCols := 19
	gridRows := 13
	segHalfPix := 3.0

	for r := 0; r < gridRows; r++ {
		y := yMin + (yMax-yMin)*(float64(r)+0.5)/float64(gridRows)
		for c := 0; c < gridCols; c++ {
			x := xMin + (xMax-xMin)*(float64(c)+0.5)/float64(gridCols)
			scope := map[string]Value{"x": ScalarReal(x), "y": ScalarReal(y)}
			v, err := e.evalAST(nd, scope, 0)
			if err != nil || !v.IsReal() {
				continue
			}
			slope := real(v.S)
			if math.IsNaN(slope) || math.IsInf(slope, 0) {
				continue
			}
			cx := (x - xMin) / (xMax - xMin) * float64(bc.pixW-1)
			cy := (yMax - y) / (yMax - yMin) * float64(bc.pixH-1)

			pixSlope := -slope * (float64(bc.pixH-1) / (yMax - yMin)) / (float64(bc.pixW-1) / (xMax - xMin))
			norm := math.Hypot(1.0, pixSlope)
			dx := segHalfPix / norm
			dy := segHalfPix * pixSlope / norm

			bc.drawLine(int(math.Round(cx-dx)), int(math.Round(cy-dy)), int(math.Round(cx+dx)), int(math.Round(cy+dy)), 0)
		}
	}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("Slope Field: dy/dx = %s   [x: %.4g to %.4g, y: %.4g to %.4g]\n",
		strings.TrimSpace(dydxExpr), xMin, xMax, yMin, yMax))
	sb.WriteString(bc.renderWithMargin(yMax, yMin, useColor))
	sb.WriteString(formatXAxisFooter(xMin, xMax, charW))
	return sb.String(), nil
}

func (e *Engine) RenderTable(exprs []string, xMin, xMax, step float64) (string, error) {
	if len(exprs) == 0 {
		return "", fmt.Errorf("table requires expression")
	}
	if step <= 0 {
		step = (xMax - xMin) / 10.0
		if step <= 0 {
			step = 1.0
		}
	}
	nodes := make([]ASTNode, len(exprs))
	for i, ex := range exprs {
		nd, err := parseMath(strings.TrimSpace(ex))
		if err != nil {
			return "", err
		}
		nodes[i] = e.expandSymbolicRefs(nd)
	}

	var rows [][]string
	head := []string{"x"}
	for _, ex := range exprs {
		head = append(head, strings.TrimSpace(ex))
	}
	rows = append(rows, head)

	count := 0
	for x := xMin; x <= xMax+1e-12*math.Max(1.0, math.Abs(step)) && count < 200; x += step {
		xv := x
		if math.Abs(xv) < 1e-13 {
			xv = 0
		}
		row := []string{formatFloatClean(xv)}
		scope := map[string]Value{"x": ScalarReal(xv), "X": ScalarReal(xv)}
		for _, nd := range nodes {
			val, err := e.evalAST(nd, scope, 0)
			if err != nil {
				row = append(row, "undef")
			} else if val.Kind == KindScalar {
				row = append(row, formatComplex(val.S, false))
			} else {
				row = append(row, val.String())
			}
		}
		rows = append(rows, row)
		count++
	}

	widths := make([]int, len(head))
	for _, r := range rows {
		for c, cell := range r {
			if len(cell) > widths[c] {
				widths[c] = len(cell)
			}
		}
	}

	var sb strings.Builder
	for idx, r := range rows {
		for c, cell := range r {
			pad := widths[c] - len(cell)
			if c > 0 {
				sb.WriteString("  |  ")
			}
			sb.WriteString(strings.Repeat(" ", pad) + cell)
		}
		sb.WriteRune('\n')
		if idx == 0 {
			totalW := 0
			for c, w := range widths {
				if c > 0 {
					totalW += 5
				}
				totalW += w
			}
			sb.WriteString(strings.Repeat("-", totalW) + "\n")
		}
	}
	return sb.String(), nil
}

func RenderGraph(expr string, xMin, xMax float64, width, height int) string {
	out, err := defaultEngine.PlotFunctions([]string{expr}, xMin, xMax, math.NaN(), math.NaN(), width, height, false)
	if err != nil {
		return fmt.Sprintf("Error plotting %s: %v\n", expr, err)
	}
	return out
}
