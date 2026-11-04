package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

func cloneMat(a [][]complex128) [][]complex128 {
	out := make([][]complex128, len(a))
	for i := range a {
		out[i] = make([]complex128, len(a[i]))
		copy(out[i], a[i])
	}
	return out
}

func matIdentity(n int) [][]complex128 {
	m := make([][]complex128, n)
	for i := 0; i < n; i++ {
		m[i] = make([]complex128, n)
		m[i][i] = 1
	}
	return m
}

func matZeros(r, c int) [][]complex128 {
	m := make([][]complex128, r)
	for i := 0; i < r; i++ {
		m[i] = make([]complex128, c)
	}
	return m
}

func matAdd(a, b [][]complex128, sign complex128) ([][]complex128, error) {
	if len(a) == 0 || len(b) == 0 || len(a) != len(b) || len(a[0]) != len(b[0]) {
		return nil, fmt.Errorf("matrix dimension mismatch")
	}
	r, c := len(a), len(a[0])
	out := matZeros(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[i][j] = cleanComplex(a[i][j] + sign*b[i][j])
		}
	}
	return out, nil
}

func matScale(a [][]complex128, s complex128) [][]complex128 {
	r := len(a)
	if r == 0 {
		return nil
	}
	c := len(a[0])
	out := matZeros(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[i][j] = cleanComplex(a[i][j] * s)
		}
	}
	return out
}

func matMul(a, b [][]complex128) ([][]complex128, error) {
	if len(a) == 0 || len(b) == 0 || len(a[0]) != len(b) {
		return nil, fmt.Errorf("matrix multiplication dimension mismatch")
	}
	r, kDim, c := len(a), len(a[0]), len(b[0])
	out := matZeros(r, c)
	for i := 0; i < r; i++ {
		for k := 0; k < kDim; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			for j := 0; j < c; j++ {
				out[i][j] += aik * b[k][j]
			}
		}
		for j := 0; j < c; j++ {
			out[i][j] = cleanComplex(out[i][j])
		}
	}
	return out, nil
}

func matVecMul(a [][]complex128, v []complex128) ([]complex128, error) {
	if len(a) == 0 || len(a[0]) != len(v) {
		return nil, fmt.Errorf("matrix-vector dimension mismatch (%dx%d vs %d)", len(a), len(a[0]), len(v))
	}
	out := make([]complex128, len(a))
	for i := range a {
		var sum complex128
		for j, z := range v {
			sum += a[i][j] * z
		}
		out[i] = cleanComplex(sum)
	}
	return out, nil
}

func vecMatMul(v []complex128, a [][]complex128) ([]complex128, error) {
	if len(a) == 0 || len(v) != len(a) {
		return nil, fmt.Errorf("vector-matrix dimension mismatch (%d vs %dx%d)", len(v), len(a), len(a[0]))
	}
	cols := len(a[0])
	out := make([]complex128, cols)
	for j := 0; j < cols; j++ {
		var sum complex128
		for i := range v {
			sum += v[i] * a[i][j]
		}
		out[j] = cleanComplex(sum)
	}
	return out, nil
}

func matTranspose(a [][]complex128) [][]complex128 {
	if len(a) == 0 {
		return nil
	}
	r, c := len(a), len(a[0])
	out := matZeros(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[j][i] = a[i][j]
		}
	}
	return out
}

func matConjTranspose(a [][]complex128) [][]complex128 {
	if len(a) == 0 {
		return nil
	}
	r, c := len(a), len(a[0])
	out := matZeros(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[j][i] = cmplx.Conj(a[i][j])
		}
	}
	return out
}

func matPow(a [][]complex128, exp int) ([][]complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return nil, fmt.Errorf("matrix power requires square matrix")
	}
	n := len(a)
	if exp == 0 {
		return matIdentity(n), nil
	}
	base := cloneMat(a)
	var err error
	if exp < 0 {
		base, err = matInverse(base)
		if err != nil {
			return nil, err
		}
		exp = -exp
	}
	res := matIdentity(n)
	for exp > 0 {
		if exp%2 == 1 {
			res, err = matMul(res, base)
			if err != nil {
				return nil, err
			}
		}
		exp /= 2
		if exp > 0 {
			base, err = matMul(base, base)
			if err != nil {
				return nil, err
			}
		}
	}
	return res, nil
}

func matNormScale(a [][]complex128) float64 {
	mx := 0.0
	for i := range a {
		for j := range a[i] {
			if v := cmplx.Abs(a[i][j]); v > mx {
				mx = v
			}
		}
	}
	if mx == 0 {
		return 1.0
	}
	return mx
}

func matDet(a [][]complex128) (complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return 0, fmt.Errorf("det requires square matrix")
	}
	n := len(a)
	m := cloneMat(a)
	tol := 1e-13 * matNormScale(m)
	det := complex(1, 0)

	for col := 0; col < n; col++ {
		piv := col
		maxVal := cmplx.Abs(m[col][col])
		for row := col + 1; row < n; row++ {
			if v := cmplx.Abs(m[row][col]); v > maxVal {
				maxVal = v
				piv = row
			}
		}
		if maxVal <= tol {
			return 0, nil
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			det = -det
		}
		pivVal := m[col][col]
		det *= pivVal
		for row := col + 1; row < n; row++ {
			factor := m[row][col] / pivVal
			if factor == 0 {
				continue
			}
			for j := col + 1; j < n; j++ {
				m[row][j] -= factor * m[col][j]
			}
		}
	}
	return cleanComplex(det), nil
}

func matInverse(a [][]complex128) ([][]complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return nil, fmt.Errorf("inv requires square matrix")
	}
	n := len(a)
	aug := matZeros(n, 2*n)
	for i := 0; i < n; i++ {
		copy(aug[i][:n], a[i])
		aug[i][n+i] = 1
	}
	tol := 1e-13 * matNormScale(a)

	for col := 0; col < n; col++ {
		piv := col
		maxVal := cmplx.Abs(aug[col][col])
		for row := col + 1; row < n; row++ {
			if v := cmplx.Abs(aug[row][col]); v > maxVal {
				maxVal = v
				piv = row
			}
		}
		if maxVal <= tol {
			return nil, fmt.Errorf("matrix is singular (non-invertible)")
		}
		if piv != col {
			aug[col], aug[piv] = aug[piv], aug[col]
		}
		pivVal := aug[col][col]
		for j := col; j < 2*n; j++ {
			aug[col][j] /= pivVal
		}
		for row := 0; row < n; row++ {
			if row == col {
				continue
			}
			factor := aug[row][col]
			if factor == 0 {
				continue
			}
			for j := col; j < 2*n; j++ {
				aug[row][j] -= factor * aug[col][j]
			}
		}
	}

	inv := matZeros(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			z := cleanComplex(aug[i][n+j])
			if cmplx.Abs(z) <= 1e-14*matNormScale(a) {
				z = 0
			}
			inv[i][j] = z
		}
	}
	return inv, nil
}

func matRREF(a [][]complex128, fullReduce bool) ([][]complex128, []int) {
	if len(a) == 0 {
		return nil, nil
	}
	rows, cols := len(a), len(a[0])
	m := cloneMat(a)
	scale := matNormScale(a)
	tol := 1e-12 * scale

	var pivots []int
	leadRow := 0
	for col := 0; col < cols && leadRow < rows; col++ {
		piv := leadRow
		maxVal := cmplx.Abs(m[leadRow][col])
		for r := leadRow + 1; r < rows; r++ {
			if v := cmplx.Abs(m[r][col]); v > maxVal {
				maxVal = v
				piv = r
			}
		}
		if maxVal <= tol {
			m[leadRow][col] = 0
			continue
		}
		if piv != leadRow {
			m[leadRow], m[piv] = m[piv], m[leadRow]
		}
		pivVal := m[leadRow][col]
		for j := col; j < cols; j++ {
			m[leadRow][j] /= pivVal
		}
		startRow := leadRow + 1
		if fullReduce {
			startRow = 0
		}
		for r := startRow; r < rows; r++ {
			if r == leadRow {
				continue
			}
			factor := m[r][col]
			if cmplx.Abs(factor) <= tol {
				m[r][col] = 0
				continue
			}
			for j := col; j < cols; j++ {
				m[r][j] -= factor * m[leadRow][j]
			}
		}
		pivots = append(pivots, col)
		leadRow++
	}

	for r := 0; r < rows; r++ {
		for c := 0; c < cols; c++ {
			if cmplx.Abs(m[r][c]) <= tol {
				m[r][c] = 0
			} else {
				m[r][c] = cleanComplex(m[r][c])
			}
		}
	}
	return m, pivots
}

func matNullspace(a [][]complex128) [][]complex128 {
	if len(a) == 0 {
		return nil
	}
	cols := len(a[0])
	rrefM, pivots := matRREF(a, true)
	isPivot := make(map[int]int, len(pivots))
	for rowIdx, colIdx := range pivots {
		isPivot[colIdx] = rowIdx
	}

	var basisCols [][]complex128
	for freeCol := 0; freeCol < cols; freeCol++ {
		if _, ok := isPivot[freeCol]; ok {
			continue
		}
		vec := make([]complex128, cols)
		vec[freeCol] = 1
		for pivCol, pivRow := range isPivot {
			vec[pivCol] = cleanComplex(-rrefM[pivRow][freeCol])
		}
		basisCols = append(basisCols, vec)
	}
	if len(basisCols) == 0 {
		return matZeros(cols, 1)
	}
	out := matZeros(cols, len(basisCols))
	for j, bVec := range basisCols {
		for i := 0; i < cols; i++ {
			out[i][j] = bVec[i]
		}
	}
	return out
}

func matLinearSolve(a [][]complex128, b [][]complex128) ([][]complex128, error) {
	m := len(a)
	if m == 0 || len(b) != m {
		return nil, fmt.Errorf("linsolve dimension mismatch")
	}
	n := len(a[0])
	rhsCols := len(b[0])

	if m == n {
		aug := matZeros(n, n+rhsCols)
		for i := 0; i < n; i++ {
			copy(aug[i][:n], a[i])
			copy(aug[i][n:], b[i])
		}
		rrefAug, pivots := matRREF(aug, true)
		if len(pivots) >= n && pivots[n-1] == n-1 {
			sol := matZeros(n, rhsCols)
			for i := 0; i < n; i++ {
				for j := 0; j < rhsCols; j++ {
					sol[i][j] = cleanComplex(rrefAug[i][n+j])
				}
			}
			return sol, nil
		}
	}

	ah := matConjTranspose(a)
	ata, err := matMul(ah, a)
	if err != nil {
		return nil, err
	}
	atb, err := matMul(ah, b)
	if err != nil {
		return nil, err
	}
	invAta, err := matInverse(ata)
	if err != nil {
		return nil, fmt.Errorf("system is singular or underdetermined")
	}
	return matMul(invAta, atb)
}

func hessenbergReduce(a [][]complex128) [][]complex128 {
	n := len(a)
	h := cloneMat(a)
	for k := 0; k < n-2; k++ {
		mLen := n - (k + 1)
		x := make([]complex128, mLen)
		var normSq float64
		for i := 0; i < mLen; i++ {
			x[i] = h[k+1+i][k]
			normSq += real(x[i])*real(x[i]) + imag(x[i])*imag(x[i])
		}
		normX := math.Sqrt(normSq)
		if normX <= 1e-15 {
			continue
		}
		var phase complex128 = 1
		if cmplx.Abs(x[0]) > 1e-30 {
			phase = x[0] / complex(cmplx.Abs(x[0]), 0)
		}
		x[0] += phase * complex(normX, 0)
		var vNormSq float64
		for _, xi := range x {
			vNormSq += real(xi)*real(xi) + imag(xi)*imag(xi)
		}
		if vNormSq <= 1e-30 {
			continue
		}

		for j := k; j < n; j++ {
			var dot complex128
			for i := 0; i < mLen; i++ {
				dot += cmplx.Conj(x[i]) * h[k+1+i][j]
			}
			scale := 2 * dot / complex(vNormSq, 0)
			for i := 0; i < mLen; i++ {
				h[k+1+i][j] -= scale * x[i]
			}
		}

		for i := 0; i < n; i++ {
			var dot complex128
			for j := 0; j < mLen; j++ {
				dot += h[i][k+1+j] * x[j]
			}
			scale := 2 * dot / complex(vNormSq, 0)
			for j := 0; j < mLen; j++ {
				h[i][k+1+j] -= scale * cmplx.Conj(x[j])
			}
		}
	}
	return h
}

func matEigenvalues(a [][]complex128) ([]complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return nil, fmt.Errorf("eigenvalues require square matrix")
	}
	n := len(a)
	if n == 1 {
		return []complex128{cleanComplex(a[0][0])}, nil
	}
	if n == 2 {
		tr := a[0][0] + a[1][1]
		det := a[0][0]*a[1][1] - a[0][1]*a[1][0]
		disc := cmplx.Sqrt(tr*tr - 4*det)
		return sortComplexSlice([]complex128{(tr - disc) / 2, (tr + disc) / 2}), nil
	}

	h := hessenbergReduce(a)
	scale := matNormScale(a)
	var eigs []complex128
	m := n - 1
	iters := 0
	maxTotalIters := 400 * n

	for m >= 0 && iters < maxTotalIters {
		if m == 0 {
			eigs = append(eigs, cleanComplex(h[0][0]))
			m--
			continue
		}

		subMag := cmplx.Abs(h[m][m-1])
		diagMag := cmplx.Abs(h[m-1][m-1]) + cmplx.Abs(h[m][m])
		if diagMag == 0 {
			diagMag = scale
		}
		if subMag <= 1e-13*diagMag {
			eigs = append(eigs, cleanComplex(h[m][m]))
			h[m][m-1] = 0
			m--
			iters = 0
			continue
		}

		if m == 1 {
			tr := h[0][0] + h[1][1]
			det := h[0][0]*h[1][1] - h[0][1]*h[1][0]
			disc := cmplx.Sqrt(tr*tr - 4*det)
			eigs = append(eigs, cleanComplex((tr-disc)/2), cleanComplex((tr+disc)/2))
			m -= 2
			continue
		}

		iters++
		var shift complex128
		if iters%15 == 0 {
			shift = h[m][m] + complex(0.75*subMag, -0.4375*subMag)
		} else {
			a11, a12 := h[m-1][m-1], h[m-1][m]
			a21, a22 := h[m][m-1], h[m][m]
			tr := a11 + a22
			det := a11*a22 - a12*a21
			disc := cmplx.Sqrt(tr*tr - 4*det)
			l1, l2 := (tr+disc)/2, (tr-disc)/2
			if cmplx.Abs(l1-a22) < cmplx.Abs(l2-a22) {
				shift = l1
			} else {
				shift = l2
			}
		}

		for i := 0; i <= m; i++ {
			h[i][i] -= shift
		}

		cs := make([]complex128, m)
		sn := make([]complex128, m)
		for k := 0; k < m; k++ {
			x, y := h[k][k], h[k+1][k]
			r := math.Hypot(cmplx.Abs(x), cmplx.Abs(y))
			if r == 0 {
				cs[k] = 1
				sn[k] = 0
				continue
			}
			c := complex(cmplx.Abs(x)/r, 0)
			var phaseX complex128 = 1
			if cmplx.Abs(x) > 0 {
				phaseX = x / complex(cmplx.Abs(x), 0)
			}
			s := phaseX * cmplx.Conj(y) / complex(r, 0)
			cs[k] = c
			sn[k] = s

			for j := k; j <= m; j++ {
				t1 := cmplx.Conj(c)*h[k][j] + s*h[k+1][j]
				t2 := -cmplx.Conj(s)*h[k][j] + c*h[k+1][j]
				h[k][j] = t1
				h[k+1][j] = t2
			}
		}

		for k := 0; k < m; k++ {
			c, s := cs[k], sn[k]
			for i := 0; i <= k+1; i++ {
				t1 := h[i][k]*c + h[i][k+1]*cmplx.Conj(s)
				t2 := -h[i][k]*s + h[i][k+1]*cmplx.Conj(c)
				h[i][k] = t1
				h[i][k+1] = t2
			}
		}

		for i := 0; i <= m; i++ {
			h[i][i] += shift
		}
	}

	for i := 0; i <= m; i++ {
		eigs = append(eigs, cleanComplex(h[i][i]))
	}

	out := sortComplexSlice(eigs)
	for i := range out {
		if cmplx.Abs(out[i]) <= 1e-13*scale {
			out[i] = 0
		}
	}
	return out, nil
}

func matEigenvectors(a [][]complex128) ([][]complex128, error) {
	eigs, err := matEigenvalues(a)
	if err != nil {
		return nil, err
	}
	n := len(a)
	var cols [][]complex128
	usedEigs := make([]complex128, 0, n)

	for _, lam := range eigs {
		already := false
		for _, prev := range usedEigs {
			if cmplx.Abs(prev-lam) <= 1e-8*math.Max(1.0, cmplx.Abs(lam)) {
				already = true
				break
			}
		}
		if already {
			continue
		}
		usedEigs = append(usedEigs, lam)

		shifted := cloneMat(a)
		for i := 0; i < n; i++ {
			shifted[i][i] -= lam
		}
		ns := matNullspace(shifted)
		hasNonZero := false
		for c := 0; c < len(ns[0]); c++ {
			vec := make([]complex128, n)
			var normSq float64
			for r := 0; r < n; r++ {
				vec[r] = ns[r][c]
				normSq += real(vec[r])*real(vec[r]) + imag(vec[r])*imag(vec[r])
			}
			if normSq > 1e-20 {
				hasNonZero = true
				nm := math.Sqrt(normSq)
				for r := 0; r < n; r++ {
					vec[r] = cleanComplex(vec[r] / complex(nm, 0))
				}
				cols = append(cols, vec)
			}
		}
		if !hasNonZero {
			pert := cloneMat(shifted)
			for i := 0; i < n; i++ {
				pert[i][i] -= complex(1e-10*(1+cmplx.Abs(lam)), 0)
			}
			rhs := matZeros(n, 1)
			for i := 0; i < n; i++ {
				rhs[i][0] = complex(1.0/float64(i+1), 0)
			}
			sol, err := matLinearSolve(pert, rhs)
			if err == nil {
				vec := make([]complex128, n)
				var normSq float64
				for r := 0; r < n; r++ {
					vec[r] = sol[r][0]
					normSq += real(vec[r])*real(vec[r]) + imag(vec[r])*imag(vec[r])
				}
				nm := math.Sqrt(normSq)
				if nm > 0 {
					for r := 0; r < n; r++ {
						vec[r] = cleanComplex(vec[r] / complex(nm, 0))
					}
				}
				cols = append(cols, vec)
			}
		}
	}

	if len(cols) == 0 {
		return matIdentity(n), nil
	}
	out := matZeros(n, len(cols))
	for j, v := range cols {
		for i := 0; i < n; i++ {
			out[i][j] = v[i]
		}
	}
	return out, nil
}

func faddeevLeVerrier(a [][]complex128) ([]complex128, [][]complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return nil, nil, fmt.Errorf("square matrix required")
	}
	n := len(a)
	coeffs := make([]complex128, n+1)
	coeffs[n] = 1
	m := matZeros(n, n)
	ident := matIdentity(n)
	var adj [][]complex128

	for k := 1; k <= n; k++ {
		if k == 1 {
			m = cloneMat(ident)
		} else {
			am, err := matMul(a, m)
			if err != nil {
				return nil, nil, err
			}
			prevC := coeffs[n-k+1]
			m, _ = matAdd(am, matScale(ident, prevC), 1)
		}
		if k == n {
			sign := complex(1, 0)
			if (n-1)%2 == 1 {
				sign = complex(-1, 0)
			}
			adj = matScale(m, sign)
		}
		am, err := matMul(a, m)
		if err != nil {
			return nil, nil, err
		}
		var tr complex128
		for i := 0; i < n; i++ {
			tr += am[i][i]
		}
		coeffs[n-k] = cleanComplex(-tr / complex(float64(k), 0))
	}
	if (n-1)%2 == 1 {
		adj = matScale(m, -1)
	} else {
		adj = cloneMat(m)
	}
	return coeffs, adj, nil
}

func matGramSchmidt(a [][]complex128) [][]complex128 {
	if len(a) == 0 {
		return nil
	}
	rows, cols := len(a), len(a[0])
	var ortho [][]complex128

	for j := 0; j < cols; j++ {
		v := make([]complex128, rows)
		for i := 0; i < rows; i++ {
			v[i] = a[i][j]
		}
		for _, u := range ortho {
			var dot complex128
			for i := 0; i < rows; i++ {
				dot += cmplx.Conj(u[i]) * v[i]
			}
			for i := 0; i < rows; i++ {
				v[i] -= dot * u[i]
			}
		}
		var normSq float64
		for i := 0; i < rows; i++ {
			normSq += real(v[i])*real(v[i]) + imag(v[i])*imag(v[i])
		}
		nm := math.Sqrt(normSq)
		if nm > 1e-12 {
			for i := 0; i < rows; i++ {
				v[i] = cleanComplex(v[i] / complex(nm, 0))
			}
			ortho = append(ortho, v)
		}
	}

	if len(ortho) == 0 {
		return matZeros(rows, 1)
	}
	out := matZeros(rows, len(ortho))
	for j, u := range ortho {
		for i := 0; i < rows; i++ {
			out[i][j] = u[i]
		}
	}
	return out
}

func matExpm(a [][]complex128) ([][]complex128, error) {
	if len(a) == 0 || len(a) != len(a[0]) {
		return nil, fmt.Errorf("expm requires square matrix")
	}
	n := len(a)
	norm1 := 0.0
	for j := 0; j < n; j++ {
		colSum := 0.0
		for i := 0; i < n; i++ {
			colSum += cmplx.Abs(a[i][j])
		}
		if colSum > norm1 {
			norm1 = colSum
		}
	}
	s := 0
	for norm1 > 0.5 {
		norm1 /= 2.0
		s++
	}
	scaled := matScale(a, complex(math.Pow(2, -float64(s)), 0))

	sum := matIdentity(n)
	term := matIdentity(n)
	var err error
	for k := 1; k <= 24; k++ {
		term, err = matMul(term, scaled)
		if err != nil {
			return nil, err
		}
		term = matScale(term, complex(1.0/float64(k), 0))
		sum, _ = matAdd(sum, term, 1)
		if matNormScale(term) <= 1e-16*matNormScale(sum) {
			break
		}
	}

	for i := 0; i < s; i++ {
		sum, err = matMul(sum, sum)
		if err != nil {
			return nil, err
		}
	}
	return sum, nil
}

func vecDot(u, v []complex128) (complex128, error) {
	if len(u) != len(v) {
		return 0, fmt.Errorf("dot requires equal-length vectors")
	}
	var sum complex128
	for i := range u {
		sum += u[i] * cmplx.Conj(v[i])
	}
	return cleanComplex(sum), nil
}

func vecNorm(u []complex128, p float64) float64 {
	if math.IsInf(p, 1) {
		mx := 0.0
		for _, z := range u {
			if a := cmplx.Abs(z); a > mx {
				mx = a
			}
		}
		return mx
	}
	var sum float64
	for _, z := range u {
		sum += math.Pow(cmplx.Abs(z), p)
	}
	return math.Pow(sum, 1.0/p)
}

func evalLinAlgFunc(fn string, args []Value) (bool, Value, error) {
	switch fn {
	case "dot":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("dot(u, v) requires 2 vectors")
		}
		d, err := vecDot(args[0].V, args[1].V)
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarComplex(d), nil

	case "cross":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("cross(u, v) requires 2 3D vectors")
		}
		u, v := args[0].V, args[1].V
		if len(u) != 3 || len(v) != 3 {
			return true, Value{}, fmt.Errorf("cross product requires 3D vectors")
		}
		res := []complex128{
			u[1]*v[2] - u[2]*v[1],
			u[2]*v[0] - u[0]*v[2],
			u[0]*v[1] - u[1]*v[0],
		}
		return true, VectorVal(res), nil

	case "outer":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("outer(u, v) requires 2 vectors")
		}
		u, v := args[0].V, args[1].V
		m := matZeros(len(u), len(v))
		for i := range u {
			for j := range v {
				m[i][j] = cleanComplex(u[i] * cmplx.Conj(v[j]))
			}
		}
		return true, MatrixVal(m), nil

	case "norm":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("norm(x [, p]) takes 1 or 2 arguments")
		}
		p := 2.0
		if len(args) == 2 {
			pr, err := args[1].Real()
			if err != nil || pr <= 0 {
				return true, Value{}, fmt.Errorf("norm p must be positive real or inf")
			}
			p = pr
		}
		if args[0].Kind == KindVector {
			return true, ScalarReal(vecNorm(args[0].V, p)), nil
		}
		if args[0].Kind == KindMatrix {
			var sum float64
			for _, row := range args[0].M {
				for _, z := range row {
					a := cmplx.Abs(z)
					sum += a * a
				}
			}
			return true, ScalarReal(math.Sqrt(sum)), nil
		}
		if args[0].Kind == KindScalar {
			return true, ScalarReal(cmplx.Abs(args[0].S)), nil
		}
		return true, Value{}, fmt.Errorf("norm requires vector or matrix")

	case "unit", "normalize":
		if len(args) != 1 || args[0].Kind != KindVector {
			return true, Value{}, fmt.Errorf("unit(v) requires 1 vector")
		}
		nm := vecNorm(args[0].V, 2)
		if nm == 0 {
			return true, Value{}, fmt.Errorf("cannot normalize zero vector")
		}
		out := make([]complex128, len(args[0].V))
		for i, z := range args[0].V {
			out[i] = z / complex(nm, 0)
		}
		return true, VectorVal(out), nil

	case "proj":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("proj(u, v) requires 2 vectors")
		}
		uv, err := vecDot(args[0].V, args[1].V)
		if err != nil {
			return true, Value{}, err
		}
		vv, _ := vecDot(args[1].V, args[1].V)
		if vv == 0 {
			return true, Value{}, fmt.Errorf("cannot project onto zero vector")
		}
		scale := uv / vv
		out := make([]complex128, len(args[1].V))
		for i, z := range args[1].V {
			out[i] = scale * z
		}
		return true, VectorVal(out), nil

	case "angle":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("angle(u, v) requires 2 vectors")
		}
		uv, err := vecDot(args[0].V, args[1].V)
		if err != nil {
			return true, Value{}, err
		}
		nu := vecNorm(args[0].V, 2)
		nv := vecNorm(args[1].V, 2)
		if nu == 0 || nv == 0 {
			return true, Value{}, fmt.Errorf("angle undefined for zero vector")
		}
		cosTheta := math.Max(-1.0, math.Min(1.0, real(uv)/(nu*nv)))
		return true, ScalarReal(math.Acos(cosTheta)), nil

	case "dim", "size", "len":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("%s(x) takes 1 argument", fn)
		}
		if args[0].Kind == KindVector {
			return true, ScalarReal(float64(len(args[0].V))), nil
		}
		if args[0].Kind == KindMatrix {
			r := len(args[0].M)
			c := 0
			if r > 0 {
				c = len(args[0].M[0])
			}
			return true, VectorVal([]complex128{complex(float64(r), 0), complex(float64(c), 0)}), nil
		}
		return true, Value{}, fmt.Errorf("%s expects vector or matrix", fn)

	case "eye", "ident":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("eye(n) takes 1 integer argument")
		}
		nr, err := args[0].Real()
		if err != nil || nr < 1 || nr > 200 {
			return true, Value{}, fmt.Errorf("invalid dimension for eye")
		}
		return true, MatrixVal(matIdentity(int(math.Round(nr)))), nil

	case "zeros", "ones":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("%s(r [, c]) takes 1 or 2 arguments", fn)
		}
		rF, err := args[0].Real()
		if err != nil || rF < 1 || rF > 200 {
			return true, Value{}, fmt.Errorf("invalid row count")
		}
		r := int(math.Round(rF))
		c := r
		if len(args) == 2 {
			cF, err := args[1].Real()
			if err != nil || cF < 1 || cF > 200 {
				return true, Value{}, fmt.Errorf("invalid col count")
			}
			c = int(math.Round(cF))
		}
		m := matZeros(r, c)
		if fn == "ones" {
			for i := 0; i < r; i++ {
				for j := 0; j < c; j++ {
					m[i][j] = 1
				}
			}
		}
		return true, MatrixVal(m), nil

	case "diag":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("diag(v_or_M) takes 1 argument")
		}
		if args[0].Kind == KindVector {
			n := len(args[0].V)
			m := matZeros(n, n)
			for i, z := range args[0].V {
				m[i][i] = z
			}
			return true, MatrixVal(m), nil
		}
		if args[0].Kind == KindMatrix {
			r := len(args[0].M)
			if r == 0 {
				return true, VectorVal(nil), nil
			}
			k := r
			if len(args[0].M[0]) < k {
				k = len(args[0].M[0])
			}
			v := make([]complex128, k)
			for i := 0; i < k; i++ {
				v[i] = args[0].M[i][i]
			}
			return true, VectorVal(v), nil
		}
		return true, Value{}, fmt.Errorf("diag expects vector or matrix")

	case "hilbert":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("hilbert(n) takes 1 integer argument")
		}
		nr, err := args[0].Real()
		if err != nil || nr < 1 || nr > 50 {
			return true, Value{}, fmt.Errorf("hilbert dimension must be 1..50")
		}
		n := int(math.Round(nr))
		m := matZeros(n, n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				m[i][j] = complex(1.0/float64(i+j+1), 0)
			}
		}
		return true, MatrixVal(m), nil

	case "vandermonde":
		if len(args) < 1 || len(args) > 2 || args[0].Kind != KindVector {
			return true, Value{}, fmt.Errorf("vandermonde(v [, cols]) requires vector")
		}
		v := args[0].V
		cols := len(v)
		if len(args) == 2 {
			cr, err := args[1].Real()
			if err != nil || cr < 1 || cr > 100 {
				return true, Value{}, fmt.Errorf("invalid column count")
			}
			cols = int(math.Round(cr))
		}
		m := matZeros(len(v), cols)
		for i, z := range v {
			p := complex(1, 0)
			for j := 0; j < cols; j++ {
				m[i][j] = cleanComplex(p)
				p *= z
			}
		}
		return true, MatrixVal(m), nil

	case "trans", "transpose":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("trans(A) requires matrix")
		}
		return true, MatrixVal(matTranspose(args[0].M)), nil

	case "ctrans", "dagger", "hermitian":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("ctrans(A) requires matrix")
		}
		return true, MatrixVal(matConjTranspose(args[0].M)), nil

	case "trace", "tr":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("trace(A) requires square matrix")
		}
		m := args[0].M
		if len(m) == 0 || len(m) != len(m[0]) {
			return true, Value{}, fmt.Errorf("trace(A) requires square matrix")
		}
		var sum complex128
		for i := range m {
			sum += m[i][i]
		}
		return true, ScalarComplex(sum), nil

	case "det":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("det(A) requires square matrix")
		}
		d, err := matDet(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarComplex(d), nil

	case "inv":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("inv(A) requires square matrix")
		}
		inv, err := matInverse(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, MatrixVal(inv), nil

	case "pinv":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("pinv(A) requires matrix")
		}
		a := args[0].M
		r, c := len(a), len(a[0])
		ah := matConjTranspose(a)
		if r >= c {
			aha, _ := matMul(ah, a)
			for i := 0; i < c; i++ {
				aha[i][i] += 1e-13
			}
			inv, err := matInverse(aha)
			if err != nil {
				return true, Value{}, err
			}
			res, _ := matMul(inv, ah)
			return true, MatrixVal(res), nil
		}
		aah, _ := matMul(a, ah)
		for i := 0; i < r; i++ {
			aah[i][i] += 1e-13
		}
		inv, err := matInverse(aah)
		if err != nil {
			return true, Value{}, err
		}
		res, _ := matMul(ah, inv)
		return true, MatrixVal(res), nil

	case "adj", "adjugate":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("adj(A) requires square matrix")
		}
		_, adj, err := faddeevLeVerrier(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, MatrixVal(adj), nil

	case "rref":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("rref(A) requires matrix")
		}
		r, _ := matRREF(args[0].M, true)
		return true, MatrixVal(r), nil

	case "ref":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("ref(A) requires matrix")
		}
		r, _ := matRREF(args[0].M, false)
		return true, MatrixVal(r), nil

	case "rank":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("rank(A) requires matrix")
		}
		_, pivots := matRREF(args[0].M, false)
		return true, ScalarReal(float64(len(pivots))), nil

	case "nullity":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("nullity(A) requires matrix")
		}
		if len(args[0].M) == 0 {
			return true, ScalarReal(0), nil
		}
		_, pivots := matRREF(args[0].M, false)
		return true, ScalarReal(float64(len(args[0].M[0]) - len(pivots))), nil

	case "nullspace", "kernel":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("nullspace(A) requires matrix")
		}
		return true, MatrixVal(matNullspace(args[0].M)), nil

	case "colspace":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("colspace(A) requires matrix")
		}
		a := args[0].M
		if len(a) == 0 {
			return true, MatrixVal(nil), nil
		}
		_, pivots := matRREF(a, true)
		if len(pivots) == 0 {
			return true, MatrixVal(matZeros(len(a), 1)), nil
		}
		out := matZeros(len(a), len(pivots))
		for j, pCol := range pivots {
			for i := 0; i < len(a); i++ {
				out[i][j] = a[i][pCol]
			}
		}
		return true, MatrixVal(out), nil

	case "rowspace":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("rowspace(A) requires matrix")
		}
		rrefM, pivots := matRREF(args[0].M, true)
		if len(pivots) == 0 {
			return true, MatrixVal(matZeros(1, len(args[0].M[0]))), nil
		}
		return true, MatrixVal(rrefM[:len(pivots)]), nil

	case "linsolve":
		if len(args) != 2 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("linsolve(A, b) requires matrix A and vector/matrix b")
		}
		if args[1].Kind == KindVector {
			bCol := matZeros(len(args[1].V), 1)
			for i, z := range args[1].V {
				bCol[i][0] = z
			}
			sol, err := matLinearSolve(args[0].M, bCol)
			if err != nil {
				return true, Value{}, err
			}
			vec := make([]complex128, len(sol))
			for i := range sol {
				vec[i] = sol[i][0]
			}
			return true, VectorVal(vec), nil
		}
		if args[1].Kind == KindMatrix {
			sol, err := matLinearSolve(args[0].M, args[1].M)
			if err != nil {
				return true, Value{}, err
			}
			return true, MatrixVal(sol), nil
		}
		return true, Value{}, fmt.Errorf("b in linsolve(A, b) must be vector or matrix")

	case "eig", "eigenvals":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("eig(A) requires square matrix")
		}
		eigs, err := matEigenvalues(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, VectorVal(eigs), nil

	case "eigvecs":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("eigvecs(A) requires square matrix")
		}
		vecs, err := matEigenvectors(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, MatrixVal(vecs), nil

	case "charpoly":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("charpoly(A) requires square matrix")
		}
		coeffs, _, err := faddeevLeVerrier(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		realCoeffs := make([]float64, len(coeffs))
		allReal := true
		for i, z := range coeffs {
			if math.Abs(imag(z)) > 1e-11 {
				allReal = false
				break
			}
			realCoeffs[i] = real(z)
		}
		if allReal {
			astP := buildPolyAST(realCoeffs, "x")
			return true, SymbolicVal(astP.Format(), astP), nil
		}
		return true, VectorVal(coeffs), nil

	case "gramschmidt":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("gramschmidt(A) requires matrix")
		}
		return true, MatrixVal(matGramSchmidt(args[0].M)), nil

	case "expm":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("expm(A) requires square matrix")
		}
		em, err := matExpm(args[0].M)
		if err != nil {
			return true, Value{}, err
		}
		return true, MatrixVal(em), nil

	case "cond":
		if len(args) != 1 || args[0].Kind != KindMatrix {
			return true, Value{}, fmt.Errorf("cond(A) requires square matrix")
		}
		inv, err := matInverse(args[0].M)
		if err != nil {
			return true, ScalarReal(math.Inf(1)), nil
		}
		var nA, nInv float64
		for _, row := range args[0].M {
			for _, z := range row {
				nA += real(z)*real(z) + imag(z)*imag(z)
			}
		}
		for _, row := range inv {
			for _, z := range row {
				nInv += real(z)*real(z) + imag(z)*imag(z)
			}
		}
		return true, ScalarReal(math.Sqrt(nA * nInv)), nil
	}

	return false, Value{}, nil
}
