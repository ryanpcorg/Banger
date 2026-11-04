package main

import (
	"fmt"
	"math"
	"math/big"
	"sort"
	"strconv"
	"strings"
)

func requireInt64(v Value, name string) (int64, error) {
	r, err := v.Real()
	if err != nil {
		return 0, fmt.Errorf("%s requires integer argument", name)
	}
	rounded := math.Round(r)
	if math.Abs(r-rounded) > 1e-9 || math.Abs(rounded) > 9e18 {
		return 0, fmt.Errorf("%s requires 64-bit integer argument", name)
	}
	return int64(rounded), nil
}

func extractRealData(args []Value, fn string) ([]float64, error) {
	if len(args) == 0 {
		return nil, fmt.Errorf("%s requires at least 1 argument", fn)
	}
	if len(args) == 1 && args[0].Kind == KindVector {
		out := make([]float64, len(args[0].V))
		for i, z := range args[0].V {
			if math.Abs(imag(z)) > 1e-13 {
				return nil, fmt.Errorf("%s requires real data", fn)
			}
			out[i] = real(z)
		}
		if len(out) == 0 {
			return nil, fmt.Errorf("%s received empty vector", fn)
		}
		return out, nil
	}
	out := make([]float64, len(args))
	for i, a := range args {
		r, err := a.Real()
		if err != nil {
			return nil, fmt.Errorf("%s requires real data", fn)
		}
		out[i] = r
	}
	return out, nil
}

func isPrime64(n int64) bool {
	if n < 2 {
		return false
	}
	if n == 2 || n == 3 {
		return true
	}
	if n%2 == 0 || n%3 == 0 {
		return false
	}
	return big.NewInt(n).ProbablyPrime(0)
}

func primeFactorList(n int64) []int64 {
	if n < 0 {
		n = -n
	}
	if n < 2 {
		return nil
	}
	var factors []int64
	for n%2 == 0 {
		factors = append(factors, 2)
		n /= 2
	}
	for n%3 == 0 {
		factors = append(factors, 3)
		n /= 3
	}
	for d := int64(5); d*d <= n; d += 6 {
		for n%d == 0 {
			factors = append(factors, d)
			n /= d
		}
		for n%(d+2) == 0 {
			factors = append(factors, d+2)
			n /= (d + 2)
		}
	}
	if n > 1 {
		factors = append(factors, n)
	}
	return factors
}

func formatPrimeFactorization(n int64) string {
	if n == 0 {
		return "0"
	}
	if n == 1 {
		return "1"
	}
	if n == -1 {
		return "-1"
	}
	prefix := ""
	if n < 0 {
		prefix = "-1 * "
		n = -n
	}
	pf := primeFactorList(n)
	var parts []string
	i := 0
	for i < len(pf) {
		p := pf[i]
		cnt := 1
		for i+cnt < len(pf) && pf[i+cnt] == p {
			cnt++
		}
		if cnt == 1 {
			parts = append(parts, strconv.FormatInt(p, 10))
		} else {
			parts = append(parts, fmt.Sprintf("%d^%d", p, cnt))
		}
		i += cnt
	}
	return prefix + strings.Join(parts, " * ")
}

func eulerTotient(n int64) int64 {
	if n <= 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	res := n
	temp := n
	for d := int64(2); d*d <= temp; d++ {
		if temp%d == 0 {
			for temp%d == 0 {
				temp /= d
			}
			res -= res / d
		}
	}
	if temp > 1 {
		res -= res / temp
	}
	return res
}

func mobiusMu(n int64) int64 {
	if n <= 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	pf := primeFactorList(n)
	for i := 1; i < len(pf); i++ {
		if pf[i] == pf[i-1] {
			return 0
		}
	}
	if len(pf)%2 == 0 {
		return 1
	}
	return -1
}

func extendedGCD(a, b int64) (int64, int64, int64) {
	oldR, r := a, b
	oldS, s := int64(1), int64(0)
	oldT, t := int64(0), int64(1)
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldS, s = s, oldS-q*s
		oldT, t = t, oldT-q*t
	}
	if oldR < 0 {
		return -oldR, -oldS, -oldT
	}
	return oldR, oldS, oldT
}

func solveCRT(rems, mods []int64) (int64, int64, error) {
	if len(rems) == 0 || len(rems) != len(mods) {
		return 0, 0, fmt.Errorf("crt requires equal-length remainder and modulus vectors")
	}
	r0 := big.NewInt(rems[0])
	m0 := big.NewInt(mods[0])
	if m0.Sign() <= 0 {
		return 0, 0, fmt.Errorf("moduli in crt must be positive")
	}
	r0.Mod(r0, m0)

	for i := 1; i < len(rems); i++ {
		r1 := big.NewInt(rems[i])
		m1 := big.NewInt(mods[i])
		if m1.Sign() <= 0 {
			return 0, 0, fmt.Errorf("moduli in crt must be positive")
		}
		r1.Mod(r1, m1)

		g := new(big.Int)
		p := new(big.Int)
		q := new(big.Int)
		g.GCD(p, q, m0, m1)

		diff := new(big.Int).Sub(r1, r0)
		remCheck := new(big.Int).Mod(diff, g)
		if remCheck.Sign() != 0 {
			return 0, 0, fmt.Errorf("system of congruences has no solution (incompatible moduli)")
		}

		lcm := new(big.Int).Mul(new(big.Int).Quo(m0, g), m1)
		step := new(big.Int).Quo(diff, g)
		term := new(big.Int).Mul(m0, new(big.Int).Mul(step, p))
		r0.Add(r0, term)
		r0.Mod(r0, lcm)
		m0 = lcm
	}
	if !r0.IsInt64() {
		return 0, 0, fmt.Errorf("crt result overflows 64-bit integer")
	}
	modOut := int64(0)
	if m0.IsInt64() {
		modOut = m0.Int64()
	}
	return r0.Int64(), modOut, nil
}

func jacobiSymbol(a, n int64) (int64, error) {
	if n <= 0 || n%2 == 0 {
		return 0, fmt.Errorf("jacobi(a, n) requires odd positive n")
	}
	a %= n
	if a < 0 {
		a += n
	}
	t := int64(1)
	for a != 0 {
		for a%2 == 0 {
			a /= 2
			r := n % 8
			if r == 3 || r == 5 {
				t = -t
			}
		}
		a, n = n, a
		if a%4 == 3 && n%4 == 3 {
			t = -t
		}
		a %= n
	}
	if n == 1 {
		return t, nil
	}
	return 0, nil
}

func nCrExact(n, r float64) float64 {
	if n < 0 || r < 0 || r > n {
		return 0
	}
	if n == math.Trunc(n) && r == math.Trunc(r) {
		ni, ri := int64(n), int64(r)
		if ri > ni-ri {
			ri = ni - ri
		}
		if ri == 0 {
			return 1
		}
		b := new(big.Int).Binomial(ni, ri)
		f, _ := new(big.Float).SetInt(b).Float64()
		return f
	}
	return math.Gamma(n+1) / (math.Gamma(r+1) * math.Gamma(n-r+1))
}

func nPrExact(n, r float64) float64 {
	if n < 0 || r < 0 || r > n {
		return 0
	}
	if n == math.Trunc(n) && r == math.Trunc(r) {
		ni, ri := int64(n), int64(r)
		res := 1.0
		for k := int64(0); k < ri; k++ {
			res *= float64(ni - k)
			if math.IsInf(res, 1) {
				return res
			}
		}
		return res
	}
	return math.Gamma(n+1) / math.Gamma(n-r+1)
}

func partitionP(n int) float64 {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	dp := make([]float64, n+1)
	dp[0] = 1
	for i := 1; i <= n; i++ {
		var sum float64
		for k := 1; ; k++ {
			g1 := k * (3*k - 1) / 2
			g2 := k * (3*k + 1) / 2
			if g1 > i {
				break
			}
			sign := 1.0
			if k%2 == 0 {
				sign = -1.0
			}
			sum += sign * dp[i-g1]
			if g2 <= i {
				sum += sign * dp[i-g2]
			}
		}
		dp[i] = sum
	}
	return dp[n]
}

func stirlingSecond(n, k int) float64 {
	if n < 0 || k < 0 || k > n {
		return 0
	}
	if n == 0 && k == 0 {
		return 1
	}
	if k == 0 {
		return 0
	}
	dp := make([][]float64, n+1)
	for i := 0; i <= n; i++ {
		dp[i] = make([]float64, k+1)
	}
	dp[0][0] = 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= i && j <= k; j++ {
			dp[i][j] = float64(j)*dp[i-1][j] + dp[i-1][j-1]
		}
	}
	return dp[n][k]
}

func stirlingFirstUnsigned(n, k int) float64 {
	if n < 0 || k < 0 || k > n {
		return 0
	}
	if n == 0 && k == 0 {
		return 1
	}
	if k == 0 {
		return 0
	}
	dp := make([][]float64, n+1)
	for i := 0; i <= n; i++ {
		dp[i] = make([]float64, k+1)
	}
	dp[0][0] = 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= i && j <= k; j++ {
			dp[i][j] = float64(i-1)*dp[i-1][j] + dp[i-1][j-1]
		}
	}
	return dp[n][k]
}

func quantileSorted(sorted []float64, p float64) float64 {
	n := len(sorted)
	if n == 0 {
		return math.NaN()
	}
	if p <= 0 {
		return sorted[0]
	}
	if p >= 1 {
		return sorted[n-1]
	}
	idx := p * float64(n-1)
	lo := int(math.Floor(idx))
	hi := int(math.Ceil(idx))
	if lo == hi {
		return sorted[lo]
	}
	frac := idx - float64(lo)
	return sorted[lo]*(1-frac) + sorted[hi]*frac
}

func normCDF(x, mu, sigma float64) float64 {
	z := (x - mu) / (sigma * math.Sqrt2)
	return 0.5 * (1 + math.Erf(z))
}

func normInv(p, mu, sigma float64) float64 {
	if p <= 0 {
		if p == 0 {
			return math.Inf(-1)
		}
		return math.NaN()
	}
	if p >= 1 {
		if p == 1 {
			return math.Inf(1)
		}
		return math.NaN()
	}
	z := math.Sqrt2 * erfInv(2*p-1)
	return mu + sigma*z
}

func studentTCDF(t, df float64) float64 {
	x := df / (df + t*t)
	ib := regBetaInc(0.5*df, 0.5, x)
	if t >= 0 {
		return 1.0 - 0.5*ib
	}
	return 0.5 * ib
}

func invertMonotonicCDF(cdf func(float64) float64, p, lo, hi float64) float64 {
	for i := 0; i < 60; i++ {
		if cdf(lo) > p {
			lo -= (hi - lo)
		} else {
			break
		}
	}
	for i := 0; i < 60; i++ {
		if cdf(hi) < p {
			hi += (hi - lo)
		} else {
			break
		}
	}
	for i := 0; i < 100; i++ {
		mid := 0.5 * (lo + hi)
		v := cdf(mid)
		if math.Abs(v-p) <= 1e-13 || math.Abs(hi-lo) <= 1e-13*(1+math.Abs(mid)) {
			return mid
		}
		if v < p {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

func evalNumStatsFunc(fn string, args []Value) (bool, Value, error) {
	switch fn {
	case "gcd":
		data, err := extractRealData(args, "gcd")
		if err != nil {
			return true, Value{}, err
		}
		var g int64
		for _, v := range data {
			g = gcdInt64(g, int64(math.Round(v)))
		}
		return true, ScalarReal(float64(g)), nil

	case "lcm":
		data, err := extractRealData(args, "lcm")
		if err != nil {
			return true, Value{}, err
		}
		if len(data) == 0 {
			return true, ScalarReal(0), nil
		}
		res := absInt64(int64(math.Round(data[0])))
		for i := 1; i < len(data); i++ {
			b := absInt64(int64(math.Round(data[i])))
			if res == 0 || b == 0 {
				res = 0
			} else {
				res = (res / gcdInt64(res, b)) * b
			}
		}
		return true, ScalarReal(float64(res)), nil

	case "xgcd":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("xgcd(a, b) requires 2 integers")
		}
		a, err1 := requireInt64(args[0], "xgcd")
		b, err2 := requireInt64(args[1], "xgcd")
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("xgcd requires integer arguments")
		}
		g, x, y := extendedGCD(a, b)
		return true, VectorVal([]complex128{complex(float64(g), 0), complex(float64(x), 0), complex(float64(y), 0)}), nil

	case "isprime":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("isprime(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "isprime")
		if err != nil {
			return true, Value{}, err
		}
		if isPrime64(n) {
			return true, ScalarReal(1), nil
		}
		return true, ScalarReal(0), nil

	case "nextprime":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("nextprime(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "nextprime")
		if err != nil {
			return true, Value{}, err
		}
		if n < 2 {
			return true, ScalarReal(2), nil
		}
		cand := n + 1
		if cand%2 == 0 {
			cand++
		}
		for !isPrime64(cand) {
			cand += 2
		}
		return true, ScalarReal(float64(cand)), nil

	case "prevprime":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("prevprime(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "prevprime")
		if err != nil || n <= 2 {
			return true, Value{}, fmt.Errorf("prevprime requires integer > 2")
		}
		if n == 3 {
			return true, ScalarReal(2), nil
		}
		cand := n - 1
		if cand%2 == 0 {
			cand--
		}
		for cand >= 2 && !isPrime64(cand) {
			cand -= 2
		}
		return true, ScalarReal(float64(cand)), nil

	case "nthprime":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("nthprime(k) takes 1 positive integer")
		}
		k, err := requireInt64(args[0], "nthprime")
		if err != nil || k < 1 || k > 200000 {
			return true, Value{}, fmt.Errorf("nthprime requires k between 1 and 200000")
		}
		if k == 1 {
			return true, ScalarReal(2), nil
		}
		cnt := int64(1)
		cand := int64(3)
		for {
			if isPrime64(cand) {
				cnt++
				if cnt == k {
					return true, ScalarReal(float64(cand)), nil
				}
			}
			cand += 2
		}

	case "primepi":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("primepi(x) takes 1 argument")
		}
		xr, err := args[0].Real()
		if err != nil || xr > 1e7 {
			return true, Value{}, fmt.Errorf("primepi requires real x <= 10,000,000")
		}
		limit := int(math.Floor(xr))
		if limit < 2 {
			return true, ScalarReal(0), nil
		}
		composite := make([]bool, limit+1)
		for i := 2; i*i <= limit; i++ {
			if !composite[i] {
				for j := i * i; j <= limit; j += i {
					composite[j] = true
				}
			}
		}
		count := 0
		for i := 2; i <= limit; i++ {
			if !composite[i] {
				count++
			}
		}
		return true, ScalarReal(float64(count)), nil

	case "factor", "factors":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("factor(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "factor")
		if err != nil {
			return true, Value{}, err
		}
		return true, SymbolicVal(formatPrimeFactorization(n), nil), nil

	case "primefactors":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("primefactors(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "primefactors")
		if err != nil {
			return true, Value{}, err
		}
		pf := primeFactorList(n)
		out := make([]complex128, len(pf))
		for i, p := range pf {
			out[i] = complex(float64(p), 0)
		}
		return true, VectorVal(out), nil

	case "divisors":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("divisors(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "divisors")
		if err != nil || n == 0 {
			return true, Value{}, fmt.Errorf("divisors requires non-zero integer")
		}
		if n < 0 {
			n = -n
		}
		var divs []int64
		for d := int64(1); d*d <= n; d++ {
			if n%d == 0 {
				divs = append(divs, d)
				if d*d != n {
					divs = append(divs, n/d)
				}
			}
		}
		sort.Slice(divs, func(i, j int) bool { return divs[i] < divs[j] })
		out := make([]complex128, len(divs))
		for i, d := range divs {
			out[i] = complex(float64(d), 0)
		}
		return true, VectorVal(out), nil

	case "sigma":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("sigma(n [, k]) takes 1 or 2 arguments")
		}
		n, err := requireInt64(args[0], "sigma")
		if err != nil || n <= 0 {
			return true, Value{}, fmt.Errorf("sigma requires positive integer n")
		}
		k := 1.0
		if len(args) == 2 {
			kr, err := args[1].Real()
			if err != nil {
				return true, Value{}, err
			}
			k = kr
		}
		var sum float64
		for d := int64(1); d*d <= n; d++ {
			if n%d == 0 {
				sum += math.Pow(float64(d), k)
				if d*d != n {
					sum += math.Pow(float64(n/d), k)
				}
			}
		}
		return true, ScalarReal(sum), nil

	case "totient", "eulerphi":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("totient(n) takes 1 positive integer")
		}
		n, err := requireInt64(args[0], "totient")
		if err != nil || n <= 0 {
			return true, Value{}, fmt.Errorf("totient requires positive integer")
		}
		return true, ScalarReal(float64(eulerTotient(n))), nil

	case "mobius":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("mobius(n) takes 1 positive integer")
		}
		n, err := requireInt64(args[0], "mobius")
		if err != nil || n <= 0 {
			return true, Value{}, fmt.Errorf("mobius requires positive integer")
		}
		return true, ScalarReal(float64(mobiusMu(n))), nil

	case "modpow":
		if len(args) != 3 {
			return true, Value{}, fmt.Errorf("modpow(base, exp, mod) requires 3 integers")
		}
		b, err1 := requireInt64(args[0], "modpow")
		e, err2 := requireInt64(args[1], "modpow")
		m, err3 := requireInt64(args[2], "modpow")
		if err1 != nil || err2 != nil || err3 != nil || m <= 0 {
			return true, Value{}, fmt.Errorf("modpow requires integer base, exp, and positive mod")
		}
		bb := big.NewInt(b)
		bm := big.NewInt(m)
		if e < 0 {
			inv := new(big.Int).ModInverse(bb, bm)
			if inv == nil {
				return true, Value{}, fmt.Errorf("%d has no modular inverse mod %d", b, m)
			}
			bb = inv
			e = -e
		}
		res := new(big.Int).Exp(bb, big.NewInt(e), bm)
		return true, ScalarReal(float64(res.Int64())), nil

	case "modinv":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("modinv(a, m) requires 2 integers")
		}
		a, err1 := requireInt64(args[0], "modinv")
		m, err2 := requireInt64(args[1], "modinv")
		if err1 != nil || err2 != nil || m <= 0 {
			return true, Value{}, fmt.Errorf("modinv requires integer a and positive modulus m")
		}
		inv := new(big.Int).ModInverse(big.NewInt(a), big.NewInt(m))
		if inv == nil {
			return true, Value{}, fmt.Errorf("%d is not invertible mod %d (gcd != 1)", a, m)
		}
		return true, ScalarReal(float64(inv.Int64())), nil

	case "crt":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("crt([r1, r2, ...], [m1, m2, ...]) requires 2 vectors")
		}
		rems := make([]int64, len(args[0].V))
		mods := make([]int64, len(args[1].V))
		if len(rems) != len(mods) {
			return true, Value{}, fmt.Errorf("crt vectors must have equal length")
		}
		for i := range rems {
			r, err1 := ScalarComplex(args[0].V[i]).Real()
			m, err2 := ScalarComplex(args[1].V[i]).Real()
			if err1 != nil || err2 != nil {
				return true, Value{}, fmt.Errorf("crt requires integer vectors")
			}
			rems[i] = int64(math.Round(r))
			mods[i] = int64(math.Round(m))
		}
		sol, _, err := solveCRT(rems, mods)
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(float64(sol)), nil

	case "znorder":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("znorder(a, n) requires 2 integers")
		}
		a, err1 := requireInt64(args[0], "znorder")
		n, err2 := requireInt64(args[1], "znorder")
		if err1 != nil || err2 != nil || n <= 1 {
			return true, Value{}, fmt.Errorf("znorder requires n > 1")
		}
		if gcdInt64(a, n) != 1 {
			return true, Value{}, fmt.Errorf("znorder requires gcd(a, n) = 1")
		}
		phi := eulerTotient(n)
		var divs []int64
		for d := int64(1); d*d <= phi; d++ {
			if phi%d == 0 {
				divs = append(divs, d)
				if d*d != phi {
					divs = append(divs, phi/d)
				}
			}
		}
		sort.Slice(divs, func(i, j int) bool { return divs[i] < divs[j] })
		ba, bn := big.NewInt(a), big.NewInt(n)
		for _, d := range divs {
			if new(big.Int).Exp(ba, big.NewInt(d), bn).Cmp(big.NewInt(1)) == 0 {
				return true, ScalarReal(float64(d)), nil
			}
		}
		return true, ScalarReal(float64(phi)), nil

	case "primroot":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("primroot(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "primroot")
		if err != nil || n <= 1 {
			return true, Value{}, fmt.Errorf("primroot requires n > 1")
		}
		if n == 2 {
			return true, ScalarReal(1), nil
		}
		phi := eulerTotient(n)
		pf := primeFactorList(phi)
		var distinctPF []int64
		for i, p := range pf {
			if i == 0 || p != pf[i-1] {
				distinctPF = append(distinctPF, p)
			}
		}
		bn := big.NewInt(n)
		for g := int64(2); g < n; g++ {
			if gcdInt64(g, n) != 1 {
				continue
			}
			bg := big.NewInt(g)
			isRoot := true
			for _, q := range distinctPF {
				if new(big.Int).Exp(bg, big.NewInt(phi/q), bn).Cmp(big.NewInt(1)) == 0 {
					isRoot = false
					break
				}
			}
			if isRoot {
				return true, ScalarReal(float64(g)), nil
			}
		}
		return true, Value{}, fmt.Errorf("no primitive root exists modulo %d", n)

	case "jacobi", "legendre":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(a, n) requires 2 integers", fn)
		}
		a, err1 := requireInt64(args[0], fn)
		n, err2 := requireInt64(args[1], fn)
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("%s requires integer arguments", fn)
		}
		if fn == "legendre" && !isPrime64(n) {
			return true, Value{}, fmt.Errorf("legendre(a, p) requires odd prime p")
		}
		j, err := jacobiSymbol(a, n)
		if err != nil {
			return true, Value{}, err
		}
		return true, ScalarReal(float64(j)), nil

	case "ncr", "binom", "choose":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("ncr(n, r) requires 2 real arguments")
		}
		n, err1 := args[0].Real()
		r, err2 := args[1].Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("ncr requires real arguments")
		}
		return true, ScalarReal(nCrExact(n, r)), nil

	case "npr", "perm":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("npr(n, r) requires 2 real arguments")
		}
		n, err1 := args[0].Real()
		r, err2 := args[1].Real()
		if err1 != nil || err2 != nil {
			return true, Value{}, fmt.Errorf("npr requires real arguments")
		}
		return true, ScalarReal(nPrExact(n, r)), nil

	case "fib":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("fib(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "fib")
		if err != nil || absInt64(n) > 1500 {
			return true, Value{}, fmt.Errorf("fib requires integer with |n| <= 1500")
		}
		neg := n < 0
		k := absInt64(n)
		a, b := big.NewInt(0), big.NewInt(1)
		for i := int64(0); i < k; i++ {
			a, b = b, new(big.Int).Add(a, b)
		}
		if neg && k%2 == 0 {
			a.Neg(a)
		}
		f, _ := new(big.Float).SetInt(a).Float64()
		return true, ScalarReal(f), nil

	case "lucas":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("lucas(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "lucas")
		if err != nil || n < 0 || n > 1500 {
			return true, Value{}, fmt.Errorf("lucas requires 0 <= n <= 1500")
		}
		if n == 0 {
			return true, ScalarReal(2), nil
		}
		a, b := 2.0, 1.0
		for i := int64(1); i < n; i++ {
			a, b = b, a+b
		}
		return true, ScalarReal(b), nil

	case "catalannum", "catalan_n":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("catalannum(n) takes 1 non-negative integer")
		}
		n, err := requireInt64(args[0], "catalannum")
		if err != nil || n < 0 {
			return true, Value{}, fmt.Errorf("catalannum requires non-negative integer")
		}
		return true, ScalarReal(nCrExact(float64(2*n), float64(n)) / float64(n+1)), nil

	case "stirling1":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("stirling1(n, k) requires 2 integers")
		}
		n, err1 := requireInt64(args[0], "stirling1")
		k, err2 := requireInt64(args[1], "stirling1")
		if err1 != nil || err2 != nil || n < 0 || n > 500 {
			return true, Value{}, fmt.Errorf("stirling1 requires 0 <= n <= 500")
		}
		return true, ScalarReal(stirlingFirstUnsigned(int(n), int(k))), nil

	case "stirling2":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("stirling2(n, k) requires 2 integers")
		}
		n, err1 := requireInt64(args[0], "stirling2")
		k, err2 := requireInt64(args[1], "stirling2")
		if err1 != nil || err2 != nil || n < 0 || n > 500 {
			return true, Value{}, fmt.Errorf("stirling2 requires 0 <= n <= 500")
		}
		return true, ScalarReal(stirlingSecond(int(n), int(k))), nil

	case "bell":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("bell(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "bell")
		if err != nil || n < 0 || n > 300 {
			return true, Value{}, fmt.Errorf("bell requires 0 <= n <= 300")
		}
		if n == 0 {
			return true, ScalarReal(1), nil
		}
		var sum float64
		for k := 0; k <= int(n); k++ {
			sum += stirlingSecond(int(n), k)
		}
		return true, ScalarReal(sum), nil

	case "partition":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("partition(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "partition")
		if err != nil || n < 0 || n > 1000 {
			return true, Value{}, fmt.Errorf("partition requires 0 <= n <= 1000")
		}
		return true, ScalarReal(partitionP(int(n))), nil

	case "subfact", "derange":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("subfact(n) takes 1 integer")
		}
		n, err := requireInt64(args[0], "subfact")
		if err != nil || n < 0 || n > 170 {
			return true, Value{}, fmt.Errorf("subfact requires 0 <= n <= 170")
		}
		if n == 0 {
			return true, ScalarReal(1), nil
		}
		if n == 1 {
			return true, ScalarReal(0), nil
		}
		a, b := 1.0, 0.0
		for i := int64(2); i <= n; i++ {
			a, b = b, float64(i-1)*(a+b)
		}
		return true, ScalarReal(b), nil

	case "frac", "rat":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("frac(x [, maxDen]) takes 1 or 2 arguments")
		}
		x, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		maxDen := int64(1000000)
		if len(args) == 2 {
			md, err := requireInt64(args[1], "frac")
			if err != nil || md < 1 {
				return true, Value{}, fmt.Errorf("maxDen must be positive integer")
			}
			maxDen = md
		}
		num, den, ok := bestRational(x, maxDen, 1e-12)
		if !ok {
			return true, SymbolicVal(formatFloatClean(x), nil), nil
		}
		if den == 1 {
			return true, SymbolicVal(strconv.FormatInt(num, 10), nil), nil
		}
		return true, SymbolicVal(fmt.Sprintf("%d/%d", num, den), nil), nil

	case "cf":
		if len(args) < 1 || len(args) > 2 {
			return true, Value{}, fmt.Errorf("cf(x [, terms]) takes 1 or 2 arguments")
		}
		x, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		maxTerms := 12
		if len(args) == 2 {
			mt, err := requireInt64(args[1], "cf")
			if err != nil || mt < 1 || mt > 50 {
				return true, Value{}, fmt.Errorf("terms must be 1..50")
			}
			maxTerms = int(mt)
		}
		var out []complex128
		b := x
		for i := 0; i < maxTerms; i++ {
			a := math.Floor(b + 1e-13)
			out = append(out, complex(a, 0))
			rem := b - a
			if math.Abs(rem) < 1e-11 {
				break
			}
			b = 1.0 / rem
			if math.Abs(b) > 1e11 {
				break
			}
		}
		return true, VectorVal(out), nil

	case "hex", "bin", "oct":
		if len(args) != 1 {
			return true, Value{}, fmt.Errorf("%s(n) takes 1 integer", fn)
		}
		n, err := requireInt64(args[0], fn)
		if err != nil {
			return true, Value{}, err
		}
		sign := ""
		if n < 0 {
			sign = "-"
			n = -n
		}
		switch fn {
		case "hex":
			return true, SymbolicVal(sign+"0x"+strings.ToUpper(strconv.FormatInt(n, 16)), nil), nil
		case "bin":
			return true, SymbolicVal(sign+"0b"+strconv.FormatInt(n, 2), nil), nil
		case "oct":
			return true, SymbolicVal(sign+"0o"+strconv.FormatInt(n, 8), nil), nil
		}

	case "base":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("base(n, b) requires integer n and base 2..36")
		}
		n, err1 := requireInt64(args[0], "base")
		b, err2 := requireInt64(args[1], "base")
		if err1 != nil || err2 != nil || b < 2 || b > 36 {
			return true, Value{}, fmt.Errorf("base requires integer n and base between 2 and 36")
		}
		return true, SymbolicVal(strings.ToUpper(strconv.FormatInt(n, int(b))), nil), nil

	case "min", "max", "sum", "prod", "mean", "avg", "geomean", "harmmean", "rms", "median", "mode", "var", "pvar", "std", "stdev", "pstd", "iqr", "mad", "skewness", "kurtosis", "sort", "zscore":
		data, err := extractRealData(args, fn)
		if err != nil {
			return true, Value{}, err
		}
		n := len(data)
		switch fn {
		case "min":
			m := data[0]
			for _, x := range data[1:] {
				if x < m {
					m = x
				}
			}
			return true, ScalarReal(m), nil

		case "max":
			m := data[0]
			for _, x := range data[1:] {
				if x > m {
					m = x
				}
			}
			return true, ScalarReal(m), nil

		case "sum":
			var s float64
			for _, x := range data {
				s += x
			}
			return true, ScalarReal(s), nil

		case "prod":
			p := 1.0
			for _, x := range data {
				p *= x
			}
			return true, ScalarReal(p), nil

		case "sort":
			cp := append([]float64(nil), data...)
			sort.Float64s(cp)
			out := make([]complex128, n)
			for i, x := range cp {
				out[i] = complex(x, 0)
			}
			return true, VectorVal(out), nil

		case "mean", "avg":
			var s float64
			for _, x := range data {
				s += x
			}
			return true, ScalarReal(s / float64(n)), nil

		case "geomean":
			var logSum float64
			for _, x := range data {
				if x <= 0 {
					return true, Value{}, fmt.Errorf("geomean requires positive numbers")
				}
				logSum += math.Log(x)
			}
			return true, ScalarReal(math.Exp(logSum / float64(n))), nil

		case "harmmean":
			var invSum float64
			for _, x := range data {
				if x == 0 {
					return true, Value{}, fmt.Errorf("harmmean undefined for zero element")
				}
				invSum += 1.0 / x
			}
			return true, ScalarReal(float64(n) / invSum), nil

		case "rms":
			var sqSum float64
			for _, x := range data {
				sqSum += x * x
			}
			return true, ScalarReal(math.Sqrt(sqSum / float64(n))), nil

		case "median":
			cp := append([]float64(nil), data...)
			sort.Float64s(cp)
			return true, ScalarReal(quantileSorted(cp, 0.5)), nil

		case "mode":
			counts := make(map[float64]int)
			bestVal := data[0]
			bestCnt := 0
			for _, x := range data {
				counts[x]++
				if counts[x] > bestCnt {
					bestCnt = counts[x]
					bestVal = x
				}
			}
			return true, ScalarReal(bestVal), nil

		case "var", "pvar", "std", "stdev", "pstd", "zscore", "skewness", "kurtosis":
			var s float64
			for _, x := range data {
				s += x
			}
			mu := s / float64(n)
			var m2, m3, m4 float64
			for _, x := range data {
				d := x - mu
				d2 := d * d
				m2 += d2
				m3 += d2 * d
				m4 += d2 * d2
			}
			switch fn {
			case "pvar":
				return true, ScalarReal(m2 / float64(n)), nil
			case "pstd":
				return true, ScalarReal(math.Sqrt(m2 / float64(n))), nil
			case "var":
				if n < 2 {
					return true, Value{}, fmt.Errorf("sample variance requires at least 2 values")
				}
				return true, ScalarReal(m2 / float64(n-1)), nil
			case "std", "stdev":
				if n < 2 {
					return true, Value{}, fmt.Errorf("sample std requires at least 2 values")
				}
				return true, ScalarReal(math.Sqrt(m2 / float64(n-1))), nil
			case "zscore":
				if n < 2 || m2 == 0 {
					return true, Value{}, fmt.Errorf("zscore requires at least 2 non-identical values")
				}
				sd := math.Sqrt(m2 / float64(n-1))
				out := make([]complex128, n)
				for i, x := range data {
					out[i] = complex((x-mu)/sd, 0)
				}
				return true, VectorVal(out), nil
			case "skewness":
				sdPop := math.Sqrt(m2 / float64(n))
				if sdPop == 0 {
					return true, ScalarReal(0), nil
				}
				return true, ScalarReal((m3 / float64(n)) / (sdPop * sdPop * sdPop)), nil
			case "kurtosis":
				varPop := m2 / float64(n)
				if varPop == 0 {
					return true, ScalarReal(0), nil
				}
				return true, ScalarReal((m4/float64(n))/(varPop*varPop) - 3.0), nil
			}

		case "iqr":
			cp := append([]float64(nil), data...)
			sort.Float64s(cp)
			return true, ScalarReal(quantileSorted(cp, 0.75) - quantileSorted(cp, 0.25)), nil

		case "mad":
			cp := append([]float64(nil), data...)
			sort.Float64s(cp)
			med := quantileSorted(cp, 0.5)
			devs := make([]float64, n)
			for i, x := range cp {
				devs[i] = math.Abs(x - med)
			}
			sort.Float64s(devs)
			return true, ScalarReal(quantileSorted(devs, 0.5)), nil
		}

	case "quantile", "percentile":
		if len(args) != 2 || args[0].Kind != KindVector {
			return true, Value{}, fmt.Errorf("%s(vector, p) requires vector and scalar p", fn)
		}
		data, err := extractRealData(args[:1], fn)
		if err != nil {
			return true, Value{}, err
		}
		p, err := args[1].Real()
		if err != nil {
			return true, Value{}, err
		}
		if fn == "percentile" {
			p /= 100.0
		}
		if p < 0 || p > 1 {
			return true, Value{}, fmt.Errorf("quantile probability must be in [0, 1]")
		}
		sort.Float64s(data)
		return true, ScalarReal(quantileSorted(data, p)), nil

	case "cov", "corr", "linreg":
		if len(args) != 2 || args[0].Kind != KindVector || args[1].Kind != KindVector {
			return true, Value{}, fmt.Errorf("%s(X, Y) requires 2 equal-length vectors", fn)
		}
		xs, err1 := extractRealData(args[:1], fn)
		ys, err2 := extractRealData(args[1:2], fn)
		if err1 != nil || err2 != nil || len(xs) != len(ys) || len(xs) < 2 {
			return true, Value{}, fmt.Errorf("%s requires 2 equal-length real vectors with len >= 2", fn)
		}
		n := float64(len(xs))
		var sx, sy float64
		for i := range xs {
			sx += xs[i]
			sy += ys[i]
		}
		mx, my := sx/n, sy/n
		var sxx, syy, sxy float64
		for i := range xs {
			dx := xs[i] - mx
			dy := ys[i] - my
			sxx += dx * dx
			syy += dy * dy
			sxy += dx * dy
		}
		if fn == "cov" {
			return true, ScalarReal(sxy / (n - 1)), nil
		}
		if fn == "corr" {
			if sxx == 0 || syy == 0 {
				return true, ScalarReal(0), nil
			}
			return true, ScalarReal(sxy / math.Sqrt(sxx*syy)), nil
		}
		if sxx == 0 {
			return true, Value{}, fmt.Errorf("linreg undefined when X has zero variance")
		}
		slope := sxy / sxx
		intercept := my - slope*mx
		r2 := 1.0
		if syy > 0 {
			r := sxy / math.Sqrt(sxx*syy)
			r2 = r * r
		}
		return true, VectorVal([]complex128{complex(slope, 0), complex(intercept, 0), complex(r2, 0)}), nil

	case "normpdf", "normcdf", "norminv":
		if len(args) != 1 && len(args) != 3 {
			return true, Value{}, fmt.Errorf("%s(x [, mu, sigma]) takes 1 or 3 arguments", fn)
		}
		x, err := args[0].Real()
		if err != nil {
			return true, Value{}, err
		}
		mu, sigma := 0.0, 1.0
		if len(args) == 3 {
			mu, err = args[1].Real()
			if err != nil {
				return true, Value{}, err
			}
			sigma, err = args[2].Real()
			if err != nil || sigma <= 0 {
				return true, Value{}, fmt.Errorf("sigma must be positive")
			}
		}
		switch fn {
		case "normpdf":
			z := (x - mu) / sigma
			return true, ScalarReal(math.Exp(-0.5*z*z) / (sigma * math.Sqrt(2*math.Pi))), nil
		case "normcdf":
			return true, ScalarReal(normCDF(x, mu, sigma)), nil
		case "norminv":
			return true, ScalarReal(normInv(x, mu, sigma)), nil
		}

	case "tpdf", "tcdf", "tinv":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(x, df) requires 2 real arguments", fn)
		}
		x, err1 := args[0].Real()
		df, err2 := args[1].Real()
		if err1 != nil || err2 != nil || df <= 0 {
			return true, Value{}, fmt.Errorf("%s requires real x and df > 0", fn)
		}
		switch fn {
		case "tpdf":
			lg1, _ := math.Lgamma(0.5 * (df + 1))
			lg2, _ := math.Lgamma(0.5 * df)
			num := math.Exp(lg1 - lg2)
			den := math.Sqrt(df*math.Pi) * math.Pow(1+x*x/df, 0.5*(df+1))
			return true, ScalarReal(num / den), nil
		case "tcdf":
			return true, ScalarReal(studentTCDF(x, df)), nil
		case "tinv":
			if x <= 0 || x >= 1 {
				return true, Value{}, fmt.Errorf("tinv requires 0 < p < 1")
			}
			return true, ScalarReal(invertMonotonicCDF(func(t float64) float64 { return studentTCDF(t, df) }, x, -10, 10)), nil
		}

	case "chipdf", "chicdf", "chiinv":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(x, df) requires 2 real arguments", fn)
		}
		x, err1 := args[0].Real()
		df, err2 := args[1].Real()
		if err1 != nil || err2 != nil || df <= 0 {
			return true, Value{}, fmt.Errorf("%s requires real x and df > 0", fn)
		}
		switch fn {
		case "chipdf":
			if x <= 0 {
				return true, ScalarReal(0), nil
			}
			k := 0.5 * df
			lg, _ := math.Lgamma(k)
			return true, ScalarReal(math.Exp((k-1)*math.Log(x) - 0.5*x - k*math.Ln2 - lg)), nil
		case "chicdf":
			if x <= 0 {
				return true, ScalarReal(0), nil
			}
			return true, ScalarReal(regGammaP(0.5*df, 0.5*x)), nil
		case "chiinv":
			if x <= 0 || x >= 1 {
				return true, Value{}, fmt.Errorf("chiinv requires 0 < p < 1")
			}
			return true, ScalarReal(invertMonotonicCDF(func(v float64) float64 {
				if v <= 0 {
					return 0
				}
				return regGammaP(0.5*df, 0.5*v)
			}, x, 0, df+10*math.Sqrt(2*df))), nil
		}

	case "fpdf", "fcdf":
		if len(args) != 3 {
			return true, Value{}, fmt.Errorf("%s(x, d1, d2) requires 3 real arguments", fn)
		}
		x, err1 := args[0].Real()
		d1, err2 := args[1].Real()
		d2, err3 := args[2].Real()
		if err1 != nil || err2 != nil || err3 != nil || d1 <= 0 || d2 <= 0 {
			return true, Value{}, fmt.Errorf("%s requires real x and positive d1, d2", fn)
		}
		if x <= 0 {
			return true, ScalarReal(0), nil
		}
		if fn == "fcdf" {
			return true, ScalarReal(regBetaInc(0.5*d1, 0.5*d2, (d1*x)/(d1*x+d2))), nil
		}
		lg1, _ := math.Lgamma(0.5 * d1)
		lg2, _ := math.Lgamma(0.5 * d2)
		lg12, _ := math.Lgamma(0.5 * (d1 + d2))
		logB := lg1 + lg2 - lg12
		logPdf := 0.5*d1*math.Log(d1) + 0.5*d2*math.Log(d2) + (0.5*d1-1)*math.Log(x) -
			0.5*(d1+d2)*math.Log(d1*x+d2) - logB
		return true, ScalarReal(math.Exp(logPdf)), nil

	case "binompmf", "binomcdf":
		if len(args) != 3 {
			return true, Value{}, fmt.Errorf("%s(k, n, p) requires 3 arguments", fn)
		}
		kF, err1 := args[0].Real()
		nF, err2 := args[1].Real()
		p, err3 := args[2].Real()
		if err1 != nil || err2 != nil || err3 != nil || p < 0 || p > 1 || nF < 0 {
			return true, Value{}, fmt.Errorf("invalid binomial parameters")
		}
		k := int(math.Floor(kF))
		n := int(math.Round(nF))
		pmf := func(i int) float64 {
			if i < 0 || i > n {
				return 0
			}
			return nCrExact(float64(n), float64(i)) * math.Pow(p, float64(i)) * math.Pow(1-p, float64(n-i))
		}
		if fn == "binompmf" {
			return true, ScalarReal(pmf(int(math.Round(kF)))), nil
		}
		if k < 0 {
			return true, ScalarReal(0), nil
		}
		if k >= n {
			return true, ScalarReal(1), nil
		}
		var sum float64
		for i := 0; i <= k; i++ {
			sum += pmf(i)
		}
		return true, ScalarReal(sum), nil

	case "poisspmf", "poisscdf":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(k, lambda) requires 2 arguments", fn)
		}
		kF, err1 := args[0].Real()
		lam, err2 := args[1].Real()
		if err1 != nil || err2 != nil || lam <= 0 {
			return true, Value{}, fmt.Errorf("lambda must be positive")
		}
		pmf := func(i int) float64 {
			if i < 0 {
				return 0
			}
			lg, _ := math.Lgamma(float64(i + 1))
			return math.Exp(float64(i)*math.Log(lam) - lam - lg)
		}
		if fn == "poisspmf" {
			return true, ScalarReal(pmf(int(math.Round(kF)))), nil
		}
		k := int(math.Floor(kF))
		if k < 0 {
			return true, ScalarReal(0), nil
		}
		var sum float64
		for i := 0; i <= k; i++ {
			sum += pmf(i)
		}
		return true, ScalarReal(sum), nil

	case "geompmf", "geomcdf":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(k, p) requires 2 arguments", fn)
		}
		kF, err1 := args[0].Real()
		p, err2 := args[1].Real()
		if err1 != nil || err2 != nil || p <= 0 || p > 1 {
			return true, Value{}, fmt.Errorf("p must be in (0, 1]")
		}
		k := math.Floor(kF)
		if k < 1 {
			return true, ScalarReal(0), nil
		}
		if fn == "geompmf" {
			return true, ScalarReal(math.Pow(1-p, k-1) * p), nil
		}
		return true, ScalarReal(1 - math.Pow(1-p, k)), nil

	case "exppdf", "expcdf":
		if len(args) != 2 {
			return true, Value{}, fmt.Errorf("%s(x, lambda) requires 2 arguments", fn)
		}
		x, err1 := args[0].Real()
		lam, err2 := args[1].Real()
		if err1 != nil || err2 != nil || lam <= 0 {
			return true, Value{}, fmt.Errorf("lambda must be positive")
		}
		if x < 0 {
			return true, ScalarReal(0), nil
		}
		if fn == "exppdf" {
			return true, ScalarReal(lam * math.Exp(-lam*x)), nil
		}
		return true, ScalarReal(1 - math.Exp(-lam*x)), nil
	}

	return false, Value{}, nil
}

