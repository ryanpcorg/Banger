package main

import (
	"go/constant"
	"go/token"
	"go/types"
)

func EvaluateExpression(expr string) float64 {

	fset := token.NewFileSet()

	tv, err := types.Eval(fset, nil, token.NoPos, expr)
	if err != nil {
		return 0
	}

	if tv.Value != nil {
		val, _ := constant.Float64Val(tv.Value)
		return val
	}

	return 0
}
