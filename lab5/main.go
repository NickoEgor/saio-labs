package main

import (
	"fmt"
	"os"

	"gonum.org/v1/gonum/mat"

	"lab5/src"
)

func main() {
	var filename string
	if len(os.Args) < 2 {
		filename = "examples/ex0.txt"
	} else {
		filename = os.Args[1]
	}

	Cost, Need, Stock := src.EnterTransportConditions(filename)

	flow, err := src.TransportPotentials(Cost, Need, Stock)
	if err != nil {
		fmt.Println(err)
		return
	}

	fmt.Println("Flow matrix:")
	src.MatPrint(flow)

	rows, cols := flow.Dims()

	fmt.Println("Consumption:")
	consumption := mat.NewVecDense(rows, nil)
	for i := 0; i < rows; i++ {
		consumption.SetVec(i, mat.Sum(flow.RowView(i)))
	}
	src.MatPrint(consumption.T())

	fmt.Println("Satisfied needs:")
	satisfaction := mat.NewVecDense(cols, nil)
	for i := 0; i < cols; i++ {
		satisfaction.SetVec(i, mat.Sum(flow.ColView(i)))
	}
	src.MatPrint(satisfaction.T())

	fmt.Println("Minimum cost:")
	tempSum := mat.NewDense(rows, cols, nil)
	tempSum.MulElem(Cost, flow)
	cost := mat.Sum(tempSum)
	fmt.Println(cost)
}
