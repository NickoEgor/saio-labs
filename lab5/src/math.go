package src

import (
	"errors"
	"math"

	"gonum.org/v1/gonum/mat"
)

type intPair struct {
	first  int
	second int
}

func fillDense(dense *mat.Dense, value float64) {
	rows, cols := dense.Dims()
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			dense.Set(i, j, value)
		}
	}
}

func minIndex(dense *mat.Dense) intPair {
	rows, cols := dense.Dims()
	minVal := mat.Min(dense)
	minInd := intPair{-1, -1}
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			if dense.At(i, j) == minVal {
				minInd = intPair{i, j}
				break
			}

		}
	}
	return minInd
}

func minIndexVec(dense *mat.VecDense) int {
	length := dense.Len()
	minVal := mat.Min(dense)
	minInd := -1
	for i := 0; i < length; i++ {
		if dense.AtVec(i) == minVal {
			minInd = i
			break
		}
	}
	return minInd
}

// TransportPotentials - solves transport problem with potentials method
func TransportPotentials(cost *mat.Dense, need, stock *mat.VecDense) (*mat.Dense, error) {
	_, cols0 := cost.Dims()

	amountStock := mat.Sum(stock)
	amountNeed := mat.Sum(need)

	if amountStock < amountNeed {
		return nil, errors.New("No solutions exists")
	}

	costs := cost.Grow(0, 1)

	needs := mat.NewVecDense(cols0+1, nil)
	for i := 0; i < cols0; i++ {
		needs.SetVec(i, need.AtVec(i))
	}
	needs.SetVec(cols0, amountStock-amountNeed)

	stock0 := stock.AtVec(0)
	need0 := needs.AtVec(0)
	rows, cols := costs.Dims()
	bFlow := mat.NewDense(rows, cols, nil)
	mFlow := mat.NewDense(rows, cols, nil)

	for i, j := 0, 0; i < rows && j < cols; {
		value := math.Min(stock0, need0)
		mFlow.Set(i, j, value)
		bFlow.Set(i, j, 1)

		stock0 -= value
		need0 -= value

		if stock0 == 0.0 && i != rows-1 {
			i++
			stock0 = stock.AtVec(i)
		} else if need0 == 0.0 && j != cols-1 {
			j++
			need0 = needs.AtVec(j)
		} else {
			i++
			j++
			if i != rows || j != cols {
				panic("i != rows or j != cols")
			}
		}
	}

	isRunning := true
	for isRunning {
		isRunning, mFlow, bFlow = transportStep(&costs, mFlow, bFlow)
	}

	return mat.DenseCopyOf(mFlow.Slice(0, rows, 0, cols-1)), nil
}

func transportStep(costs *mat.Matrix, mFlow, bFlow *mat.Dense) (bool, *mat.Dense, *mat.Dense) {
	rows, cols := (*costs).Dims()

	// fmt.Println("bFlow")
	// io.MatPrint(bFlow)
	potentials := countPotentials(costs, bFlow)
	// fmt.Println("paul pot")
	// io.MatPrint(potentials)

	nFlow := mat.NewDense(rows, cols, nil)
	fillDense(nFlow, 1.0)
	nFlow.Sub(nFlow, bFlow)

	tempMat := mat.NewDense(rows, cols, nil)
	tempMat.MulElem(potentials, nFlow)

	minInd := minIndex(tempMat)
	minVal := mat.Min(tempMat)

	if minVal >= 0 {
		return false, mFlow, bFlow
	}

	// fmt.Println("min ind", minInd)
	// io.MatPrint(bFlow)

	historyPath := make([]intPair, 0)
	path := findPath(minInd, 0, bFlow, &historyPath, minInd)

	// fmt.Println("path len", len(path))

	cycleValuesMinus := mat.NewVecDense(len(path), nil)
	for i := 0; i < len(path); i++ {
		cycleValuesMinus.SetVec(i, mFlow.At(path[i].first, path[i].second))
	}

	minIndeces := make([]float64, 0)
	for i := 1; i < cycleValuesMinus.Len(); i += 2 {
		minIndeces = append(minIndeces, cycleValuesMinus.AtVec(i))
	}
	minDense := mat.NewVecDense(len(minIndeces), minIndeces)
	minIndexID := minIndexVec(minDense)*2 + 1
	swapIndex := path[minIndexID]
	minFlow := cycleValuesMinus.AtVec(minIndexID)

	for i := 0; i < len(path); i += 2 {
		ind := path[i]
		mFlow.Set(ind.first, ind.second,
			mFlow.At(ind.first, ind.second)+minFlow)
	}

	for i := 1; i < len(path); i += 2 {
		ind := path[i]
		mFlow.Set(ind.first, ind.second,
			mFlow.At(ind.first, ind.second)-minFlow)
	}

	// fmt.Println(swapIndex)
	// fmt.Println(minInd)
	bFlow.Set(swapIndex.first, swapIndex.second, 0)
	bFlow.Set(minInd.first, minInd.second, 1)

	return true, mFlow, bFlow
}

func findPath(curIndex intPair, direction int, bFlow *mat.Dense, historyPath *[]intPair, startIndex intPair) []intPair {
	rows, cols := bFlow.Dims()
	x, y := curIndex.first, curIndex.second

	var parent intPair
	if len(*historyPath) != 0 {
		parent = (*historyPath)[len(*historyPath)-1]
	} else {
		parent = intPair{-1, -1}
	}

	if curIndex == startIndex && len(*historyPath) != 0 {
		path := make([]intPair, len(*historyPath))
		copy(path, *historyPath)
		return path
	}

	(*historyPath) = append((*historyPath), curIndex)

	var result []intPair
	if direction == 0 {
		for i := 0; i < rows; i++ {
			nextIndex := intPair{i, y}
			if nextIndex != parent && nextIndex != curIndex &&
				(bFlow.At(nextIndex.first, nextIndex.second) == 1 ||
					nextIndex == startIndex) {
				result = findPath(nextIndex, 1, bFlow, historyPath, startIndex)
			}

			if len(result) > 0 {
				break
			}
		}
	} else {
		for i := 0; i < cols; i++ {
			nextIndex := intPair{x, i}
			if nextIndex != parent && nextIndex != curIndex &&
				(bFlow.At(nextIndex.first, nextIndex.second) == 1 ||
					nextIndex == startIndex) {
				result = findPath(nextIndex, 0, bFlow, historyPath, startIndex)
			}

			if len(result) > 0 {
				break
			}
		}
	}

	(*historyPath) = (*historyPath)[:len(*historyPath)-1]

	return result
}

func countPotentials(costs *mat.Matrix, bFlow *mat.Dense) *mat.Dense {
	rows, cols := (*costs).Dims()
	u := mat.NewVecDense(rows, nil)
	v := mat.NewVecDense(cols, nil)
	used := mat.NewDense(rows, cols, nil)

	queue := make([]intPair, 1)
	queue[0] = intPair{0, 0}

	for len(queue) != 0 {
		pair := queue[0]
		queue = queue[1:]
		corr, direction := pair.first, pair.second
		if direction == 0 {
			for i := 0; i < cols; i++ {
				if bFlow.At(corr, i) == 1 && used.At(corr, i) == 0 {
					used.Set(corr, i, 1)
					v.SetVec(i, (*costs).At(corr, i)-u.AtVec(corr))
					queue = append(queue, intPair{i, 1})
				}
			}
		} else {
			for i := 0; i < rows; i++ {
				if bFlow.At(i, corr) == 1 && used.At(i, corr) == 0 {
					used.Set(i, corr, 1)
					u.SetVec(i, (*costs).At(i, corr)-v.AtVec(corr))
					queue = append(queue, intPair{i, 0})
				}
			}
		}
	}

	result := mat.NewDense(rows, cols, nil)

	result.MulElem(*costs, bFlow)
	result.Add(result, *costs)

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			value := result.At(i, j)
			value -= u.AtVec(i)
			value -= v.AtVec(j)
			result.Set(i, j, value)
		}
	}

	return result
}
