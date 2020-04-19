const sorts = require('./sorts');
const { runTest } = require('./tests');

runTest(sorts.bubbleSort)
runTest(sorts.selectionSort)
runTest(sorts.insertionSort)
runTest(sorts.mergeSort)
runTest(sorts.quickSort)
runTest(sorts.heapSort)
runTest(sorts.countingSort)
runTest(sorts.bucketSort)
runTest(sorts.radixSort)
runTest(sorts.shellSort)
runTest(sorts.combSort)
runTest(sorts.bogoSort, 10)

