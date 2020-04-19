// Bubble Sort
function bubbleSort(array) {
    for (let i = array.length - 1; i > 0; --i) {
        for (let j = 0; j < i; ++j) {
            if (array[i] < array[j]) {
                swapArrayElements(array, i, j)
            }
        }
    }
}

// Selection Sort
function selectionSort(array) {

    for (let i = 0; i < array.length - 1; ++i) {
        let indexOfMin = i;

        for (let j = i + 1; j < array.length; ++j) {
            if (array[j] < array[indexOfMin]) {
                indexOfMin = j;
            }
        }

        if (indexOfMin !== i) {
            swapArrayElements(array, indexOfMin, i)
        }
    }
}

// Insertion sort
function insertionSort(array) {

    for (let i = 1; i < array.length; ++i) {
        let j = i;
        let tmp = array[j];

        while (j > 0 && array[j - 1] > tmp) {
            array[j] = array[j - 1];
            --j;
        }

        array[j] = tmp;
    }
}

// Merge sort
function mergeSort(array) {
    let result = mergeSortRec(array);

    for (let i = 0; i < array.length; ++i) {
        array[i] = result[i];
    }
}

function mergeSortRec(array) {

    if (array.length === 1) {
        return array;
    }

    let mid = Math.floor(array.length / 2);
    let left = array.slice(0, mid);
    let right = array.slice(mid, array.length);

    return merge(mergeSortRec(left), mergeSortRec(right));
}

function merge(left, right) {
    let result = [];
    let indexLeft = 0;
    let indexRight = 0;

    while (indexLeft < left.length && indexRight < right.length) {
        if (left[indexLeft] < right[indexRight]) {
            result.push(left[indexLeft]);
            ++indexLeft;
        } else {
            result.push(right[indexRight]);
            ++indexRight;
        }
    }

    while (indexLeft < left.length) {
        result.push(left[indexLeft]);
        ++indexLeft;
    }

    while (indexRight < right.length) {
        result.push(right[indexRight]);
        ++indexRight;
    }

    return result;
}

// Quick Sort
function quickSort(array) {
    quick(array, 0, array.length - 1);
}

function quick(array, left, right) {

    if (left < right) {
        let index = partition(array, left, right);

        quick(array, left, index - 1);
        quick(array, index, right);
    }
}

function partition(array, left, right) {
    let pivot = array[Math.floor((left + right) / 2)];

    while (left <= right) {

        while (array[left] < pivot) {
            ++left;
        }

        while (array[right] > pivot) {
            --right;
        }

        if (left <= right) {
            let tmp = array[left];
            array[left] = array[right];
            array[right] = tmp;

            ++left;
            --right;
        }
    }

    return left;
}

// Heap Sort
function heapSort(array) {
    let heapSize = array.length;

    buildHeap(array);

    while (heapSize > 1) {
        --heapSize;

        swapArrayElements(array, 0, heapSize);

        heapify(array, heapSize, 0);
    }
}

function buildHeap(array) {

    for (let i = Math.floor(array.length / 2); i >= 0; --i) {
        heapify(array, array.length, i);
    }
}

function heapify(array, heapSize, i) {
    let left = 2 * i + 1;
    let right = 2 * i + 2;
    let largest = i;

    if (left < heapSize && array[left] > array[largest]) {
        largest = left;
    }

    if (right < heapSize && array[right] > array[largest]) {
        largest = right;
    }

    if (i !== largest) {
        swapArrayElements(array, i, largest);

        heapify(array, heapSize, largest);
    }
}

// Counting Sort
function countingSort(array) {
    let maxValue = getMaxElement(array);
    let counts = new Array(maxValue + 1);
    let sortedIndex = 0;

    for (let i = 0; i < array.length; ++i) {
        if (counts[array[i]] === undefined) {
            counts[array[i]] = 0;
        }

        ++counts[array[i]];
    }

    for (let i = 0; i < counts.length; ++i) {
        while (counts[i] > 0) {
            array[sortedIndex] = i;
            ++sortedIndex;
            --counts[i];
        }
    }
}

// Bucket Sort
function bucketSort(array) {
    let bucketSize = 5;
    let minVal = getMinElement(array);
    let maxVal = getMaxElement(array);
    let bucketCount = Math.floor((maxVal - minVal) / bucketSize) + 1;

    let buckets = new Array(bucketCount);

    for (let i = 0; i < buckets.length; ++i) {
        buckets[i] = [];
    }

    for (let i = 0; i < array.length; ++i) {
        buckets[Math.floor((array[i] - minVal) / bucketCount)].push(array[i]);
    }

    let result = [];

    for (let i = 0; i < buckets.length; ++i) {
        insertionSort(buckets[i]);

        for (let j = 0; j < buckets[i].length; ++j) {
            result.push(buckets[i][j]);
        }
    }

    for (let i = 0; i < array.length; ++i) {
        array[i] = result[i];
    }
}

// Radix Sort
function radixSort(array) {
    let minVal = getMinElement(array);
    let maxVal = getMaxElement(array);

    let radixBase = 10;
    let significantDigit = 1;

    while ((maxVal - minVal) / significantDigit >= 2) {
        countingSortForRadix(array, radixBase, significantDigit, minVal);
        significantDigit *= radixBase;
    }
}

function countingSortForRadix(array, radixBase, significantDigit, minVal) {
    let counts = new Array(radixBase);
    let tmpArray = new Array(array.length);

    let countsIndex;

    for (let i = 0; i < counts.length; ++i) {
        counts[i] = 0;
    }

    for (let i = 0; i < array.length; ++i) {
        countsIndex = Math.floor(((array[i] - minVal) / significantDigit) % radixBase);
        ++counts[countsIndex];
    }

    for (let i = 1; i < counts.length; ++i) {
        counts[i] += counts[i - 1];
    }

    for (let i = array.length - 1; i >= 0; --i) {
        countsIndex = Math.floor(((array[i] - minVal) / significantDigit) % radixBase);
        tmpArray[--counts[countsIndex]] = array[i];
    }

    for (let i = 0; i < array.length; ++i) {
        array[i] = tmpArray[i];
    }
}

// Shell Sort
function shellSort(array) {

    for (let gap = Math.floor(array.length / 2); gap > 0; gap = Math.floor(gap / 2)) {
        for (let i = gap; i < array.length; ++i) {
            let tmp = array[i];

            let j;
            for (j = i; j >= gap && array[j - gap] > tmp; j -= gap) {
                array[j] = array[j - gap];
            }

            array[j] = tmp;
        }
    }
}

// Comb Sort
function combSort(array) {
    let gap = array.length;
    let swapped = true;

    while (gap !== 1 || swapped) {
        gap = getNextGap(gap);

        swapped = false;

        for (let i = 0; i < array.length - gap; ++i) {
            if (array[i] > array[i + gap]) {
                swapArrayElements(array, i, i + gap);
                swapped = true;
            }
        }
    }
}

function getNextGap(gap) {
    const shrinkFacor = 1.247;

    gap = Math.floor(gap / shrinkFacor);

    return gap < 1 ? 1 : gap;
}

// Bogo Sort
function bogoSort(array) {
    let sorted = false;

    while (sorted === false) {
        shuffle(array);
        sorted = isSorted(array);
    }
}

function shuffle(array) {
    // for (let i = array.length, j, tmp; i > 0; --i) {
    for (let i = 1, j, tmp; i < array.length + 1; ++i) {
        j = Math.floor(Math.random() * i);
        tmp = array[i - 1];
        array[i - 1] = array[j];
        array[j] = tmp;
    }
}

function isSorted(array) {

    for (let i = 0; i < array.length - 1; ++i) {
        if (array[i] > array[i + 1]) {
            return false;
        }
    }

    return true;
}

// utils
function swapArrayElements(array, index1, index2) {
    console.assert(index1 < array.length && index2 < array.length);

    let tmp = array[index1];
    array[index1] = array[index2];
    array[index2] = tmp;
}

function getMinElement(array) {
    let min = array[0];

    for (let i = 0; i < array.length; ++i) {
        if (min > array[i]) {
            min = array[i];
        }
    }

    return min;
}

function getMaxElement(array) {
    let max = array[0];

    for (let i = 0; i < array.length; ++i) {
        if (max < array[i]) {
            max = array[i];
        }
    }

    return max;
}

module.exports = {
    'bubbleSort': bubbleSort,
    'selectionSort': selectionSort,
    'insertionSort': insertionSort,
    'mergeSort': mergeSort,
    'quickSort': quickSort,
    'heapSort': heapSort,
    'countingSort': countingSort,
    'bucketSort': bucketSort,
    'radixSort': radixSort,
    'shellSort': shellSort,
    'combSort': combSort,
    'bogoSort': bogoSort
};