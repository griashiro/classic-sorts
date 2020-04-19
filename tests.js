function genRandArray(maxArrayLength = 100, maxNumber = 1000) {
    let array = [];

    function getRandomNumber(min, max) {
        return Math.random() * (max - min) ^ 0;
    }

    for (let i = 0; i < maxArrayLength; ++i) {
        array.push(getRandomNumber(0, maxNumber));
    }

    return array;
}

function runTest(sortingFunction, arrayLen) {
    let array = genRandArray(arrayLen);
    let standard = array.slice();
    standard.sort((a, b) => a - b);

    sortingFunction(array);

    for (let i = 0; i < array.length; ++i) {
        if (standard[i] !== array[i]) {
            console.log('fault', sortingFunction.name);
            return;
        }
    }

    console.log('pass', sortingFunction.name);
}

module.exports = {
    runTest,
    genRandArray,
};
