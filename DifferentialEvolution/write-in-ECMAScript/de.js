function de(F, CR, D, fnum) {          //Amplifactor Rate, Crossover Rate, Dimension


    var d1 = new Date().getTime();          //startTime

    var Np = D * 10,            //Population Size
        maximum = 100,          //Low Bound
        minimum = -100,         //Height Bound
        nfes = 0,               //number of function evaluates
        Max_nfes = Np * 1000,   //Max of nfes
        stop_con = false,       //Stop condition
        arr = [],           //for array translation , amount to all individuals of all population
        pop = [],               //Initialized population
        cpop,                   //Population after Crossover operation
        value,                  //an array to record the mimimum value and its index in a population
        bestIndex,              //index of minimum value in population
        result = [];            //record the result

    //initialize population
    for(var i = 0; i < Np; i++){
        for(var l =0; l < D; l++){
            arr.push(rand(minimum, maximum));
        }
        pop.push(arr);
        arr = [];
    }

    cpop = pop;

    var fitness = ackley(pop);
    nfes = nfes + Np;

    // bestIndex = getMinimum(fitness)[0];

    if (nfes > Max_nfes) stop_con = true;

    while (!stop_con) {

        //mutate + crossover
        for (var j = 0; j < Np; j++) {
            var a, b, c, randj;
            do { a = Math.floor(rand(0, Np)) } while (a == j);
            do { b = Math.floor(rand(0, Np)) } while (b == j || b == a);
            do { c = Math.floor(rand(0, Np)) } while (c == j || c == a || c == b);
            randj = Math.ceil(rand(0, D));
            for (var k = 0; k < D; k++) {
                if (k == randj || Math.random() <= CR) {
                    cpop[j][k] = pop[a][k] + F * (pop[b][k] - pop[c][k]);
                    // cpop[j][k] = pop[a][k] + F * (pop[bestIndex][k] - pop[c][k]);
                    while (cpop[j][k] > maximum || cpop[j][k] < minimum) {
                        cpop[j][k] = rand(minimum, maximum);
                    }
                }
            }
        }

        //select
        var cfitness = ackley(cpop);
        nfes = nfes + Np;
        for (var n = 0; n < fitness.length; n++) {
            if (cfitness[n] < fitness[n]) {
                fitness[n] = cfitness[n];
                pop[n] = cpop[n];
            }
        }

        value = getMinimum(fitness);
        bestIndex = value[0];
        cpop = pop;

        if (nfes > Max_nfes) stop_con = true;

        result.push(value[1]);
    }

    var d2 = new Date().getTime();
    var timeCount = (d2 - d1) / 1000;

    //time record, value collection
    return [timeCount, result];
}

function getMinimum(arr) {
    var index = 0;
    var bestIndex = 0;
    var bestValue = arr[index];
    while (++index < arr.length) {
        if (arr[index] < bestValue) {
            bestValue = arr[index];
            bestIndex = index;
        }
    }
    return [bestIndex, bestValue];
}

function rand(min, max) {
    return min + Math.random() * (max - min);
}

function sphere(matrix) {
    var sum = 0;
    var arr = [];
    for (var i = 0; i < matrix.length; i++) {
        for (var j = 0; j < matrix[i].length; j++) {
            sum += (matrix[i][j] * matrix[i][j]);
        }
        arr.push(sum);
        sum = 0;
    }
    return arr;
}

function ackley(matrix) {
    var arr = [];
    for (let i = 0; i < matrix.length; i++) {
        var sum1 = 0,
            sum2 = 0,
            f = 0;
        for (let j = 0; j < matrix[i].length; j++) {
            sum1 += matrix[i][j] * matrix[i][j];
            sum2 += Math.cos(2 * Math.PI * matrix[i][j]);
        }
        sum1 = -0.2 * Math.sqrt(sum1 / matrix[i].length);
        sum2 /= matrix[i].length;
        f = Math.E - 20 * Math.exp(sum1) - Math.exp(sum2) + 20;
        arr.push(f);
    }
    return arr;
}

function rastrigin(matrix) {
    var arr = [];
    for (let i = 0; i < matrix.length; i++) {
        var sum = 0;
        for (let j = 0; j < matrix[i].length; j++) {
            sum += (matrix[i][j] * matrix[i][j] - 10 * Math.cos(2 * Math.PI * matrix[i][j]) + 10);
        }
        arr.push(sum);
    }
    return arr;
}

module.exports = de;