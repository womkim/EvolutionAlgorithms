/**
 * Created by Administrator on 2016/10/14.
 */
!function () {
    'use strict';
    var fs = require('fs');
    var de = require('./de');

    var str = '';
    var resultArr = [];
    for(var i = 0; i < 51; i++){
        var result = de(0.5, 0.3, 2);      //params: F, CR, strategy ; return [time, results]
        resultArr.push(result[1][result.length - 1]);
        console.log(i+1 + " The total time is : " + result[0] + "s");
        console.log(i+1 + " The best value is : " + result[1][result.length - 1]);
        fs.writeFile('./value/v_'+i+'.txt', result[1].join('\n'));
        fs.writeFile('result.txt', result[1][result.length - 1] + '\n', {flag : 'a'});
    }
    var resultMin = Math.min.apply(Math, resultArr);
    console.log('\n The Minimum result is : ' + resultMin);

}();