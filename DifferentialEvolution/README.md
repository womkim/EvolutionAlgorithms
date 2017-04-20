## Evolution Algorithm: Differential Evolution

---

Differential Evolution (DE) is a method that optimizes a problem by iteratively trying to improve a candidate solution with regard to a given measure of quality. Such methods are commonly known as metaheuristics as they make few or no assumptions about the problem being optimized and can search very large spaces of candidate solutions. However, metaheuristics such as DE do not guarantee an optimal solution is ever found.

this is a simple code write in java, more advances need to complete...

(add some code with ecmascript, to test envorinment in node)

### How to run

---

#### Run ClassicalDE20160323.jar file

```markdown
- Make sure the folder "input_data" and .jar file are in the same directory.
- using Java environment to run the .jar file with command line.
- for example: 
	Run command: java -jar 0.5 0.5 1 1 30 10 30 5 results
          the first parameter is the scaling factor, value of F;
          the second parameter is the crossover rate, value of CR;
          the third parameter is DE mutation strategy, 
                  there are 7 strategies: 
                              < DE/rand/1 > , 
                              < DE/best/1 > , 
                              < DE/rand/2 > , 
                              < DE/best/2 > ,
                              < DE/rand-To-best/1 > , 
                              < DE/current-To-best/1 > ,
                              < DE/rand-To-best/2 >;
          the fourth parameter is the start function number, value of 1 to 30;
          the fifth parameter is the end function number, value of 1 to 30, and must be bigger than start function number or be equal to it;
          the following parameter(s) is(are) the dimension of the problem (must be set as 10,30,50,100 while function number is within 1 to 15);
          the penultimate parameter is the compute times under current conditions, at least 1 time and should not be more than 51 times;
          the last parameter is the file directory for saving results, see details in the directory.
```

#### Run code in IDE

```markdown
- copy the folder "input_data" to the project directory.
- set the parameters in program arguments of IDE configurations. with same as above.
- run Main.java
```
- About function numbers
- *function number 1 to 15 is cec15's function, more details see (Liang J J, Qu B Y, Suganthan P N, et al. Problem definitions and evaluation criteria for the CEC 2015 competition on learning-based real-parameter single objective optimization[J]. Technical Report201411A, Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou China and Technical Report, Nanyang Technological University, Singapore, 2014.).
- *And the other function numbers indicate:
  - 16 Sphere Function
  - 17 High Conditioned Elliptic Function
  - 18 Cigar Function
  - 19 Discus Function
  - 20 Rosenbrock's Function
  - 21 Ackley's Function
  - 22 Weierstrass Function
  - 23 Griewank's Function
  - 24 Rastrigin's Function
  - 25 Modified Schwefel's Function
  - 26 Katsuura Function
  - 27 HappyCat Function
  - 28 HGBat Function
  - 29 Expanded Griewank's plus Rosenbrock's function
  - 30 Expanded Scaffer's F6 Function

#### Queries

Have any concerns, please contact us on 337668563@qq.com