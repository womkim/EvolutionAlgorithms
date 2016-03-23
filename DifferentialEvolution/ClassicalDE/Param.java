
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class Param {
	
	final static double E  = 2.7182818284590452353602874713526625;
	final static double PI = 3.1415926535897932384626433832795029;
	
	static double Xmin = -100.0;
	static double Xmax = 100.0;
	
	static Random random = new Random();

	/**目标函数序号
	 * （序号1-15分别是CEC15的函数函数集序号）
	 * （序号16-30分别是	16	Sphere Function
	 * 					17	High Conditioned Elliptic Function
	 * 					18	Cigar Function
	 * 					19	Discus Function
	 * 					20	Rosenbrock's Function
	 * 					21	Ackley's Function
	 * 					22	Weierstrass Function
	 * 					23	Griewank's Function
	 * 					24	Rastrigin's Function
	 * 					25	Modified Schwefel's Function
	 * 					26	Katsuura Function
	 * 					27	HappyCat Function
	 * 					28	HGBat Function
	 * 					29	Expanded Griewank's plus Rosenbrock's function
	 * 					30	Expanded Scaffer's F6 Function
	 * ）
	 */
	
	/**-----------------------------------------------------------------------**/
	/**--------------------------------测试函数--------------------------------**/
	/**--------------------具体函数说明请参考CEC2015相关文件---------------------**/
	/**----------------------------------------------------------------------**/
	/**
	 * 目标函数（Sphere函数）
	 * 表达式：y = sum(v.^2)
	 * @param v		每一个个体
	 * @return		每个个体求出的解
	 */
	public static double sphere(double[] v){
		double sum = 0;
		for(int i=0; i<v.length; i++){
			sum += v[i]*v[i];
		}
		return sum;
	}	
	/**
	 * Elliptic函数
	 * @param v
	 * @return
	 */
	public static double ellips(double[] v){
		double f = 0;
		for (int i=0; i<v.length; i++)
	    {
	        f += Math.pow(10.0,6.0*i/(v.length-1))*v[i]*v[i];
	    }
	    return f;
	}
	/**
	 * Cigar函数
	 * @param v
	 * @return
	 */
	public static double cigar(double[] v){
		double sum = 0;
		sum = v[0]*v[0];
		for (int i=1; i<v.length; i++)
	    {
	        sum += Math.pow(10.0,6.0)*v[i]*v[i];
	    }
	    return sum;
	}
	/**
	 * Discus函数
	 * @param v
	 * @return
	 */
	public static double discus(double[] v){
		double f = 0;
		f = Math.pow(10.0,6.0)*v[0]*v[0];
	    for (int i=1; i<v.length; i++)
	    {
	        f += v[i]*v[i];
	    }
	    return f;
	}
	/**
	 * 目标函数（Rosenbrock Valley函数）<<此函数与matlab精度不同，算出的结果会有些小的偏差>>
	 * 表达式：y = sum( 100*{ [ v(xi+1) - v(xi).^2 ].^2 } + [ 1 - v(xi) ].^2)
	 * @param v
	 * @return
	 */
	public static double rosenbrock(double[] v){
		double f = 0;
		double tmp1,tmp2;
		for (int i=0; i<v.length-1; i++)
	    {
	    	tmp1=v[i]*v[i]-v[i+1];
			tmp2=v[i]-1.0;
	        f += 100.0*tmp1*tmp1 +tmp2*tmp2;
	    }
	    return f;
	}
	/**
	 * 目标函数（Ackley函数）
	 * 表达式： y = -20*exp(-0.2*sqrt( (1/D) * sum(v.^2) )) - exp( (1/D) * sum(cos(2*pi*v(i))) ) + 20 + exp(1);
	 * @param v
	 * @return
	 */
	public static double ackley(double[] v){
		double sum1 = 0, sum2 = 0;
		int n = v.length;
		for (int i=0; i<n; i++)
		{
			sum1 += v[i]*v[i];
			sum2 += Math.cos(2.0*PI*v[i]);
		}
		sum1 = -0.2*Math.sqrt(sum1/n);
		sum2 /= n;
		double f =  E - 20.0*Math.exp(sum1) - Math.exp(sum2) +20.0;
		return f;
	}
	/**
	 * Weierstrass函数
	 * @param v
	 * @return
	 */
	public static double weierstrass(double[] v){
		int k_max = 20;
		double sum, sum2 = 0;
		double a = 0.5;
		double b = 3.0;
		double f = 0;
		for (int i=0; i<v.length; i++)
	    {
	        sum = 0.0;
			sum2 = 0.0;
	        for (int j=0; j<=k_max; j++)
	        {
	            sum += Math.pow(a,j)*Math.cos(2.0*PI*Math.pow(b,j)*(v[i]+0.5));
				sum2 += Math.pow(a,j)*Math.cos(2.0*PI*Math.pow(b,j)*0.5);
	        }
	        f += sum;
	    }
		f -= v.length*sum2;
		return f;
	}
	/**
	 * 目标函数（Griewank函数）
	 * 表达式：y = sum((v.^2)/4000) - prod(cos(v./sqrt(xi))) + 1
	 * @param v		个体
	 * @return		求出的解
	 */
	public static double griewank(double[] v){
		double v1 = 0, v2 = 1, f = 0;
		for (int i=0; i<v.length; i++)
	    {
	        v1 += v[i]*v[i];
	        v2 *= Math.cos(v[i]/Math.sqrt(1.0+i));
	    }
	    f = 1.0 + v1/4000.0 - v2;
	    return f;
	}
	/**
	 * 目标函数（Rastrigin函数）
	 * 表达式：y = sum(v.^2 - 10.*cos(2.*pi.*v) + 10)
	 * @param v		个体
	 * @return		求出的解
	 */
	public static double rastrigin(double[] v){
		double sum = 0;
		for (int i = 0; i < v.length; i++) {
			sum += (v[i]*v[i] - 10 * Math.cos(2 * PI * v[i]) + 10);
		}
		return sum;
	}
	/**
	 * Schwefel 函数
	 * @param v
	 * @return
	 */
	public static double schwefel(double[] v){
		double f = 0;
		double tmp;
		for (int i=0; i<v.length; i++)
		{
	    	v[i] += 4.209687462275036e+002;
	    	if (v[i]>500)
			{
				f-=(500.0-(v[i]%500))*Math.sin(Math.pow(500.0-(v[i]%500),0.5));
				tmp=(v[i]-500.0)/100;
				f+= tmp*tmp/v.length;
			}
			else if (v[i]<-500)
			{
				f-=(-500.0+(Math.abs(v[i])%500))*Math.sin(Math.pow(500.0-(Math.abs(v[i])%500),0.5));
				tmp=(v[i]+500.0)/100;
				f+= tmp*tmp/v.length;
			}
			else
				f-=v[i]*Math.sin(Math.pow(Math.abs(v[i]),0.5));
	    }
	    f=4.189828872724338e+002*v.length+f;
	    return f;
	}
	/**
	 * Katsuura 函数
	 * @param z
	 * @return
	 */
	public static double katsuura(double[] z){
		int i,j;
		double temp,tmp1,tmp2,tmp3;
		int nx = z.length;
		tmp3=Math.pow(1.0*nx,1.2);
	    double f=1.0;
	    for (i=0; i<nx; i++)
		{
			temp=0.0;
			for (j=1; j<=32; j++)
			{
				tmp1=Math.pow(2.0,j);
				tmp2=tmp1*z[i];
				temp += Math.abs(tmp2-Math.floor(tmp2+0.5))/tmp1;
			}
			f *= Math.pow(1.0+(i+1)*temp,10.0/tmp3);
	    }
		tmp1=10.0/nx/nx;
	    f=f*tmp1-tmp1;
	    return f;
	}
	/**
	 * HappyCat 函数
	 * @param z
	 * @return
	 */
	public static double happycat(double[] z){
		int i;
		double alpha,r2,sum_z;
		alpha = 1.0/8.0;
		r2 = 0.0;
		sum_z = 0.0;
		double f = 0.0;
		int nx = z.length;
		for (i=0;i<nx;i++)
		{
			r2 += z[i]*z[i];
			sum_z += z[i];
			
		}
		f = Math.pow(Math.abs(r2-nx), 2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
		return f;
	}
	/**
	 * HGBat 函数
	 * @param z
	 * @return
	 */
	public static double hgbat(double[] z){
		int i;
		double alpha,r2,sum_z;
		alpha=1.0/4.0;
		r2 = 0.0;
		sum_z=0.0;
		double f = 0;
		int nx = z.length;
	    for (i=0; i<nx; i++)
	    {
	        r2 += z[i]*z[i];
			sum_z += z[i];
	    }
	    f=Math.pow(Math.abs(Math.pow(r2,2.0)-Math.pow(sum_z,2.0)),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
	    return f;
	}
	/**
	 * Expanded Griewank's plus Rosenbrock's 函数
	 * @param z
	 * @return
	 */
	public static double grie_rosen(double[] z){
		int i;
	    double temp,tmp1,tmp2;
	    double f=0.0;
	    int nx = z.length;
	    for (i=0; i<nx-1; i++)
	    {
			tmp1 = z[i]*z[i]-z[i+1];
			tmp2 = z[i]-1.0;
	        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
	        f += (temp*temp)/4000.0 - Math.cos(temp) + 1.0;
	    }
		tmp1 = z[nx-1]*z[nx-1]-z[0];
		tmp2 = z[nx-1]-1.0;
	    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
	    f += (temp*temp)/4000.0 - Math.cos(temp) + 1.0 ;
	    return f;
	}
	/**
	 * Expanded Scaffer's F6 Function
	 * @param z
	 * @return
	 */
	public static double escaffer6(double[] z){
		int i;
	    double temp1, temp2;
	    double f = 0.0;
	    int nx = z.length;
	    for (i=0; i<nx-1; i++)
	    {
	        temp1 = Math.sin(Math.sqrt(z[i]*z[i]+z[i+1]*z[i+1]));
			temp1 =temp1*temp1;
	        temp2 = 1.0 + 0.001*(z[i]*z[i]+z[i+1]*z[i+1]);
	        f += 0.5 + (temp1-0.5)/(temp2*temp2);
	    }
	    temp1 = Math.sin(Math.sqrt(z[nx-1]*z[nx-1]+z[0]*z[0]));
		temp1 =temp1*temp1;
	    temp2 = 1.0 + 0.001*(z[nx-1]*z[nx-1]+z[0]*z[0]);
	    f += 0.5 + (temp1-0.5)/(temp2*temp2);
	    return f;
	}
	/**-----------------------------------------------------------------------**/
	/**
	 * 所有函数
	 * 求<B>单个个体</B>的目标函数值
	 * @param indiv			传入的个体
	 * @param D				个体的维度（CEC2015的函数维度只能是10，30，50，100；其他函数维度不限制）
	 * @param func_num		目标函数序号
	 * 					（序号1-15分别是CEC15的函数函数集序号）
	 * 					（序号16-30分别是	
	 * 									16	Sphere Function
	 * 									17	High Conditioned Elliptic Function
	 * 									18	Cigar Function
	 * 									19	Discus Function
	 * 									20	Rosenbrock's Function
	 * 									21	Ackley's Function
	 * 									22	Weierstrass Function
	 * 									23	Griewank's Function
	 * 									24	Rastrigin's Function
	 * 									25	Modified Schwefel's Function
	 * 									26	Katsuura Function
	 * 									27	HappyCat Function
	 * 									28	HGBat Function
	 * 									29	Expanded Griewank's plus Rosenbrock's function
	 * 									30	Expanded Scaffer's F6 Function
	 * ）
	 * @return				返回该个体的目标函数值
	 * @throws Exception
	 */
	public static double getFuncValue(double[] indiv, int D, int func_num) throws Exception{
		if(func_num > 15){
			switch(func_num){
				case 16 : return sphere(indiv); 
				case 17 : return ellips(indiv); 
				case 18 : return cigar(indiv); 
				case 19 : return discus(indiv);
				case 20 : return rosenbrock(indiv); 
				case 21 : return ackley(indiv);
				case 22 : return weierstrass(indiv);
				case 23 : return griewank(indiv);
				case 24 : return rastrigin(indiv);
				case 25 : return schwefel(indiv);
				case 26 : return katsuura(indiv);
				case 27 : return happycat(indiv);
				case 28 : return hgbat(indiv);
				case 29 : return grie_rosen(indiv);
				case 30 : return escaffer6(indiv);
				default : System.out.println("没有指定函数可用！"); return 0;
			}
		}else{
			double[] f = new double[1];
			testfunc tf = new testfunc();		
			tf.test_func(indiv,f,D,1,func_num);
			return f[0];
		}
	}
	/**
	 * CEC2015的所有函数（限定内容和上面一样）
	 * 求<B>种群中所有个体</B>的函数值
	 * @param pop		种群
	 * @param func_num	函数序号
	 * @return			返回种群内所有个体的目标函数值
	 * @throws Exception
	 */
	public static double[] getFuncAllValue(double[][] pop, int func_num) throws Exception{
		int m = pop.length;
		double[] f = new double[m];
		int n = pop[0].length;
		if(func_num>15){
			for (int i = 0; i < m; i++) {
				switch(func_num){
					case 16 : f[i] = sphere(pop[i]); break;
					case 17 : f[i] = ellips(pop[i]); break;
					case 18 : f[i] = cigar(pop[i]); break;
					case 19 : f[i] = discus(pop[i]); break;
					case 20 : f[i] = rosenbrock(pop[i]); break;
					case 21 : f[i] = ackley(pop[i]); break;
					case 22 : f[i] = weierstrass(pop[i]); break;
					case 23 : f[i] = griewank(pop[i]); break;
					case 24 : f[i] = rastrigin(pop[i]); break;
					case 25 : f[i] = schwefel(pop[i]); break;
					case 26 : f[i] = katsuura(pop[i]); break;
					case 27 : f[i] = happycat(pop[i]); break;
					case 28 : f[i] = hgbat(pop[i]); break;
					case 29 : f[i] = grie_rosen(pop[i]); break;
					case 30 : f[i] = escaffer6(pop[i]); break;
					default : System.out.println("没有指定函数可用！"); break;
				}			
			}
		}else{
			double[] x = new double[m * n];	
			testfunc tf = new testfunc ();
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					x[i*n + j] = pop[i][j];
				}
			}
			tf.test_func(x,f,n,m,func_num);
		}
		return f;
	}
	/**-----------------------------------------------------------------------**/
	/**-----------------------------------------------------------------------**/
	
	
	/** 
	 * 随机指定范围内N个不重复的数 
	 * 最简单最基本的方法 
	 * @param min 指定范围最小值 
	 * @param max 指定范围最大值 
	 * @param n 随机数个数 
	 */  
	public static int[] randomCommon(int min, int max, int n){  
	    if (n > (max - min + 1) || max < min) {  
	    	return null;  
	    }  
	    int[] result = new int[n];  
	    int count = 0;  
	    while(count < n) {  
	        int num = (int) (Math.random() * (max - min)) + min;
	        boolean flag = true;  
	        for (int j = 0; j < n; j++) {  
	            if(num == result[j]){  
	                flag = false;  
	                break;  
	            }  
	        }  
	        if(flag){  
	            result[count] = num;  
	            count++;  
	        }  
	    }
	    return result;  
	}
	
	/***************************************************************/
	//求一维数组的最小值位置
	public static int formin(double[] value){			
		int i = 0;
		int index = i;
		double min = value[i];
		while(++i < value.length){
			if(value[i] < min){
				min = value[i];
				index = i;
			}
		}
		return index;			
	}
	
	/**------------------------------------------------------------------------**/
	/**------------------将相关数据写入到文件writeInFile函数重载--------------------**/
	/**------------------------------------------------------------------------**/
	/**
	 * 将字符串写入到txt文件，每次会写入一行
	 * @param filename	定义的txt文件名
	 * @param str	需要写入的txt字符串
	 */
	public static void writeInFile(String filename, String str){
		try 
		{
			File file = new File(filename);
			if(!file.exists()){
				System.out.println(filename + " 文件不存在！正在创建文件...");
				file.createNewFile();
			}
			FileWriter fileWritter = new FileWriter(file.getAbsolutePath(),true);
			fileWritter.write(str + "\n");
			fileWritter.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * 将一个一维数组写入txt文件，每行写入一个数组
	 * @param filename	定义的txt文件名
	 * @param arr	需写入的一维数组
	 */
	public static void writeInFile(String filename, double[] arr){
		try 
		{
			File file = new File(filename);
			if(!file.exists()){
				System.out.println(filename + " 文件不存在！正在创建文件...");
				file.createNewFile();
			}
			FileWriter fileWritter = new FileWriter(file.getAbsolutePath(),true);
			for (int i = 0; i < arr.length; i++) {
				String str = String.valueOf(arr[i]);
				fileWritter.write(str + "\t");
			}
			fileWritter.write("\n");
			fileWritter.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * 将一个二维数组写入txt文件，按矩阵形式写入
	 * @param filename	定义的txt文件名
	 * @param matrix	需写入的二维矩阵
	 */
	public static void writeInFile(String filename, double[][] matrix){
		try 
		{
			File txt = new File(filename);
			if(!txt.exists()){
				System.out.println(filename + " 文件不存在！正在创建文件...");
				txt.createNewFile();
			}
			FileWriter fw = new FileWriter(txt.getAbsolutePath(),true);
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[i].length; j++) {
					String v3 = String.valueOf(matrix[i][j]);
					fw.write(v3 + "\t");
				}
				fw.write("\n");
			}
			fw.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
	

}
