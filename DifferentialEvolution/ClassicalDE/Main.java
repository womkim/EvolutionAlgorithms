import java.io.File;
import java.util.Random;


public class Main {
	
	
	public static double Xmin = Param.Xmin;
	public static double Xmax = Param.Xmax;	
	public static Random rand =  new Random();

	public static void main(String[] args) throws Exception {
		if(args.length < 8){
			System.out.println("\t!------------------------------------------------------!\n"
					+ "\t------------------- 请输入正确的参数!------------------\n"
					+ "\t--------- 第【一】个参数为【变异因子F】的值------------\n"
					+ "\t--------- 第【二】个参数为【交叉因子CR】的值-----------\n"
					+ "\t--------- 第【三】个参数为【变异策略】，为数字1 To 7---\n"
					+ "\t------------------ 1 -- DE/rand/1 ---------------------\n"
					+ "\t------------------ 2 -- DE/best/1 ---------------------\n"
					+ "\t------------------ 3 -- DE/rand/2 ---------------------\n"
					+ "\t------------------ 4 -- DE/best/2 ---------------------\n"
					+ "\t------------------ 5 -- DE/rand-To-best/1 -------------\n"
					+ "\t------------------ 6 -- DE/current-To-best/1 ----------\n"
					+ "\t------------------ 7 -- DE/rand-To-best/2 -------------\n"
					+ "\t-------- 第【四】个参数为【起始函数】，为数字1 To 30---\n"
					+ "\t-------- 第【五】个参数为【结束函数】，为数字1 To 30---\n"
					+ "\t-------- 接下来的参数为【问题的维度】，可以设置多个 ---\n"
					+ "\t--起始函数<15时，问题的维度只能设置为 10、30、50、100--\n"
					+ "\t-------- 倒数第二个参数为【函数在当前条件的计算次数】--\n"
					+ "\t-------- 最后一个参数为【目录名】，放置生成的文件 -----\n"
					+ "\t-------------- 种群大小默认为维度的10倍 ---------------\n"
					+ "\t!------------------------------------------------------!");
			System.exit(0);
		}
		double F = Double.parseDouble(args[0]);
		if(F > 1 || F < 0) System.out.println("变异因子的值应该在[0,1]范围内 :)");
		double CR = Double.parseDouble(args[1]);
		if(CR > 1 || CR < 0) System.out.println("交叉因子的值应该在[0,1]范围内 :)");
		int strategy = Integer.parseInt(args[2]);
		if(strategy > 7 || strategy < 1) System.out.println("变异的策略只设置七种，请输入1-7的整数 :)");
		int fstart = Integer.parseInt(args[3]);
		if(fstart > 30 || fstart < 1) System.out.println("开始的函数序号范围是[1,30]，请输入1-30的整数 :)");
		int fend = Integer.parseInt(args[4]);
		if(fstart > fend || fend > 30 || fend < 1) System.out.println("结束的函数序号在1-30之间，且不能小于开始的函数序号，请输入正确的整数 :)");
		int[] d = new int[args.length-7];
		for(int i = 5; i<args.length-2; i++){
			d[i-5]  = Integer.parseInt(args[i]);
		}
		int ll = Integer.parseInt(args[args.length-2]);
		if(ll < 1){
			System.out.println("不能计算少于1次!");
			System.exit(0);
		}else if(ll > 51){
			System.out.println("最大计算次数限制为51次！");
			System.exit(0);
		}
		String dir = args[args.length-1];
		
		//判断文件是否存在
		File dirfile = new File("./"+dir);
		if(!dirfile.exists() || !dirfile.isDirectory()){
			System.out.println(dir+"文件夹不存在！创建文件夹中");
			dirfile.mkdir();
		}		
		/*写入标题，只用写一次就够了*/
		String filename = "./"+dir+"/F"+F+"_CR"+CR+"_ResultRecord";
		String title = "func_num\t Np \t D \t G \t F \t CR \t Xmin \t Xmax \t bestvalue \t totalTime \t AvgTime\t";

		//开始测试DE
		
		//计算各种维度
		for (int k = 0; k < d.length; k++) {
			Param.writeInFile(filename+"_D_"+d[k]+".txt", title);
			Param.writeInFile(filename+"_valueCollect.txt", "D="+d[k]+"\tfunc_num\tMin\tMax\tMedian\tAvg\tStd");
			//计算各个函数
			for (int i = fstart; i <= fend; i++) {			
				double[] vv = new double[ll];	//用于求最好值fmin、最差值fmax、中值fmedia、平均值favg、标准差fstd
				for (int j = 0; j < ll; j++) {
					//！--- 主测试程序在这里 ---！
					vv[j] = classicalDE(d[k], d[k]*10, F, CR, d[k]*10000, i, strategy, filename);
				}
				
				//开始为每个维度的每个函数计算值分析
				double fmin = vv[0], fmax = vv[0], favg = 0, fmedia = vv[(int)(ll/2)], fstd = 0;
				for (int t = 0; t < ll; t++) {
					if(vv[t] < fmin)
						fmin = vv[t];
					if(vv[t] > fmax)
						fmax = vv[t];
					favg += vv[t]/ll;
				}
				for (int t = 0; t < ll; t++) {
					fstd += Math.pow((vv[t] - favg), 2);
				}
				fstd = Math.sqrt(fstd/ll);
				//写入计算结果
				String valueStr = "\tf"+i+"\t"+fmin+"\t"+fmax+"\t"+fmedia+"\t"+favg+"\t"+fstd;
				Param.writeInFile(filename+"_valueCollect.txt", valueStr);
			}
		}
	}
	
	/**
	 * 标准DE
	 * @param D				问题维度
	 * @param Np			种群大小
	 * @param Max_nfes		最大函数评估次数
	 * @param func_num		函数序号
	 * @param strategy		变异策略（1-7）
	 * @param filename		写入的文件名
	 * @throws Exception	抛出计算函数相关异常
	 */
	public static double classicalDE(int D, int Np, double F, double CR, int Max_nfes, int func_num, int strategy, String filename) throws Exception{
		
		double startTime = System.currentTimeMillis();		//记录程序开始时间
		
		double[] recordValue = new double[17];		//用于记录进化过程的最优值
		
		int G = 0;				//迭代次数
		int nfes = 0;			//函数评估次数
		double[] recp = {0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		double[] rec = new double[17];
		int recIndex = 0;
		for (int i = 0; i < recp.length; i++) {
			rec[i] = recp[i] * Max_nfes;	//用于当nfes运行到rec[i]次时，记录函数最优值
		}
		boolean stop_con = true;		//进化停止条件
		
		double bestvalue = 0;	//每代的最优值
		int bestIndex = 0;		//最优值位置
		
		//第G代种群
		double[][] XG = new double[Np][D];
		
		/************************************初始化种群**************************************/
		for(int i=0; i<Np; i++){
			for(int j=0; j<D; j++){
				XG[i][j] = (Xmax -Xmin) * rand.nextDouble() + Xmin;
			}
		}
		
		/**************************************评价初始种群*********************************/
		double[] origValue = Param.getFuncAllValue(XG, func_num);	//记录当代目标函数值
		bestIndex = Param.formin(origValue);	//记录初始种群的最优值位置
		
		//变异和交叉后的个体
		double[][] Xnext = new double[Np][D];
		double[] newValue = new double[Np];			//记录交叉后的目标函数值
		
		for (int i = 0; i < newValue.length; i++) {
			nfes++;
			for (int j = recIndex; j < recp.length; j++) {
				if(nfes == rec[j]){
					recordValue[recIndex] = (func_num < 16 ? (origValue[i]-func_num*100) : origValue[i]);
					recIndex++;
				}
			}
		}
		
		if(nfes < Max_nfes){
			stop_con = false;
		}
		
		/**
		 * 开始进化
		 */
		while(!stop_con){
			
			/*************************************变异操作+交叉操作********************************/
			for(int i=0; i<Np; i++){
				int[] ds = Param.randomCommon(0,Np,6);	//生成6个种群内的随机数，用于参与变异的个体序号
				
				while(ds[5] == i){		//保证最后一个随机数与当前个体序号不同。
					int[] q = Param.randomCommon(0,Np,1);
					ds[5] = q[0];
				}
				
				int j=ds[0], p=ds[1], k=ds[2], l=ds[3], n=ds[4];
				//保证选择参与变异的个体与当前个体不同
				if(i==j)	j=ds[5];
				else if(p==i)	p=ds[5];
				else if(k==i)	k=ds[5];
				else if(l==i)	l=ds[5];
				else if(n==i)	n=ds[5];
				
				int rd = rand.nextInt(D);	//生成当前维度范围的随机数
				//开始变异，交叉
				for (int a = 0; a < D; a++) {
					//二项交叉（这里跟书上的表达式是反的，不影响；XG是原个体，Xnext是变异和交叉后的个体）
					if(Math.random() > CR && rd != a){
						Xnext[i][a] = XG[i][a];	//原个体
					}else{
						switch(strategy){
							case 1 : Xnext[i][a] = XG[j][a] + F * (XG[p][a] - XG[k][a]); break;					// DE/rand/1
							case 2 : Xnext[i][a] = XG[bestIndex][a] + F * (XG[p][a] - XG[k][a]); break;			// DE/best/1
							case 3 : Xnext[i][a] = XG[j][a] + F * (XG[p][a] - XG[k][a]) + F * (XG[l][a] - XG[n][a]); break;					// DE/rand/2
							case 4 : Xnext[i][a] = XG[bestIndex][a] + F * (XG[p][a] - XG[k][a]) + F * (XG[l][a] - XG[n][a]); break;			// DE/best/2
							case 5 : Xnext[i][a] = XG[j][a] + F * (XG[bestIndex][a] - XG[j][a]) + F * (XG[p][a] - XG[k][a]); break;			// DE/rand-To-best/1
							case 6 : Xnext[i][a] = XG[i][a] + F * (XG[bestIndex][a] - XG[i][a]) + F * (XG[j][a] - XG[p][a]); break;			// DE/current-To-best/1
							case 7 : Xnext[i][a] = XG[i][a] + F * (XG[bestIndex][a] - XG[i][a]) + F * (XG[j][a] - XG[k][a]) + F * (XG[l][a] - XG[n][a]); break;			// DE/rand-To-best/2
							default : System.out.println("Strategy Wrong!!! Can not Evolution!!! Right Strategy Number is 1-7"); System.exit(0);
						}						
						
						//对超出边界处理，重新初始化
						if(Xnext[i][a] > Xmax || Xnext[i][a] < Xmin){
							Xnext[i][a] = (Xmax -Xmin) * rand.nextDouble() + Xmin;
						}
					}
				}
			}
			
			/*************************************评价函数计算*********************************/
			newValue = Param.getFuncAllValue(Xnext, func_num);
			for (int i = 0; i < newValue.length; i++) {
				nfes++;
				for (int j = recIndex; j < recp.length; j++) {
					if(nfes == rec[j]){
						//这里记录的是nfes值为rec[j]时的评价函数值，不一定是每一代的最优值
						recordValue[recIndex] = (func_num < 16 ? (newValue[i]-func_num*100) : newValue[i]);
						recIndex++;
					}
				}			
			}
			
			/*************************************选择操作**************************************/
			for (int i = 0; i < Np; i++) {
				if(newValue[i] < origValue[i]){
					origValue[i] = newValue[i];			//新的评价值替换
					for (int j = 0; j < D; j++) {
						XG[i][j] = Xnext[i][j];			//个体替换
					}
				}
			}
			
			/**
			 * 记录最优值
			 */
			bestIndex = Param.formin(origValue);
			bestvalue = (func_num < 16 ? (origValue[bestIndex]-func_num*100) : origValue[bestIndex]);
//			bestvalue = (func_num < 16 ? (bestvalue < 1e-08 ? 0.0 : bestvalue) : bestvalue);	//Error value smaller than 1e-08 will be taken as zero!
			
			//代数增加1
			G++;
			//当达到最大函数评估次数或者已经找到最优值时，结束循环
			if(nfes > Max_nfes || bestvalue == 0){
				stop_con = true;
			}
			
		}//进化结束
		
		double endTime = System.currentTimeMillis();
		double totalTime = (endTime - startTime)/1000;		//计算总时间
		System.out.println("函数f"+func_num+" 维度为："+D+" 的进化过程总时间花费："+totalTime+"s");
		
		//写入当前进化信息，包括<函数序号>、<种群大小>、<维度>、<最大迭代次数>、<变异因子F的值>、<交叉因子CR的值>、<函数下界>、<函数上界>、<进化的最优值>、<进化过程总时间>、<评价每代进化时间>
		String outstr = "f"+func_num+"\t"+Np+"\t"+D+"\t"+G+"\t"+F+"\t"+CR+"\t"+Xmin+"\t"+Xmax+"\t"+bestvalue+"\t"+totalTime+"\t"+(totalTime/G)+"\t";
		Param.writeInFile(filename+"_D_"+D+".txt", outstr);
		//写入进化过程中记录的评价函数值（设置的17个截点）
		Param.writeInFile(filename+"_f"+func_num+"_D"+D+"_BestValue.txt", recordValue);
		
		//返回当前进化的最优值
		return bestvalue;
	}
	
	
}
