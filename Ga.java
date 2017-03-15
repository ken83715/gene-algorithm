package ga;

import java.io.IOException;
import java.io.File;

class geneAlgo{
	//遺傳演算法初始數據
	int geneLength = 21;
	int populationCnt = 30;
	int iteration = 100;
	double crossoverRate = 0.8;
	double mutationRate = 0.2;
	
	int time = 0;
	int gene[][] = new int[populationCnt][geneLength];
	int generation[][] = new int[populationCnt][geneLength];
	double sd[] = new double[populationCnt];
	double cj[] = new double[populationCnt];
	double error[] = new double[populationCnt];
	
	double sd0[] = new double[populationCnt];
	double cj0[] = new double[populationCnt];
	double sd1[] = new double[populationCnt];
	double cj1[] = new double[populationCnt];
	double error0[] = new double[populationCnt]; //imax
	double error1[] = new double[populationCnt]; //rte
	
	int best0 = 1;
	int best1 = 1;
	int skyline0[] = new int[populationCnt];
	int skyline1[] = new int[populationCnt];
	
	double rtecj = 1106;
	double rtesd = 0.591;
	double imaxcj = 98;
	double imaxsd = 0.561;
	double rteerror;
	double imaxerror;
	
	//取得資料	
	String Data[][] = new String[270][3];
	double doubleData[][] = new double[270][3];
	
	geneAlgo() {
		
		for(int i = 0;i<270;i++)
			for(int j = 0;j<3;j++)
				Data[i][j] = "0";
		for(int i = 0;i<270;i++)
			for(int j = 0;j<3;j++)
				doubleData[i][j] = Double.parseDouble(Data[i][j]);
	}
	//產生初始族群
	void initial(){
		double temp;
		for(int i = 0;i<geneLength;i++){
		    for(int j = 0;j<populationCnt;j++){
		    	temp = Math.random();
		        if(temp>0.5)
		            gene[j][i] = 1;
		        else
		            gene[j][i] = 0;
		    }
		}
	}	
	//解碼
	void decode(int imaxrte){
		System.out.println("解碼");
		for(int i = 0;i<populationCnt;i++){
			sd[i] = 0;
			cj[i] = 0;
		}
		for(int i = 0;i<populationCnt;i++){
			for(int j = 0;j<10;j++) //第一變數長度
		        sd[i] = sd[i]+gene[i][j]*Math.pow(2,j);
		    for(int k = 0;k<11;k++) //第二變數長度
		        cj[i] = cj[i]+gene[i][k+10]*Math.pow(2,k);
		}
		for(int i = 0;i<populationCnt;i++)
			sd[i] = sd[i]*0.001;
	}
	//計算適應程度
	void fitness(int imaxrte){
		System.out.println("計算適應程度");
	    double temp = 0;
	    double temp2 = 0;
	    for(int i = 0;i<populationCnt;i++){
	    	if(imaxrte==0)
	    		error0[i] = 0;
	    	else
	    		error1[i] = 0;
		}
	    for(int i = 0;i<populationCnt;i++){	
	        for(int j =0;j<270;j++){ //資料數目全部帶入
	        	double t = Math.log(doubleData[j][imaxrte]/cj[i])/(sd[i]*Math.sqrt(2));
	        	double d = 2*t/10000;
	        	double x = -t;
	        	for(int k = 0;k<10000;k++){
	        		temp = temp+Math.exp(-x*x);
	        		x = x+d;
	        	}
	        	temp = temp*d;
	        	temp = 1/2+temp/2;
		        if(temp>doubleData[j][2])
		        	temp2 = temp-doubleData[j][2];
		        else
		        	temp2 = doubleData[j][2]-temp;
		        if(imaxrte==0)
		        	error0[i] = error0[i]+temp2;	
		        else
		        	error1[i] = error1[i]+temp2;
	        }
	        if(imaxrte==0)
	        	error0[i] = Math.sqrt(error0[i]/270); //方均根
	        else
	        	error1[i] = Math.sqrt(error1[i]/270); //方均根
	    }
	}
	//競爭式選擇  複製
	void copy(int imaxrte){		
	    System.out.println("競爭式選擇  複製");
	    int comp1,comp2;
	    int i = 0;
	    while(i!=populationCnt){
	    	comp1 = (int)(Math.random()*populationCnt);
	        comp2 = (int)(Math.random()*populationCnt);
	        while(comp1==comp2){
	            comp2 = (int)(Math.random()*populationCnt);
	        }
	        if(imaxrte==0){
	        	if(error0[comp1]>error0[comp2]){
		        	for(int j = 0;j<geneLength;j++){
		                generation[i][j] = gene[comp2][j];	           
		            }
		        	i++;
		        }
	        	else{
	        		for(int j = 0;j<geneLength;j++){
		                generation[i][j] = gene[comp1][j];	           
		            }
		        	i++;
	        	}	        		
	        }
	        else{
	        	if(error1[comp1]>error1[comp2]){
		        	for(int j = 0;j<geneLength;j++){
		                generation[i][j] = gene[comp2][j];	           
		            }
		        	i++;
		        }
	        	else{
	        		for(int j = 0;j<geneLength;j++){
		                generation[i][j] = gene[comp1][j];	           
		            }
		        	i++;
	        	}
	        }
	    }
	}
	//交配  交換資訊
	void crossover(){
	    System.out.println("交配  交換資訊");
	    double crossover = populationCnt*crossoverRate/2;
	    int crossoverChk[] = new int[populationCnt];
	    for(int i = 0;i<populationCnt;i++)
	    	crossoverChk[i] = 0;
	    for(int i = 0;i<crossover;i++){
	        //隨機選取2條未交配之相異染色體
	        int cross1 = (int)(Math.random()*populationCnt);
	        while(crossoverChk[cross1]!=0){
	            cross1 = (int)(Math.random()*populationCnt);
	        }
	        int cross2 = (int)(Math.random()*populationCnt);
	        while(cross1==cross2||crossoverChk[cross2]!=0){
	            cross2 = (int)(Math.random()*populationCnt);
	        }
	        //隨機選取2個相異交配點
	        int crosspoint1 = (int)(Math.random()*geneLength);
	        int crosspoint2 = (int)(Math.random()*geneLength);
	        while(crosspoint1==crosspoint2){
	            crosspoint2 = (int)(Math.random()*geneLength);
	        }
	        //交換資訊
	        int t;
	        t = generation[cross1][crosspoint1];
	        generation[cross1][crosspoint1] = generation[cross2][crosspoint1];
	        generation[cross2][crosspoint1] = t;
	        t = generation[cross1][crosspoint2];
	        generation[cross1][crosspoint2] = generation[cross2][crosspoint2];
	        generation[cross2][crosspoint2] = t;
	        //紀錄已交配之染色體
	        crossoverChk[cross1] = 1;
	        crossoverChk[cross2] = 1;
	    }
	}
	//突變
	void mutation(){
	    System.out.println("突變");
	    if(time!=iteration){ //最後一世代不突變
	        double mutation = populationCnt*mutationRate;
	        int mutationChk[] = new int[populationCnt];
	        for(int i = 0;i<populationCnt;i++)
	        	mutationChk[i] = 0;
	        for(int i = 0;i<mutation;i++){
	            //隨機選取1條未突變之染色體
	            int mu = (int)(Math.random()*populationCnt);
	            while(mutationChk[mu]!=0)
	                mu = (int)(Math.random()*populationCnt);		            
	            //隨機選取1個突變點
	            int mupoint = (int)(Math.random()*geneLength);
	            //突變
	            if(generation[mu][mupoint]==0)
	                generation[mu][mupoint] = 1;
	            else
	                generation[mu][mupoint] = 0;
	            //紀錄已突變之染色體
	            mutationChk[mu] = 1;
	        }
	    }
	}
	//取代上一世代
	void change(){
	    for(int i = 0;i<populationCnt;i++)
	        for(int j = 0;j<geneLength;j++)
	            gene[i][j] = generation[i][j];
	}
	//最後一代之間競爭
	void lastGen(int imaxrte){
		System.out.println("最後一代之間競爭");
		
		//解碼
		System.out.println("解碼");
		for(int i = 0;i<populationCnt;i++){
			sd[i] = 0;
			cj[i] = 0;
		}
		for(int i = 0;i<populationCnt;i++){
			for(int j = 0;j<10;j++) //第一變數長度
		        sd[i] = sd[i]+gene[i][j]*Math.pow(2,j);
		    for(int k = 0;k<11;k++) //第二變數長度
		        cj[i] = cj[i]+gene[i][k+10]*Math.pow(2,k);
		}
		for(int i = 0;i<populationCnt;i++)
			sd[i] = sd[i]*0.001;

		
		//找出最佳
		System.out.println("計算適應程度");
	    double temp = 0;
	    double temp2 = 0;
	    for(int i = 0;i<populationCnt;i++){
	    	if(imaxrte==0)
	    		error0[i] = 0;
	    	else
	    		error1[i] = 0;
		}
	    for(int i = 0;i<populationCnt;i++){	
	    	for(int j =0;j<270;j++){ //資料數目全部帶入
	        	double t = Math.log(doubleData[j][imaxrte]/cj[i])/(sd[i]*Math.sqrt(2));
	        	double d = 2*t/10000;
	        	double x = -t;
	        	for(int k = 0;k<10000;k++){
	        		temp = temp+Math.exp(-x*x);
	        		x = x+d;
	        	}
	        	temp = temp*d;
	        	temp = 1/2+temp/2;
		        if(temp>doubleData[j][2])
		        	temp2 = temp-doubleData[j][2];
		        else
		        	temp2 = doubleData[j][2]-temp;
		        if(imaxrte==0)
		        	error0[i] = error0[i]+temp2;	
		        else
		        	error1[i] = error1[i]+temp2;
	        }
	        if(imaxrte==0)
	        	error0[i] = Math.sqrt(error0[i]/270); //方均根
	        else
	        	error1[i] = Math.sqrt(error1[i]/270); //方均根
	    }
		for(int i = 0;i<populationCnt;i++){
			if(imaxrte==0){
				if(error0[i]<error0[best0])
					best0 = i;
			}
			else{
			    if(error1[i]<error1[best1])
			    	best1 = i;
			}
		}
		//紀錄結果
		for(int i = 0;i<populationCnt;i++){
			if(imaxrte==0){
				sd0[i] = sd[i];
				cj0[i] = cj[i];
			}
			else{
				sd1[i] = sd[i];
				cj1[i] = cj[i];
			}
		}
		//比較兩者
		System.out.println("比較兩者");
		temp = 0;
		temp2 = 0;
		for(int j =0;j<270;j++){ //資料數目全部帶入
        	double t = Math.log(doubleData[j][0]/imaxcj)/(imaxsd*Math.sqrt(2));
        	double d = 2*t/10000;
        	double x = -t;
        	for(int k = 0;k<10000;k++){		
        		temp = temp+Math.exp(-x*x);
        		x = x+d;
        	}
        	temp = temp*d;
        	temp = 1/2+temp/2;
	        if(temp>doubleData[j][2])
	        	temp2 = temp-doubleData[j][2];
	        else
	        	temp2 = doubleData[j][2]-temp;
	        imaxerror = imaxerror+temp2;
        }
		imaxerror = Math.sqrt(imaxerror/270); //方均根	
		for(int j =0;j<270;j++){ //資料數目全部帶入
        	double t = Math.log(doubleData[j][1]/rtecj)/(rtesd*Math.sqrt(2));
        	double d = 2*t/10000;
        	double x = -t;
        	for(int k = 0;k<10000;k++){		
        		temp = temp+Math.exp(-x*x);
        		x = x+d;
        	}
        	temp = temp*d;
        	temp = 1/2+temp/2;
	        if(temp>doubleData[j][2])
	        	temp2 = temp-doubleData[j][2];
	        else
	        	temp2 = doubleData[j][2]-temp;
	        rteerror = rteerror+temp2;
        }
		rteerror = Math.sqrt(rteerror/270); //方均根
	}
	//SKYLINE
	void skyline(){
		for(int i = 0;i<populationCnt;i++){ 
			if(error0[i]>error0[best0]&&error1[i]<error1[best0])
				skyline0[i] = 1;
			if(error1[i]>error1[best1]&&error0[i]<error0[best1])
				skyline1[i] = 1;
		}
	}
	void train(){
		initial();
		for(int i = 0;i<iteration;i++){
			System.out.println("gen: "+(i+1));
			decode(0);
			fitness(0);
			copy(0);
			crossover();
			mutation();
			change();
		}
		lastGen(0);
		
		initial();
		for(int i = 0;i<iteration;i++){
			System.out.println("gen: "+(i+1));
			decode(1);
			fitness(1);
			copy(1);
			crossover();
			mutation();
			change();
		}
		lastGen(1);
		
		skyline();
	}
	//顯示結果
	void result(){
		System.out.println("原始曲線imax誤差： "+imaxerror);
		System.out.println("原始曲線rte誤差： "+rteerror);
		System.out.println("實驗曲線imax誤差： "+error0[best0]);
		System.out.println("實驗曲線rte誤差： "+error1[best1]);
		System.out.println("best imax cj： "+cj0[best0]);
		System.out.println("best imax sd： "+sd0[best0]);
		System.out.println("best rte cj： "+cj1[best1]);
		System.out.println("best rte sd： "+sd1[best1]);
		System.out.println(" ");
		for(int i = 0;i<populationCnt;i++)
			System.out.println(cj0[i]+" "+sd0[i]+" "+error0[i]+" "+skyline0[i]);
		System.out.println(" ");
		for(int i = 0;i<populationCnt;i++)
			System.out.println(cj1[i]+" "+sd1[i]+" "+error1[i]+" "+skyline1[i]);
	}
}

public class Ga {
	public static void main(String[] args) {
		geneAlgo ga = new geneAlgo();
		ga.train();
		ga.result();
	}
}