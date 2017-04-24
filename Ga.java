package ga;

class geneAlgo {
	//gene algorithm initiate
	int geneLength = 20;
	int populationCnt = 30;
	int iteration = 50;
	double crossoverRate = 0.8;
	double mutationRate = 0.2;
	
	int gene[][] = new int[populationCnt][geneLength];
	int generation[][] = new int[populationCnt][geneLength];
	double decode[] = new double[populationCnt];
	double error[] = new double[populationCnt];
	
	int bestGene[] = new int [geneLength];
	double best  = 0;
	
	geneAlgo() {
				
	}
	
	//產生初始族群
	void initial() {
		double temp;
		for(int i = 0; i < geneLength ; i++) {
		    for(int j = 0; j < populationCnt; j++) {
		    	temp = Math.random();
		        if(temp > 0.5)
		            gene[j][i] = 1;
		        else
		            gene[j][i] = 0;
		    }
		}
	}
	
	//解碼
	void decode() {
		System.out.println("解碼");
		for(int i = 0; i < populationCnt; i++) {
			decode[i] = 0;
		}
		for(int i = 0; i < populationCnt; i++)
			for(int j = 0; j < geneLength; j++) //變數長度
		        decode[i] = decode[i] + gene[i][j] * Math.pow(2,j);
		
		for(int i = 0; i < populationCnt; i++)
			decode[i] = decode[i] * 0.001;
	}
	
	//計算適應程度
	void fitness(){
		System.out.println("計算適應程度");
		
	    for(int i = 0; i < populationCnt; i++) {    	
	    	error[i] = 0;    	
		}
	    for(int i = 0; i < populationCnt; i++) {	
	        error[i] = 1 / decode[i];
	    }
	}
	
	//競爭式選擇  複製
	void copy() {
	    System.out.println("競爭式選擇  複製");
	    int comp1, comp2;
	    int i = 0;
	    while(i != populationCnt) {
	    	comp1 = (int)(Math.random() * populationCnt);
	        comp2 = (int)(Math.random() * populationCnt);
	        //不可與自己競爭
	        while(comp1 == comp2) {
	            comp2 = (int)(Math.random() * populationCnt);
	        }
	        
        	if(error[comp1] > error[comp2]) {
	        	for(int j = 0; j < geneLength; j++) {
	                generation[i][j] = gene[comp2][j];	           
	            }
	        	i++;
	        } else {
        		for(int j = 0; j < geneLength; j++) {
	                generation[i][j] = gene[comp1][j];	           
	            }
	        	i++;
        	}
	    }
	}
	
	//交配  交換資訊
	void crossover() {
	    System.out.println("交配  交換資訊");
	    
	    //計算需交配的次數
	    double crossover = populationCnt * crossoverRate / 2;
	    
	    int crossoverChk[] = new int[populationCnt];
	    for(int i = 0; i < populationCnt; i++)
	    	crossoverChk[i] = 0;
	    
	    for(int i = 0; i < crossover; i++) {
	        //隨機選取2條未交配之相異染色體
	        int cross1 = (int)(Math.random() * populationCnt);
	        while(crossoverChk[cross1] != 0) {
	            cross1 = (int)(Math.random() * populationCnt);
	        }
	        
	        int cross2 = (int)(Math.random() * populationCnt);
	        while(cross1 == cross2 || crossoverChk[cross2] != 0) {
	            cross2 = (int)(Math.random() * populationCnt);
	        }
	        
	        //隨機選取2個相異交配點
	        int crosspoint1 = (int)(Math.random() * geneLength);
	        int crosspoint2 = (int)(Math.random() * geneLength);	        
	        while(crosspoint1 == crosspoint2) {
	            crosspoint2 = (int)(Math.random() * geneLength);
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
	void mutation(int iter) {
	    System.out.println("突變");
	    
	    if(iter!=iteration) { //最後一世代不突變
	    	
	    	//計算突變數量
	        double mutation = populationCnt * mutationRate;
	        
	        int mutationChk[] = new int[populationCnt];
	        for(int i = 0; i < populationCnt; i++)
	        	mutationChk[i] = 0;
	        
	        for(int i = 0; i < mutation; i++) {
	        	
	            //隨機選取1條未突變之染色體
	            int mu = (int)(Math.random() * populationCnt);
	            while(mutationChk[mu]!=0)
	                mu = (int)(Math.random() * populationCnt);
	            
	            //隨機選取1個突變點
	            int mupoint = (int)(Math.random() * geneLength);
	            
	            //突變
	            if(generation[mu][mupoint] == 0)
	                generation[mu][mupoint] = 1;
	            else
	                generation[mu][mupoint] = 0;
	            
	            //紀錄已突變之染色體
	            mutationChk[mu] = 1;
	        }
	    }
	}
	
	//取代上一世代
	void change() {
	    for(int i = 0;i<populationCnt;i++)
	        for(int j = 0;j<geneLength;j++)
	            gene[i][j] = generation[i][j];
	}
	
	//最後一代之間競爭
	void lastGen() {
		System.out.println("最後一代之間競爭");	
		
		decode();	
		
		fitness();
		
		double besterror = Integer.MAX_VALUE;
		
		//紀錄結果
		for(int i = 0; i < populationCnt; i++) {
			if(error[i] < besterror) {
				for(int j = 0; j < geneLength; j++)
					bestGene[j] = gene[i][j];
				best = decode[i];
				besterror = error[i];
			}							
		}
	}
	
	//訓練函式
	void train() {
		initial();
		for(int i = 0; i < iteration; i++) {
			System.out.println("gen: " + (i + 1));
			decode();
			fitness();
			copy();
			crossover();
			mutation(i + 1);
			change();
		}
		lastGen();		
	}
	
	//顯示結果
	void result() {
		System.out.println("best:" + best);
		for(int i = 0; i < geneLength; i++)
			System.out.print(bestGene[i] + " ");
	}
}

public class Ga {
	public static void main(String[] args) {
		geneAlgo ga = new geneAlgo();
		ga.train();
		ga.result();
	}
}