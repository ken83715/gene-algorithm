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
	
	//���ͪ�l�ڸs
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
	
	//�ѽX
	void decode() {
		System.out.println("�ѽX");
		for(int i = 0; i < populationCnt; i++) {
			decode[i] = 0;
		}
		for(int i = 0; i < populationCnt; i++)
			for(int j = 0; j < geneLength; j++) //�ܼƪ���
		        decode[i] = decode[i] + gene[i][j] * Math.pow(2,j);
		
		for(int i = 0; i < populationCnt; i++)
			decode[i] = decode[i] * 0.001;
	}
	
	//�p��A���{��
	void fitness(){
		System.out.println("�p��A���{��");
		
	    for(int i = 0; i < populationCnt; i++) {    	
	    	error[i] = 0;    	
		}
	    for(int i = 0; i < populationCnt; i++) {	
	        error[i] = 1 / decode[i];
	    }
	}
	
	//�v�������  �ƻs
	void copy() {
	    System.out.println("�v�������  �ƻs");
	    int comp1, comp2;
	    int i = 0;
	    while(i != populationCnt) {
	    	comp1 = (int)(Math.random() * populationCnt);
	        comp2 = (int)(Math.random() * populationCnt);
	        //���i�P�ۤv�v��
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
	
	//��t  �洫��T
	void crossover() {
	    System.out.println("��t  �洫��T");
	    
	    //�p��ݥ�t������
	    double crossover = populationCnt * crossoverRate / 2;
	    
	    int crossoverChk[] = new int[populationCnt];
	    for(int i = 0; i < populationCnt; i++)
	    	crossoverChk[i] = 0;
	    
	    for(int i = 0; i < crossover; i++) {
	        //�H�����2������t���۲��V����
	        int cross1 = (int)(Math.random() * populationCnt);
	        while(crossoverChk[cross1] != 0) {
	            cross1 = (int)(Math.random() * populationCnt);
	        }
	        
	        int cross2 = (int)(Math.random() * populationCnt);
	        while(cross1 == cross2 || crossoverChk[cross2] != 0) {
	            cross2 = (int)(Math.random() * populationCnt);
	        }
	        
	        //�H�����2�Ӭ۲���t�I
	        int crosspoint1 = (int)(Math.random() * geneLength);
	        int crosspoint2 = (int)(Math.random() * geneLength);	        
	        while(crosspoint1 == crosspoint2) {
	            crosspoint2 = (int)(Math.random() * geneLength);
	        }
	        
	        //�洫��T
	        int t;
	        t = generation[cross1][crosspoint1];
	        generation[cross1][crosspoint1] = generation[cross2][crosspoint1];
	        generation[cross2][crosspoint1] = t;
	        
	        t = generation[cross1][crosspoint2];
	        generation[cross1][crosspoint2] = generation[cross2][crosspoint2];
	        generation[cross2][crosspoint2] = t;
	        
	        //�����w��t���V����
	        crossoverChk[cross1] = 1;
	        crossoverChk[cross2] = 1;
	    }
	}
	
	//����
	void mutation(int iter) {
	    System.out.println("����");
	    
	    if(iter!=iteration) { //�̫�@�@�N������
	    	
	    	//�p����ܼƶq
	        double mutation = populationCnt * mutationRate;
	        
	        int mutationChk[] = new int[populationCnt];
	        for(int i = 0; i < populationCnt; i++)
	        	mutationChk[i] = 0;
	        
	        for(int i = 0; i < mutation; i++) {
	        	
	            //�H�����1�������ܤ��V����
	            int mu = (int)(Math.random() * populationCnt);
	            while(mutationChk[mu]!=0)
	                mu = (int)(Math.random() * populationCnt);
	            
	            //�H�����1�Ӭ����I
	            int mupoint = (int)(Math.random() * geneLength);
	            
	            //����
	            if(generation[mu][mupoint] == 0)
	                generation[mu][mupoint] = 1;
	            else
	                generation[mu][mupoint] = 0;
	            
	            //�����w���ܤ��V����
	            mutationChk[mu] = 1;
	        }
	    }
	}
	
	//���N�W�@�@�N
	void change() {
	    for(int i = 0;i<populationCnt;i++)
	        for(int j = 0;j<geneLength;j++)
	            gene[i][j] = generation[i][j];
	}
	
	//�̫�@�N�����v��
	void lastGen() {
		System.out.println("�̫�@�N�����v��");	
		
		decode();	
		
		fitness();
		
		double besterror = Integer.MAX_VALUE;
		
		//�������G
		for(int i = 0; i < populationCnt; i++) {
			if(error[i] < besterror) {
				for(int j = 0; j < geneLength; j++)
					bestGene[j] = gene[i][j];
				best = decode[i];
				besterror = error[i];
			}							
		}
	}
	
	//�V�m�禡
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
	
	//��ܵ��G
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