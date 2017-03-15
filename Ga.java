package ga;

import java.io.IOException;
import java.io.File;

class geneAlgo{
	//��Ǻt��k��l�ƾ�
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
	
	//���o���	
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
	//���ͪ�l�ڸs
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
	//�ѽX
	void decode(int imaxrte){
		System.out.println("�ѽX");
		for(int i = 0;i<populationCnt;i++){
			sd[i] = 0;
			cj[i] = 0;
		}
		for(int i = 0;i<populationCnt;i++){
			for(int j = 0;j<10;j++) //�Ĥ@�ܼƪ���
		        sd[i] = sd[i]+gene[i][j]*Math.pow(2,j);
		    for(int k = 0;k<11;k++) //�ĤG�ܼƪ���
		        cj[i] = cj[i]+gene[i][k+10]*Math.pow(2,k);
		}
		for(int i = 0;i<populationCnt;i++)
			sd[i] = sd[i]*0.001;
	}
	//�p��A���{��
	void fitness(int imaxrte){
		System.out.println("�p��A���{��");
	    double temp = 0;
	    double temp2 = 0;
	    for(int i = 0;i<populationCnt;i++){
	    	if(imaxrte==0)
	    		error0[i] = 0;
	    	else
	    		error1[i] = 0;
		}
	    for(int i = 0;i<populationCnt;i++){	
	        for(int j =0;j<270;j++){ //��Ƽƥإ����a�J
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
	        	error0[i] = Math.sqrt(error0[i]/270); //�觡��
	        else
	        	error1[i] = Math.sqrt(error1[i]/270); //�觡��
	    }
	}
	//�v�������  �ƻs
	void copy(int imaxrte){		
	    System.out.println("�v�������  �ƻs");
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
	//��t  �洫��T
	void crossover(){
	    System.out.println("��t  �洫��T");
	    double crossover = populationCnt*crossoverRate/2;
	    int crossoverChk[] = new int[populationCnt];
	    for(int i = 0;i<populationCnt;i++)
	    	crossoverChk[i] = 0;
	    for(int i = 0;i<crossover;i++){
	        //�H�����2������t���۲��V����
	        int cross1 = (int)(Math.random()*populationCnt);
	        while(crossoverChk[cross1]!=0){
	            cross1 = (int)(Math.random()*populationCnt);
	        }
	        int cross2 = (int)(Math.random()*populationCnt);
	        while(cross1==cross2||crossoverChk[cross2]!=0){
	            cross2 = (int)(Math.random()*populationCnt);
	        }
	        //�H�����2�Ӭ۲���t�I
	        int crosspoint1 = (int)(Math.random()*geneLength);
	        int crosspoint2 = (int)(Math.random()*geneLength);
	        while(crosspoint1==crosspoint2){
	            crosspoint2 = (int)(Math.random()*geneLength);
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
	void mutation(){
	    System.out.println("����");
	    if(time!=iteration){ //�̫�@�@�N������
	        double mutation = populationCnt*mutationRate;
	        int mutationChk[] = new int[populationCnt];
	        for(int i = 0;i<populationCnt;i++)
	        	mutationChk[i] = 0;
	        for(int i = 0;i<mutation;i++){
	            //�H�����1�������ܤ��V����
	            int mu = (int)(Math.random()*populationCnt);
	            while(mutationChk[mu]!=0)
	                mu = (int)(Math.random()*populationCnt);		            
	            //�H�����1�Ӭ����I
	            int mupoint = (int)(Math.random()*geneLength);
	            //����
	            if(generation[mu][mupoint]==0)
	                generation[mu][mupoint] = 1;
	            else
	                generation[mu][mupoint] = 0;
	            //�����w���ܤ��V����
	            mutationChk[mu] = 1;
	        }
	    }
	}
	//���N�W�@�@�N
	void change(){
	    for(int i = 0;i<populationCnt;i++)
	        for(int j = 0;j<geneLength;j++)
	            gene[i][j] = generation[i][j];
	}
	//�̫�@�N�����v��
	void lastGen(int imaxrte){
		System.out.println("�̫�@�N�����v��");
		
		//�ѽX
		System.out.println("�ѽX");
		for(int i = 0;i<populationCnt;i++){
			sd[i] = 0;
			cj[i] = 0;
		}
		for(int i = 0;i<populationCnt;i++){
			for(int j = 0;j<10;j++) //�Ĥ@�ܼƪ���
		        sd[i] = sd[i]+gene[i][j]*Math.pow(2,j);
		    for(int k = 0;k<11;k++) //�ĤG�ܼƪ���
		        cj[i] = cj[i]+gene[i][k+10]*Math.pow(2,k);
		}
		for(int i = 0;i<populationCnt;i++)
			sd[i] = sd[i]*0.001;

		
		//��X�̨�
		System.out.println("�p��A���{��");
	    double temp = 0;
	    double temp2 = 0;
	    for(int i = 0;i<populationCnt;i++){
	    	if(imaxrte==0)
	    		error0[i] = 0;
	    	else
	    		error1[i] = 0;
		}
	    for(int i = 0;i<populationCnt;i++){	
	    	for(int j =0;j<270;j++){ //��Ƽƥإ����a�J
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
	        	error0[i] = Math.sqrt(error0[i]/270); //�觡��
	        else
	        	error1[i] = Math.sqrt(error1[i]/270); //�觡��
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
		//�������G
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
		//������
		System.out.println("������");
		temp = 0;
		temp2 = 0;
		for(int j =0;j<270;j++){ //��Ƽƥإ����a�J
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
		imaxerror = Math.sqrt(imaxerror/270); //�觡��	
		for(int j =0;j<270;j++){ //��Ƽƥإ����a�J
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
		rteerror = Math.sqrt(rteerror/270); //�觡��
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
	//��ܵ��G
	void result(){
		System.out.println("��l���uimax�~�t�G "+imaxerror);
		System.out.println("��l���urte�~�t�G "+rteerror);
		System.out.println("���禱�uimax�~�t�G "+error0[best0]);
		System.out.println("���禱�urte�~�t�G "+error1[best1]);
		System.out.println("best imax cj�G "+cj0[best0]);
		System.out.println("best imax sd�G "+sd0[best0]);
		System.out.println("best rte cj�G "+cj1[best1]);
		System.out.println("best rte sd�G "+sd1[best1]);
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