/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.util.Arrays;

public class Chromo implements Comparable<Chromo>
{
/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	public String chromo;
	public double rawFitness;
	public double sclFitness;
	public double proFitness;
	public int[] selections;
        private double[] randomArray = new double[Parameters.geneSize];
        private double[] sortedRandomArray = new double[Parameters.geneSize];
        private char[] newChromo = new char[Parameters.geneSize];
       
     

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	private static double randnum;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public Chromo(){

		//  Set gene values to a randum sequence of 1's and 0's
		char geneBit;
		chromo = "";
		for (int i=0; i<Parameters.numGenes; i++){
			for (int j=0; j<Parameters.geneSize; j++){
				randnum = Search.r.nextDouble();
				if (randnum > 0.5) geneBit = '0';
				else geneBit = '1';
				this.chromo = chromo + geneBit;
			}
		}

		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized
	}


/*******************************************************************************
*                                MEMBER METHODS                                *
     * @param chromo
*******************************************************************************/
  
	
        //  Get Alpha Represenation of a Gene **************************************
       
            
	public String getGeneAlpha(int geneID){
		int start = geneID * Parameters.geneSize;
		int end = (geneID+1) * Parameters.geneSize;
		String geneAlpha = this.chromo.substring(start, end);
		return (geneAlpha);
	}

	//  Get Integer Value of a Gene (Positive or Negative, 2's Compliment) ****

	public int getIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneSign;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=1; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		geneSign = geneAlpha.charAt(0);
		if (geneSign == '1') geneValue = geneValue - (int)Math.pow(2.0, Parameters.geneSize-1);
		return (geneValue);
	}

	//  Get Integer Value of a Gene (Positive only) ****************************

	public int getPosIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=0; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		return (geneValue);
	}

	//  Mutate a Chromosome Based on Mutation Type *****************************

	public void doMutation(){

		String mutChromo = "";
		char x;

		switch (Parameters.mutationType){

		case 1:     //  Replace with new random number

			for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
				x = this.chromo.charAt(j);
				randnum = Search.r.nextDouble();
				if (randnum < Parameters.mutationRate){
					if (x == '1') x = '0';
					else x = '1';
				}
				mutChromo = mutChromo + x;
			}
			this.chromo = mutChromo;
			break;

		default:
			System.out.println("ERROR - No mutation method selected");
		}
	}
        
        public double[] randomizeChromo()
        {
        double[] randomArray = new double[Parameters.geneSize];

            for(int i=0; i< this.chromo.length(); i++)
            {
                randomArray[i] = Search.r.nextDouble();
                
            }
            return randomArray;
            
        }
        


/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	//  Select a parent for crossover ******************************************

	public static int selectParent(){

		double rWheel = 0;
		int j = 0;
		int k = 0;
                Chromo[] rank;
              
		switch (Parameters.selectType){

		case 1:     // Proportional Selection
			randnum = Search.r.nextDouble();
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				if (randnum < rWheel) return(j);
			}
			break;

		case 3:     // Random Selection
			randnum = Search.r.nextDouble();
			j = (int) (randnum * Parameters.popSize);
			return(j);

		case 2:     //  Tournament Selection     
                        randnum = Search.r.nextDouble();
                        int firstNum = (int)(randnum * Parameters.popSize);
                        randnum = Search.r.nextDouble();
                        int secondNum = (int)(randnum * Parameters.popSize);
                        randnum = Search.r.nextDouble();
                        
                        int moreFit = 0;
                        int lessFit = 0;
                        if(Search.member[firstNum].rawFitness > Search.member[secondNum].rawFitness)
                        {
                            moreFit = firstNum;
                            lessFit = secondNum;
                        }
                        else
                        {
                            moreFit = secondNum;
                            lessFit = firstNum;
                        }
                        
                        double probability = (double).025/10;
                        if(randnum > probability){
                            return moreFit;
                        }
                        else
                        {
                            return lessFit;
                        }
                        
                    
                case 4: //Rank Selection
                     double rankSum = 0;
                     rWheel = 0;
                     double[] probabilityArray = new double[Search.member.length+1];
                     Arrays.sort(Search.member);
                    for (j=1; j<=Search.member.length; j++){
				rankSum = rankSum+j;
                                probabilityArray[j] = 0;
			}
                    
                    for(int i = 1; i<=Search.member.length; i++)
                    {
                        probabilityArray[i] = (double)(i)/rankSum;
                    }
                        randnum = Search.r.nextDouble();
			for (j=0; j<Search.member.length; j++){
				rWheel = rWheel + probabilityArray[j+1];
				if (randnum < rWheel) {
                                    return(j);
                                }
			}
                break;

		default:
			System.out.println("ERROR - No selection method selected");
		}
	return(-1);
	}

	//  Produce a new child from two parents  **********************************

	public static void mateParents(int pnum1, int pnum2, Chromo parent1, Chromo parent2, Chromo child1, Chromo child2){

		int xoverPoint1;
		int xoverPoint2;

		switch (Parameters.xoverType){

		case 1:     //  Single Point Crossover

			//  Select crossover point
			xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));

			//  Create child chromosome from parental material
			child1.chromo = parent1.chromo.substring(0,xoverPoint1) + parent2.chromo.substring(xoverPoint1);
			child2.chromo = parent2.chromo.substring(0,xoverPoint1) + parent1.chromo.substring(xoverPoint1);
			break;

		case 2:     //  Two Point Crossover

		case 3:     //  Uniform Crossover

		default:
			System.out.println("ERROR - Bad crossover method selected");
		}

		//  Set fitness values back to zero
		child1.rawFitness = -1;   //  Fitness not yet evaluated
		child1.sclFitness = -1;   //  Fitness not yet scaled
		child1.proFitness = -1;   //  Fitness not yet proportionalized
		child2.rawFitness = -1;   //  Fitness not yet evaluated
		child2.sclFitness = -1;   //  Fitness not yet scaled
		child2.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Produce a new child from a single parent  ******************************

	public static void mateParents(int pnum, Chromo parent, Chromo child){

		//  Create child chromosome from parental material
		child.chromo = parent.chromo;

		//  Set fitness values back to zero
		child.rawFitness = -1;   //  Fitness not yet evaluated
		child.sclFitness = -1;   //  Fitness not yet scaled
		child.proFitness = -1;   //  Fitness not yet proportionalized
	}

        
        public static void mateParents(int pnum1, int pnum2, Chromo parent1, Chromo parent2, Chromo child1, Chromo child2, double[] RandomKeysParent1, double[] RandomKeysParent2){

		int xoverPoint1;
		int xoverPoint2;

		switch (Parameters.xoverType){

		case 1:     //  Single Point Crossover
                        double[] randomChild1 = new double[RandomKeysParent1.length];
                        double[] randomChild2 = new double[RandomKeysParent2.length];
			//  Select crossover point
			xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));

                        
                        System.arraycopy(RandomKeysParent1, 0, randomChild1, 0, xoverPoint1);

                        System.arraycopy(RandomKeysParent2, xoverPoint1, randomChild1, xoverPoint1, 1);
                        
                         System.arraycopy(RandomKeysParent2, 0, randomChild2, 0, xoverPoint1);

                        System.arraycopy(RandomKeysParent1, xoverPoint1, randomChild2, xoverPoint1, 1);
                        
                        
                       
                        
                        
			//  Create child chromosome from parental material
			child1.chromo = determineChromo(parent1, randomChild1);
			child2.chromo = determineChromo(parent2, randomChild2);
			break;

		case 2:     //  Two Point Crossover

		case 3:     //  Uniform Crossover

		default:
			System.out.println("ERROR - Bad crossover method selected");
		}

		//  Set fitness values back to zero
		child1.rawFitness = -1;   //  Fitness not yet evaluated
		child1.sclFitness = -1;   //  Fitness not yet scaled
		child1.proFitness = -1;   //  Fitness not yet proportionalized
		child2.rawFitness = -1;   //  Fitness not yet evaluated
		child2.sclFitness = -1;   //  Fitness not yet scaled
		child2.proFitness = -1;   //  Fitness not yet proportionalized
	}
	//  Copy one chromosome to another  ***************************************

	public static void copyB2A (Chromo targetA, Chromo sourceB){

		targetA.chromo = sourceB.chromo;

		targetA.rawFitness = sourceB.rawFitness;
		targetA.sclFitness = sourceB.sclFitness;
		targetA.proFitness = sourceB.proFitness;
		return;
	}
        
                public static String determineChromo(Chromo oldChromo, double[] randomKey)
        {
        double[] sortedRandomArray = new double[Parameters.geneSize];
         char[] newChromo = new char[Parameters.geneSize];
        
            for(int i=0; i< oldChromo.chromo.length(); i++)
            {
                sortedRandomArray[i] = randomKey[i];
            }
            
            Arrays.sort(sortedRandomArray);
            
            for(int j=0; j<oldChromo.chromo.length(); j++)
            {
                for (int k =0; k<oldChromo.chromo.length(); k++)
                {
                    if(sortedRandomArray[j] == randomKey[k])
                    {
                        newChromo[j] = oldChromo.chromo.charAt(k);
                    }
                }
            }
            String returnChromo = "";
            for(int z=0; z<newChromo.length; z++)
            {
                returnChromo += newChromo[z];
            }
            return returnChromo;
        }
    @Override
    public int compareTo(Chromo chrome) {
        
        return (int)(rawFitness - chrome.rawFitness);
    }

        
    

}   // End of Chromo.java ******************************************************
