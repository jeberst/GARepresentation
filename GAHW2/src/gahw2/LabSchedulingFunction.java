
import java.io.FileWriter;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Chathika
 */
class LabSchedulingFunction extends FitnessFunction {
 
	int[][][] preferences; // [Person][Day][Shift]
	static int p1 = 100, 
			p2 = 50, 
			p3 = 25, 
			p4 = 5, 
			p_invalid = -3400, 
			shift_penalty = -700, // TODO: We need to discuss a value for this penalty
			shifts_per_day = 5;

  public LabSchedulingFunction() {
      name = "Lab Scheduling Function";
      initialize_preferences();
  }  
  
  public void doRawFitness(Chromo X){
	  	
	  	// Start with 0 as the raw fitness
		X.rawFitness = 0;
		
		// Initialize the number of shifts to 0 for each employee
		int[] numShifts = new int[8];
		
		for (int i = 0; i < numShifts.length; i ++){
			numShifts[i] = 0;
		}
		
		// Loop through the genes in the chromosome
		for (int z = 0; z < X.selections.length; z++){
			
			// Set the array-indexing variables
			int worker = X.selections[z]; // integer representation of gene from Chromo
			int day = z / shifts_per_day; // gets the day of the week
			int shift = z % shifts_per_day; // gets the shift of the day
			
			// Check the preference value, reward or penalize the score accordingly
			switch (preferences[worker][day][shift]) {
				case 0:
					X.rawFitness += p_invalid;
					break;
				case 1:
					X.rawFitness += p1;
					break;
				case 2:
					X.rawFitness += p2;
					break;
				case 3:
					X.rawFitness += p3;
					break;
				case 4:
					X.rawFitness += p4;
					break;
				default:
					break;
			}
			
			// Track the shifts for this worker
			numShifts[worker]++;
		}
		
		int shift_discrepencies = numShifts[0]; // adds the shifts assigned to "Employee 0", aka Nobody
		
		for(int k = 1; k < numShifts.length; k++){
			shift_discrepencies += Math.abs(shifts_per_day - numShifts[k]); // adds the absolute difference between actual and expected shift numbers
		}
		
		 X.rawFitness += shift_discrepencies * shift_penalty; // Penalizes the score based on the number of shifts
	}

//PRINT OUT AN INDIVIDUAL GENE TO THE SUMMARY FILE *********************************

	public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{
		// TODO: adjust this function to handle various representations
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getGeneAlpha(i),11,output);
		}
		output.write("   RawFitness");
		output.write("\n        ");
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getPosIntGeneValue(i),11,output);
		}
		Hwrite.right((int) X.rawFitness,13,output);
		output.write("\n\n");
		return;
	}
	
	private void initialize_preferences(){
		//TODO: write code to initialize preference array
	}
	}
 
