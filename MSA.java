import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;


public class MSA {

public MSA(){
	
}
	
	public void runMSA(ArrayList<String> sequencesArray, ArrayList<String> seqDescriptorsArray){
		
		Scanner scanner = new Scanner(System.in);
		
		ArrayList<String> sequences = new ArrayList<String>(sequencesArray);
		ArrayList<String> sequencesStored = new ArrayList<String>(sequences);
		ArrayList<String> seqDescriptors = new ArrayList<String>(seqDescriptorsArray);
		ArrayList<Pairs> combos = new ArrayList<Pairs>();	
		ArrayList<Pairs> treePairs = new ArrayList<Pairs>();
		
		double[][] upgmaTable = new double[sequencesStored.size()][sequencesStored.size()];
		
		int matrixScore = 0;
		int highscore = 0;
		int tempHighscore = 0;
		int location = 0;
		int matchScore;
		int mismatchScore;
		int gapOpen;
		int gapExt;
		int row = 0;
		int column = 0;
		int larger;
		int smaller;
		int score1;
		int score2;
		double rand;
		double lowTreeScore;
		
		
		boolean first;
		boolean agreement;
		char nt;
		
		String seqTitle1 = new String();
		String seqTitle2 = new String();
		String tempSeq1 = new String();
		String tempSeq2 = new String();
		String refSeq1 = new String();
		String refSeq2 = new String();
		String alignedSeq1 = new String();
		String alignedSeq2 = new String();
		String consensus = new String("");
		String UPGMA = new String("");
		
		//Ask the user for parameters
		System.out.println("Enter match score: ");
		matchScore = scanner.nextInt();
		System.out.println("Enter mismatch score: ");
		mismatchScore = scanner.nextInt();
		System.out.println("Enter gap-opening penalty: ");
		gapOpen = scanner.nextInt();
		System.out.println("Enter gap extension penalty: ");
		gapExt = scanner.nextInt();
		
		//Create pairs of sequences (objects) and add them to the combos ArrayList
		for(int i = 0; i < sequences.size()-1; i++){
			for(int j = i+1; j < sequences.size(); j++){
				tempSeq1 = sequences.get(i);
				tempSeq2 = sequences.get(j);
				
				Pairs p = new Pairs(tempSeq1, tempSeq2);
				combos.add(p);
				
				matrixScore = alignSequences(tempSeq1, tempSeq2, matchScore, mismatchScore, gapOpen, gapExt, p);
				
			}
		}
		
		while(sequences.size() > 1){
		
			highscore = 0;
		//Find the highest scoring sequence pair and store the location within combos
		for(int i = 0; i < combos.size(); i++){
			tempHighscore = combos.get(i).getHighscore();
			if(tempHighscore > highscore){
				highscore = tempHighscore;
				location = i;
			}
		}
		
		//Get the two sequences and two aligned sequences from the corresponding Pairs object
		refSeq1 = combos.get(location).getSeq1();
		refSeq2 = combos.get(location).getSeq2();
		alignedSeq1 = combos.get(location).getAlignedSeq1();
		alignedSeq2 = combos.get(location).getAlignedSeq2();
		
		
		//Create a consensus from the two aligned sequences giving preference to non-gaps
		consensus = "";
		for(int i = 0; i < alignedSeq1.length(); i++){
			if(alignedSeq1.charAt(i) == alignedSeq2.charAt(i)){
				consensus += alignedSeq1.charAt(i);
			}
			else if(alignedSeq1.charAt(i) != alignedSeq2.charAt(i)){
				if(alignedSeq1.charAt(i) == '-'){
					consensus += alignedSeq2.charAt(i);
				}
				else if(alignedSeq2.charAt(i) == '-'){
					consensus += alignedSeq1.charAt(i);
				}
				else{
					rand = Math.random();
					if(rand <= 0.5){
						consensus += alignedSeq1.charAt(i);
					}
					else if(rand > 0.5){
						consensus += alignedSeq2.charAt(i);
					}
				}
			}
		}

		
		sequences.remove(refSeq1);
		sequences.remove(refSeq2);
		
		//Remove pairs which contained the two highscoring sequences
		for(int i = 0; i < combos.size(); i++){
			if(combos.get(i).getSeq1() == refSeq1 || combos.get(i).getSeq1() == refSeq2
					|| combos.get(i).getSeq2() == refSeq1 || combos.get(i).getSeq2() == refSeq2){
				combos.remove(i);
			}
		}
		
		//Run sequence pairs with the new consensus sequence.
		for(int i = 0; i < sequences.size(); i++){
			tempSeq1 = sequences.get(i);
			
			Pairs p = new Pairs(tempSeq1, consensus);
			combos.add(p);
			
			matrixScore = alignSequences(tempSeq1, consensus, matchScore, mismatchScore, gapOpen, gapExt, p);

		}
		//Add consensus back to the sequences ArrayList
		sequences.add(consensus);
		
		}//End while
		
		combos.clear();
		
		//Compare all original sequences with the final consensus sequence and store as new objects.
		for(int i = 0; i < sequencesStored.size(); i++){
			tempSeq1 = sequencesStored.get(i);
			
			Pairs p = new Pairs(tempSeq1, consensus);
			combos.add(p);
			
			matrixScore = alignSequences(tempSeq1, consensus, matchScore, mismatchScore, gapOpen, gapExt, p);		
			
		}
		
		for(int i = 0; i < combos.size(); i++){
			System.out.println(combos.get(i).getAlignedSeq1());
		}
		System.out.println("=================================================================");
		System.out.println(consensus);
		
		for(int place = 0; place < consensus.length(); place++){
			agreement = true;
			nt = consensus.charAt(place);
			
			for(int j = 0; j < combos.size(); j++){
				if(combos.get(j).getAlignedSeq1().charAt(place) != nt){
					agreement = false;
					break;
				}
				
			}
			if(agreement == true){
				System.out.print("*");
			} else {
				System.out.print(" ");
			}
		}
		
	//TREE GENERATION
		
		System.out.println();
		System.out.println("UPGMA Tree Generation:");
		
		ArrayList<Integer> clusters = new ArrayList<Integer>();
		int tableSize = sequencesStored.size();
		boolean firstPass = true;
		UPGMA newick = new UPGMA();
		
		//-------------------------------------------------------------------------------------------------------
		//Populate UPGMA table with distance values
		for(int i = 0; i < tableSize; i++){
			for(int j = 0; j <= i; j++){
				if(i == j){
					upgmaTable[i][j] = 0;
				}
				else{
				tempSeq1 = sequencesStored.get(i);
				tempSeq2 = sequencesStored.get(j);
				
				Pairs p = new Pairs(tempSeq1, tempSeq2);
				treePairs.add(p);
				
				matrixScore = alignSequences(tempSeq1, tempSeq2, matchScore, mismatchScore, gapOpen, gapExt, p);
				
				upgmaTable[i][j] = p.getTreeScore();
				}
			}
		}
				
		
		while(tableSize >= 2){
		
		//-------------------------------------------------------------------------------------------------------
		//Find location (x,y) of lowest scoring (non zero) tree
		lowTreeScore = Double.MAX_VALUE;
		for(int i = 0; i < tableSize; i++){
			for(int j = 0; j <= i; j++){
				if(upgmaTable[i][j] < lowTreeScore && i != j){
					lowTreeScore = upgmaTable[i][j];
					row = i;
					column = j;
				}
			}
		}
		larger = Math.max(row, column);
		smaller = Math.min(row, column);
		
		
		//---------------------------------------------------------------------
		//Send sequence title to UPGMA class, if firstPass = True create a new class
		String tempTitle = new String("");
		if( firstPass == true ){
			newick.startTree(seqDescriptors.get(smaller), seqDescriptors.get(larger), lowTreeScore);
			tempTitle = seqDescriptors.get(smaller) + " + " + seqDescriptors.get(larger);
			seqDescriptors.set(smaller, tempTitle);
			seqDescriptors.remove(larger);
			sequencesStored.remove(larger);
			clusters.add(smaller);
			firstPass = false;
		}
		else if( clusters.contains(smaller) ){
			newick.addSequence(seqDescriptors.get(larger), lowTreeScore);
			tempTitle = seqDescriptors.get(smaller) + " + " + seqDescriptors.get(larger);
			seqDescriptors.set(smaller, tempTitle);
			seqDescriptors.remove(larger);
			sequencesStored.remove(larger);
		}
		else if( clusters.contains(larger) && smaller == 0){
			newick.addSequence(seqDescriptors.get(smaller), lowTreeScore);
			
			
		}
		else{
			newick.addCluster(seqDescriptors.get(smaller), seqDescriptors.get(larger), lowTreeScore);
			tempTitle = seqDescriptors.get(smaller) + " + " + seqDescriptors.get(larger);
			seqDescriptors.set(smaller, tempTitle);
			seqDescriptors.remove(larger);
			sequencesStored.remove(larger);
			clusters.add(smaller);
		}
		
		
		//---------------------------------------------------------------------
		//Copy UPGMA Table into a temporary table for reference values
		double[][] tempTable = new double[tableSize][tableSize];
		for(int i = 0; i < tableSize; i++){
			for(int j = 0; j <= i; j++){
				tempTable[i][j] = upgmaTable[i][j];
			}
		}
		
		//"Remove" one element of combined cluster by moving it to the end and reducing table size by 1
		for(int i = 1; i < tableSize; i++){
			for(int j = 0; j <= i; j++){
				if(i == j){
					upgmaTable[i][j] = 0;
				}
				if(i == row && j == column){
					continue;
				}
				else if(i < larger && j < larger && i != j){
					upgmaTable[i][j] = upgmaTable[i][j];
				}
				else if(i > larger && j < larger && i != j){
					upgmaTable[i-1][j] = upgmaTable[i][j];
				}
				else if(i >= larger && j >= larger && i != j){
					upgmaTable[i-1][j-1] = upgmaTable[i][j];

				}
			}
		}
		
		
		//Remove the last column and row from the UPGMA table
		tableSize -= 1;
		//---------------------------------------------------------------------
		//Recalculate table values for combinations with selected sequences
		
		for(int i = 1; i < tableSize; i++){
			for(int j = 0; j < i; j++){
				if(i == smaller){
					upgmaTable[i][j] = (tempTable[i][j] + tempTable[larger][j])/2;
				}
				else if(j == smaller){
					upgmaTable[i][j] = ((tempTable[i+1][j] + tempTable[i+1][larger])/2);
				}
			}
		}
		
		}//End while
		newick.printUPGMA();
	}// End runMSA method

//----------------------------------------------------------------------------------------------------------------------

	public int alignSequences(String seq1, String seq2, int matchScore, int mismatchScore, int gapOpen, int gapExt, Pairs p){
		
		int diagScore;
		int leftScore;
		int vertScore;
		int alignScore;
		int treeScore;
		double rand;

		String aliSeq1 = new String();
		String aliSeq2 = new String();
		
		final int DIAGONAL = 1;
		final int VERTICAL = 2;
		final int LEFT = 3;
		final int VERT_DIAG_TIE = 4;
		final int LEFT_DIAG_TIE = 5;
		final int VERT_LEFT_TIE = 6;
		final int ALL_TIE = 7;
		
		int seq1Length = seq1.length();
		int seq2Length = seq2.length();
		
		int[][] scoringMatrix = new int[seq2Length+1][seq1Length+1];
		int[][] tracebackMatrix = new int[seq2Length+1][seq1Length+1]; 
		
		for(int i = 0; i < seq2Length + 1; i++){
			for(int j = 0; j < seq1Length + 1; j++){
				scoringMatrix[i][j] = 0;
				if(i == 0 && j >= 1){
					tracebackMatrix[i][j] = 3;
				}
				else if(i >= 1 && j == 0){
					tracebackMatrix[i][j] = 2;
				}
				
			}
		}

		
		//Populate the matrix with actual scores
		for(int i = 1; i < seq2Length + 1; i++) {
			for(int j = 1; j < seq1Length + 1; j++){
				
				
				//If match +1 to all score methods
				if(seq1.charAt(j-1) == seq2.charAt(i-1)){
					diagScore = matchScore;
					leftScore = matchScore;
					vertScore = matchScore;
				}
				//Else instantiate scores to 0 (mismatch)
				else{
					diagScore = mismatchScore;
					leftScore = mismatchScore;
					vertScore = mismatchScore;
				}
				//Diagonal score = match/mismatch + score from diagonal
				//Left score = match/mismatch + gap penalty on left
				//Vertical score = match/mismatch + gap penalty from top
				diagScore += (scoringMatrix[i-1][j-1]);
				
				//Compute best score from cell to the left
				if( i == seq2Length){
					leftScore += scoringMatrix[i][j-1];
				}
				else if(tracebackMatrix[i][j-1] == 3 || tracebackMatrix[i][j-1] == 5
						|| tracebackMatrix[i][j-1] == 6 || tracebackMatrix[i][j-1] == 7){
					leftScore += gapExt + scoringMatrix[i][j-1];
				}
				else{
					leftScore += (gapOpen + scoringMatrix[i][j-1]);
				}
				
				//Compute best score from cell above
				if( j == seq1Length){
					vertScore += scoringMatrix[i-1][j];
				}
				else if(tracebackMatrix[i-1][j] == 2 || tracebackMatrix[i-1][j] == 4 
						|| tracebackMatrix[i-1][j] == 6 || tracebackMatrix[i-1][j] == 7){
					vertScore += (gapExt + scoringMatrix[i-1][j]);
				}
				else{
					vertScore += (gapOpen + scoringMatrix[i-1][j]);
				}
				
				//Assigning the highest score of the three possibilities to the matrix
				//At this point we aren't concerned with ties yet
				int highScore = Math.max(diagScore,Math.max(leftScore, vertScore));
				scoringMatrix[i][j] = highScore;
				
				//Determining if any ties are present, assign traceback value to the traceback matrix
				if(diagScore > leftScore && diagScore > vertScore){
					tracebackMatrix[i][j] = DIAGONAL;
				}
				else if(vertScore > diagScore && vertScore > leftScore){
					tracebackMatrix[i][j] = VERTICAL;
				}
				else if(leftScore > diagScore && leftScore > vertScore){
					tracebackMatrix[i][j] = LEFT;
				}
				else if(diagScore == vertScore && diagScore > leftScore){
					tracebackMatrix[i][j] = VERT_DIAG_TIE;
				}
				else if(diagScore == leftScore && diagScore > vertScore){
					tracebackMatrix[i][j] = LEFT_DIAG_TIE;
				}
				else if(leftScore == vertScore && leftScore > diagScore){
					tracebackMatrix[i][j] = VERT_LEFT_TIE;
				}
				else if(diagScore == vertScore && diagScore == leftScore){
					tracebackMatrix[i][j] = ALL_TIE;
				}
			}
		}
		
		alignScore = scoringMatrix[seq2Length][seq1Length];
		
		//-------------------------------------------------------------------------------------------------
		//Follow tracebackMatrix to create aligned sequences
		
		int j = seq1Length;
		int i = seq2Length;
		while(i != 0 || j != 0){
				if(tracebackMatrix[i][j] == 1){
					aliSeq1 = seq1.charAt(j-1) + aliSeq1;
					aliSeq2 = seq2.charAt(i-1) + aliSeq2;
					i += -1;
					j += -1;
				}
				else if(tracebackMatrix[i][j] == 2){
					aliSeq1 = "-" + aliSeq1;
					aliSeq2 = seq2.charAt(i-1) + aliSeq2;
					i += -1;
				}
				else if(tracebackMatrix[i][j] == 3){
					aliSeq1 = seq1.charAt(j-1) + aliSeq1;
					aliSeq2 = "-" + aliSeq2;
					j += -1;
				}
				else if(tracebackMatrix[i][j] == 4){
					rand = Math.random();
					if(rand >= 0.5){
						aliSeq1 = seq1.charAt(j-1) + aliSeq1;
						aliSeq2 = seq2.charAt(i-1) + aliSeq2;
						i += -1;
						j += -1;
					}
					else{
						aliSeq1 = "-" + aliSeq1;
						aliSeq2 = seq2.charAt(i-1) + aliSeq2;
						i += -1;
					}
				}
				else if(tracebackMatrix[i][j] == 5){
					rand = Math.random();
					if(rand >= 0.5){
						aliSeq1 = seq1.charAt(j-1) + aliSeq1;
						aliSeq2 = seq2.charAt(i-1) + aliSeq2;
						i += -1;
						j += -1;
					}
					else{
						aliSeq1 = seq1.charAt(j-1) + aliSeq1;
						aliSeq2 = "-" + aliSeq2;
						j += -1;
					}
				}
				else if(tracebackMatrix[i][j] == 6){
					rand = Math.random();
					if(rand >= 0.5){
						aliSeq1 = "-" + aliSeq1;
						aliSeq2 = seq2.charAt(i-1) + aliSeq2;
						i += -1;
					}
					else{
						aliSeq1 = seq1.charAt(j-1) + aliSeq1;
						aliSeq2 = "-" + aliSeq2;
						j += -1;
					}
				}
				else if(tracebackMatrix[i][j] == 7){
					if(i == 1 && j == 1){
						aliSeq1 = seq1.charAt(j-1) + aliSeq1;
						aliSeq2 = seq2.charAt(i-1) + aliSeq2;
						i += -1;
						j += -1;
					}
					else{
						rand = Math.random();
						if(rand <= 0.33){
							aliSeq1 = seq1.charAt(j-1) + aliSeq1;
							aliSeq2 = seq2.charAt(i-1) + aliSeq2;
							i += -1;
							j += -1;
						}
						else if(rand > 0.33 && rand < 0.66){
							aliSeq1 = "-" + aliSeq1;
							aliSeq2 = seq2.charAt(i-1) + aliSeq2;
							i += -1;
						}
						else if(rand >= 0.66){
							aliSeq1 = seq1.charAt(j-1) + aliSeq1;
							aliSeq2 = "-" + aliSeq2;
							j += -1;
						}
					}
				}
			
		}
		
		
		String consensus = "";
		for(int i1 = 0; i1 < aliSeq1.length(); i1++){
			if(aliSeq1.charAt(i1) == aliSeq2.charAt(i1)){
				consensus += aliSeq1.charAt(i1);
			}
			else if(aliSeq1.charAt(i1) != aliSeq2.charAt(i1)){
				if(aliSeq1.charAt(i1) == '-'){
					consensus += aliSeq2.charAt(i1);
				}
				else if(aliSeq2.charAt(i1) == '-'){
					consensus += aliSeq1.charAt(i1);
				}
				else{
					rand = Math.random();
					if(rand <= 0.5){
						consensus += aliSeq1.charAt(i1);
					}
					else if(rand > 0.5){
						consensus += aliSeq2.charAt(i1);
					}
				}
			}
		}
		
		treeScore = 0;
		for(int k = 0; k < aliSeq1.length(); k++){
			if(aliSeq1.charAt(k) != aliSeq2.charAt(k)){
				treeScore += 1;
			}
		}
		p.setConsensus(consensus);
		p.setTreeScore(treeScore);
		p.setHighscore(alignScore);
		p.setAlignedSeq1(aliSeq1);
		p.setAlignedSeq2(aliSeq2);
		
		return alignScore;
		
	}	//End allignSequences method
}	//End Class
