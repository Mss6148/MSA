import java.util.ArrayList;
import java.util.Scanner;


public class PairwiseAlignment {
	
	//public static void main(String[] args) {
		
		public PairwiseAlignment(){
			
		}	
			public void Align(ArrayList<String> sequences, ArrayList<String> descriptors){
				
				String seq1 = new String();
				String seq2 = new String();
				
				seq1 = sequences.get(0);
				seq2 = sequences.get(1);
			
		
		Scanner scanner = new Scanner(System.in);
		
		int gapOpen;
		int gapExt;
		int matchScore;
		int mismatchScore;
		double rand;
		
		//Variables for scores
		int diagScore = 0,
			leftScore = 0,
			vertScore = 0;
		int highScore = 0;
		int highestScore = 0;
		
		
		String aliSeq1 = new String("");
		String aliSeq2 = new String("");
		
		//Sequence Lengths
		int seq1Length = seq1.length();
		int seq2Length = seq2.length();
		
		
		//Scores for traceback matrix
		final int DIAGONAL = 1;
		final int VERTICAL = 2;
		final int LEFT = 3;
		final int VERT_DIAG_TIE = 4;
		final int LEFT_DIAG_TIE = 5;
		final int VERT_LEFT_TIE = 6;
		final int ALL_TIE = 7;
		
		//Scoring and traceback Matrix
		int[][] scoringMatrix = new int[seq2Length+1][seq1Length+1];
		int[][] tracebackMatrix = new int[seq2Length+1][seq1Length+1]; 	
	
		
		System.out.print("Enter match score: ");
			matchScore = scanner.nextInt();
		System.out.print("Enter mismatch score: ");
			mismatchScore = scanner.nextInt();
		System.out.print("Enter gap-opening penalty: ");
			gapOpen = scanner.nextInt();
		System.out.print("Enter gap extension penalty: ");
			gapExt = scanner.nextInt();
		
		//Initialize the matrix with gap penalty scores
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
				//Diaganol score = match/mismatch + score from diaganol
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
				highScore = Math.max(diagScore,Math.max(leftScore, vertScore));
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
		//-------------------------------------------------------------------------------------------------
		//Error checking* Print complete matrix
		/*for(int i = 0; i < seq2Length + 1; i++){
			for(int j = 0; j < seq1Length + 1; j++){
				System.out.print(scoringMatrix[i][j] + " ");
			}
			System.out.print("\n");
		}
		System.out.print("\n");
		
		//-------------------------------------------------------------------------------------------------
		//Error checking* Print traceback matrix.
		for(int i = 0; i < seq2Length+1; i++){
			for(int j = 0; j < seq1Length+1; j++){
				System.out.print(tracebackMatrix[i][j] + " ");
			}
			System.out.print("\n");
		}
		System.out.print("\n");
		
		*/
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
			
		}	//End while-loop
		
		System.out.println(aliSeq1);
		System.out.println(aliSeq2);
		
		for(int x = 0; x < seq2Length; x++){
			for( int y = 0; y < seq1Length; y++){
				System.out.print(tracebackMatrix[x][y] + " ");
			}
			System.out.println();
		}
		
		
		//-------------------------------------------------------------------------------------------------
	}	//End main method

}	//End class
