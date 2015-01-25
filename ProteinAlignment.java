import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;


public class ProteinAlignment {

	//public static void main(String[] args) throws FileNotFoundException{
	
	public ProteinAlignment() {
		
	}
	
	public void proteinAlign(ArrayList<String> sequences, ArrayList<String> descriptors) throws FileNotFoundException{
		
		
		String seq1 = new String();
		String seq2 = new String();
		
		seq1 = sequences.get(0);
		seq2 = sequences.get(1);
		
		File in = new File("PAM100.txt");
		Scanner input = new Scanner(in);
		
		//----------------------------------------------------------------------------
		int[][] PAM = new int[20][20];
		int tempInt;
		
		for(int i = 0; i < 20; i++){
			for(int j = 0; j < 20; j++){
				tempInt = input.nextInt();
				PAM[i][j] = tempInt;
				
			}
		}
		
		//----------------------------------------------------------------------------
		//Create a hashtable of PAM amino-acid locations.
		Hashtable<Character,Integer> aaIndex = new Hashtable<Character,Integer>(){{
			put('A',	0);
			put('R',	1);
			put('N',	2);
			put('D',	3);
			put('C',	4);
			put('Q',    5);
			put('E',	6);
			put('G',    7);
			put('H',    8);
			put('I',    9);
			put('L',    10);
			put('K',    11);
			put('M',    12);
			put('F',    13);
			put('P',    14);
			put('S',    15);
			put('T',    16);
			put('W',    17);
			put('Y',    18);
			put('V',    19);
		}};
		//----------------------------------------------------------------------------

		int subScore;
		int tempIndex1;
		int tempIndex2;
		int gapOpen;
		int gapExt;
		int matchScore;
		int mismatchScore;
		double rand;
		char c1;
		char c2;
		
		//Variables for scores
		int diagScore = 0,
			leftScore = 0,
			vertScore = 0;
		int highScore = 0;
		int highestScore = 0;
		
		//Sequence Strings - Will come through method call
		
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
	
		Scanner scanner = new Scanner(System.in);
		

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
		//-------------------------------------------------------------------------------------------------
		
		//Populate the matrix with actual scores
		for(int i = 1; i < seq2Length + 1; i++) {
			for(int j = 1; j < seq1Length + 1; j++){
				
				c1 = seq1.charAt(j-1);
				c2 = seq2.charAt(i-1);
				
				tempIndex1 = aaIndex.get(c1);
				tempIndex2 = aaIndex.get(c2);
				
				//Determine score for match/mismatch
				subScore = PAM[tempIndex1][tempIndex2];
				
				diagScore = subScore;
				leftScore = subScore;
				vertScore = subScore;
				
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
					leftScore = leftScore + gapOpen + scoringMatrix[i][j-1];
				}
				
				//Compute best score from cell above
				if( j == seq1Length){
					vertScore += scoringMatrix[i-1][j];
				}
				else if(tracebackMatrix[i-1][j] == 2 || tracebackMatrix[i-1][j] == 4 
						|| tracebackMatrix[i-1][j] == 6 || tracebackMatrix[i-1][j] == 7){
					vertScore += gapExt + scoringMatrix[i-1][j];
				}
				else{
					vertScore = vertScore + gapOpen + scoringMatrix[i-1][j];
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
		
		int j = seq1Length;
		int i = seq2Length;
		//-------------------------------------------------------------------------------------------------
		//Follow tracebackMatrix to create aligned sequences
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
	}
}
