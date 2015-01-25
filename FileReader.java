import java.awt.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;


public class FileReader {
	
	File fileName;
	static String tempSequence = new String();
	
	int numSequences = 0;
	int userChoice;
	
	ArrayList<String> seqDescriptors = new ArrayList<String>();
	ArrayList<String> sequences = new ArrayList<String>();
	
	
	public FileReader(String fileName, int choice) throws FileNotFoundException{
		userChoice = choice;
		File in = new File(fileName);
		tempSequence = "";
		
		try {
			Scanner scanner = new Scanner(in);
			
			while(scanner.hasNextLine()){
				String line = scanner.nextLine().trim();
				if(line.charAt(0) == '>'){
					seqDescriptors.add(line);
					if(tempSequence != ""){
						sequences.add(tempSequence);
						tempSequence = "";
					}

				}

				if(line.charAt(0) != '>'){
					tempSequence = tempSequence + line;

				}
				
			}
			sequences.add(tempSequence);
			scanner.close();
			
		}
		finally{
			
			switch (userChoice){
			case 1: 
				PairwiseAlignment psa = new PairwiseAlignment();
				psa.Align(sequences, seqDescriptors);
				break;
				
			case 2:
				ProteinAlignment protAlign = new ProteinAlignment();
				protAlign.proteinAlign(sequences, seqDescriptors);
				break;
			
			case 3: MSA msa = new MSA();
					msa.runMSA(sequences, seqDescriptors);
					break;
			case 4: MSAProtein msap = new MSAProtein();
					msap.runMSA(sequences, seqDescriptors);
					break;
				
			}
			
		}
	}


}
