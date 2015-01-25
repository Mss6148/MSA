import java.io.*;
import java.util.*;
import java.lang.*;

import javax.swing.JFrame;

public class SequenceAlignment {

	public static void main(String[] args) throws Exception{
		
		int choice = 0;
		
		String fileName ="";
		Scanner scanner = new Scanner(System.in);
		
		//JFrame psaWindow = new SequenceInterface();
		//psaWindow.setVisible(true);
		
		System.out.println("Please enter your FastA filename: ");
		fileName = scanner.next();
		
		System.out.println("Filename: " + fileName);
		
		System.out.print("Make a selection: " + "\n" +
						 "1: Pairwise sequence alignment with nucleotides" + "\n" +
						 "2: Pairwise sequence alignment with proteins" + "\n" +
						 "3: Multiple sequence alignment with nucleotides" + "\n" +
						 "4: Multiple sequence alignment with proteins" + "\n");
		choice = scanner.nextInt();
		
		switch (choice) {
		case 1: new FileReader(fileName, 1);
				break;
		case 2: new FileReader(fileName, 2);
				break;
        case 3: new FileReader(fileName, 3);
                break;
        case 4: new FileReader(fileName, 4);
        		break;
		}
			
		
		//FastAReader fastReader = new FastAReader(fileName);
		
		
	}
}
