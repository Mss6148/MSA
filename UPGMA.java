//Class to store and print a Newick format tree.  
//The first pass through the UPGMA table will always be a cluster.
//Then based on the lowest scores in the table it may be a sequence + cluster or a new cluster forming.
public class UPGMA {
	
	public String upgma = new String("");

	public UPGMA(){
		
	}
	public void startTree(String seq1, String seq2, double distance){
		upgma = "(" + seq1 + " " + distance + " : " + seq2 + " " + distance + ")";
	}
	
	public void addSequence(String seq3, double distance){
		upgma = "(" + upgma + " " + seq3 + " " + distance + ")";
	}
	
	public  void addCluster(String seq4, String seq5, double distance){
		upgma = "(" + upgma + ")" + "(" + seq4 + " " + distance + " : " + seq5 + " " + distance + "))";
	}
	public  void printUPGMA(){
		System.out.println(upgma);
	}
}
