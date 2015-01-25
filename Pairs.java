
public class Pairs {

	int alignScore;
	double treeScore;
	
	String seq1 = new String();
	String seq2 = new String();
	String aliSeq1 = new String();
	String aliSeq2 = new String();
	String consensus = new String();
	

	public Pairs(String tempSeq1, String tempSeq2) {
		
		seq1 = tempSeq1;
		seq2 = tempSeq2;

		
	}

	
	public int getHighscore(){
		return alignScore;
	}
	
	public void setHighscore(int highscore){
		alignScore = highscore;
	}
	
	public String getSeq1(){
		return seq1;
	}
	
	public void setSeq1(String s1){
		seq1 = s1;
	}
	
	public String getSeq2(){
		return seq2;
	}
	
	public void setSeq2(String s2){
		seq2 = s2;
	}
	
	public String getAlignedSeq1(){
		return aliSeq1;
	}
	
	public void setAlignedSeq1(String as1){
		aliSeq1 = as1;
	}
	
	public String getAlignedSeq2(){
		return aliSeq2;
	}
	
	public void setAlignedSeq2(String as2){
		aliSeq2 = as2;
	}
	
	public double getTreeScore(){
		return treeScore;
	}
	
	public void setTreeScore(double score){
		treeScore = score;
	}
	
	public void setConsensus(String conSeq){
		consensus = conSeq;
	}
	public String getConsensus(){
		return consensus;
	}
}
