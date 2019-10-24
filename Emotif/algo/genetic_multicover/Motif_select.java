import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.BufferedReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.Map;
import java.util.*;


import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.util.Vector;
import org.moeaframework.util.io.CommentedLineReader;

/**
 * Multiobjective motif selection problem. Need a motif hit file
 */
public class Motif_select implements Problem 
{
	/**
	 * How many objectives to be used 
	 */ 
	private int numObjectives;
	/**
	 * How many constraints
	 */ 
	public int numConstraints;
	/**
	 *
	 */
	private int binaryStringSize; 
	/**
	* Total number of motifs as provided in the input file which has positive hits
	*/
	public int totalNumMotifs;
	/**
	* Total number of seqs that were hit by the motifs as found by FIMO in the input file
	*/
	//public int totalNumSeqsHit; 
	/**
	* List to store number of pos seqs hits
	*/
	public ArrayList<String> posSeqListHits; 
	/**
	* List to store total number of neg seqs hits
	*/
	public ArrayList<String> negSeqListHits;
	
	/**
	* Total number of seqs in the original input fasta file (foreground)
	*/  
	public int totalPosNumSeqsFasta;
	/**
	* Total number of seqs in the original input fasta file (background)
	*/  
	public int totalNegNumSeqsFasta;
	
	/**
	* A map between motif name and list/vector of seqs it hits in the positive data set
	*/
	public Map< String, ArrayList<String> > posMotifMap = new HashMap< String, ArrayList<String> > ();
	
	/**
	* A map between motif name and list/vector of seqs it hits in the negative data set
	*/
	public Map< String, ArrayList<String> > negMotifMap = new HashMap< String, ArrayList<String> > ();
	
	/**
	 * A map between motif name and binary Id
	 */ 
	 public Map< String, Integer > nameIdMap = new HashMap< String, Integer > ();
	 
	 /**
	 * A map between binary ID and motif name
	 */ 
	 public Map< Integer, String > idNameMap = new HashMap< Integer, String> (); 
	
	/**
	* Number of evaluations for the GA algoritm
	*/  
	public int numEvals;
	/**
	* The alpha (penalty value)
	*/  
	public double alpha;
	/**
	* The max neg coverage per motif
	*/  
	public double maxNegCovMotif;
	
	
	// The min pos cov for each motif
	public double minPosCov_threshold;
	
	
	/**
	* The filtering value between the motifs
	*/  
	public double filterCutoff;
	/**
	 * Constructs Motif-select problem instance loaded from
	 * the specified file.
	 * 
	 * @param file the file containing the motif instances
	 * @throws IOException if an I/O error occurred
	 */
	public Motif_select(String motifFileName) throws IOException 
	{
		int i = 1;
	}
	
	
	/**
	 * Second constructor with the name and ID maps 
	 */ 
	public Motif_select(String motifFileName, Map< String, Integer > inNameIdMap, Map< Integer, String > inIdNameMap) throws IOException 
	{
		//set the maps
		//nameIdMap = inNameIdMap;
		//idNameMap = inIdNameMap;
		//readFile(motifFileName);
	}
	/**
	 * Constructor taking an InputData object
	 */ 
	public Motif_select(InputData inData)
	{
		numObjectives = inData.numObjectives;
		System.out.println("num of objectives:" + numObjectives);
		numConstraints = inData.numConstraints;
		System.out.println("num of constraints:" + numConstraints);
		binaryStringSize = inData.totalNumMotifs;
		totalNumMotifs = inData.totalNumMotifs;
		//totalNumSeqsHit = inData.totalNumSeqsHit;
		posSeqListHits = inData.posSeqListHits;
		negSeqListHits = inData.negSeqListHits;
		totalPosNumSeqsFasta = inData.totalPosNumSeqsFasta;
		totalNegNumSeqsFasta = inData.totalNegNumSeqsFasta;
		posMotifMap = inData.posMotifMap;
		negMotifMap = inData.negMotifMap;
		nameIdMap = inData.nameIdMap;
		idNameMap = inData.idNameMap;
		numEvals = inData.numEvals;
		maxNegCovMotif = inData.maxNegCovMotif;
		filterCutoff = inData.filterCutoff;
		
		// new added yichao li 2-29-2015
		minPosCov_threshold = inData.minPosCov_threshold;
		
		
		alpha = inData.alpha;
		System.out.println("posSeqFa:" + totalPosNumSeqsFasta + " negSeqFa:" + totalNegNumSeqsFasta + " posHitSeq:" + posSeqListHits.size() 
		+ " negHitSeq:" + negSeqListHits.size() + " numEvals:" + numEvals + " alpha:" + alpha  + " maxNegCovMotif:" + maxNegCovMotif +
		" filterCutoff:" + filterCutoff + "\n\n");
	}
	
	
	/**
	 * Find the union of a set 
	 * @d the binary vector of motifs. 1 if motif in solution 0 otherwise
	 * @i the id of the motif not to be included in the union 
	 * @posMotifMap map between motif names and list of sequences
	 */ 
	 public HashSet<String> find_union(boolean[] d, int i)
	 {
		 HashSet<String> unionSeqSet = new HashSet<String>();
		 int counter = 0;
		 String motifName = "";
		 //System.out.println("in union i:" + i);
		 //for (Boolean value : d) 
		 //{
			 //System.out.print("\t" + value);
		 //}
		 //System.out.println();
		 ArrayList<String> posSeqList = new ArrayList<String>();
		 //for (Boolean value : d) 
		 //{
			////skip the motif i 
			//if (counter == i){counter ++;continue;}
			//if (value.booleanValue()) 
			//{
				//motifName = idNameMap.get(counter);
				////System.out.println("Mname:" + motifName);
				////get the sequences 
				////for (String seqName : posMotifMap.get(motifName))
				////{
					////if (posSeqList.contains(seqName)){continue;}
					////posSeqList.add(seqName);
					////unionSeqSet.add(seqName);
				////}
				//HashSet<String> tmpSeqSet = new HashSet<String>(posMotifMap.get(motifName));
				//posSeqSet.addAll(posMotifSeqSet);
				//unionSeqSet.addAll(posMotifSeqSet);
				
			//}
			//counter ++;
		  //}
		  
		  //get the union of all the seqs and then just remove the specific motif
		  String selectMotifName = "";
		  for (Boolean value : d) 
		  {
			motifName = idNameMap.get(counter);
			//skip the motif i 
			if (counter == i){counter ++; selectMotifName = motifName; continue;}
			if (value.booleanValue()) 
			{
				HashSet<String> tmpSeqSet = new HashSet<String>(posMotifMap.get(motifName));
				unionSeqSet.addAll(tmpSeqSet);
			}
			counter ++;
		   }
		  
		  //find the seqs in selectMotifName
		  //HashSet<String> selectSeqSet = new HashSet<String>(posMotifMap.get(selectMotifName));
		  //find the difference
		  //unionSeqSet.removeAll(selectSeqSet);
		  
		  //return posSeqList;
		  return unionSeqSet;
	 }
	 
	
	/**
	 * 
	 */ 
	@Override
	public void evaluate(Solution solution) 
	{
		//the binary solution
		boolean[] d = EncodingUtils.getBinary(solution.getVariable(0));
		//System.out.println(Arrays.toString(d));
		//the objective function values
		double posCumCov = 0.0;
		int setSize = 0;
		double negCumCov = 0.0;
		double overlap = 0.0;
		//lists that store the unique seqs added
		//ArrayList<String> posSeqList = new ArrayList<String>();
		//ArrayList<String> negSeqList = new ArrayList<String>();
		//sets that store the unique seqs added per solution
		HashSet<String> posSeqSet = new HashSet<String>();
		HashSet<String> negSeqSet = new HashSet<String>();
		String motifName = "";
		//map between seq Ids and their depth 
		Map< String, Integer > seqDepthMap = new HashMap< String, Integer> (); 
		//the max depth, for now it is hard coded 
		int maxDepth = 3; 
		//System.out.println("num motifs:" + totalNumMotifs);
		double negMotifCov = 0;
		// the threshold (filter) of added sequences
		/*
		 * 
		 */ 
		double newAddedCut = filterCutoff;//0.05 = 5%
		int addedFlag = 0;
		/*
		 */
		double maxMotifNegCov = maxNegCovMotif;//0.2 = 20%
		int negCovFlag = 0;
		/*
		 * Constriant for min pos coverage per single motif in percentage 
		 */
		 
		// double minMotifPosCov = 0.5;
		
		// new added by yichao li to control this parameter outside the program
		double minMotifPosCov = minPosCov_threshold;
		int minMotifPosCovFlag = 0;
		
		//for (Boolean value : d) 
		//{
			//int v = 0;
			//if (value.booleanValue()) {v=1;}
			//if (!value.booleanValue()) {v=0;}
			
			//System.out.print(" d:" + v);
		//}
		//System.out.println();
		//find the pos cum coverage
		for (int i = 0; i < totalNumMotifs; i++) 
		{
			//System.out.println("d is:" + d[i]);
			if (d[i]) 
			{
				//how new added seqs
				//int newAddedSeqs = 0;
				//get the motif name
				motifName = idNameMap.get(i);
				
				HashSet<String> dSeqSet = new HashSet<String>(posMotifMap.get(motifName));
				HashSet<String> posMotifSeqSet = new HashSet<String>(posMotifMap.get(motifName));
				posSeqSet.addAll(posMotifSeqSet);
				//processing the hits for the seqs
				for (String seqName : posMotifMap.get(motifName))
				{
					if (seqDepthMap.containsKey(seqName))
					{
						int temp = seqDepthMap.get(seqName);
						temp++;
						if (temp < maxDepth)
						{
							seqDepthMap.put(seqName, temp);
						}else
						{
							seqDepthMap.put(seqName, maxDepth);
						}
					}else
					{
						seqDepthMap.put(seqName, 1);
					}
				}
				double posMotifCov = (double)posMotifMap.get(motifName).size()/totalPosNumSeqsFasta;
				//the min pos cov per motif
				if (posMotifCov < minMotifPosCov)
				{
					minMotifPosCovFlag = 1;
				}
				//check if the motif exists in the neg hits file or not first
				HashSet<String> negMotifSeqSet = new HashSet<String>();
				if (negMotifMap.containsKey(motifName))
				{
					negMotifSeqSet = new HashSet<String>(negMotifMap.get(motifName));
					negMotifCov = (double)negMotifMap.get(motifName).size()/totalNegNumSeqsFasta;
					if (negMotifCov > maxMotifNegCov)
					{
						negCovFlag = 1;
					}  
				}else
				{
					negMotifCov = 0;
				}
				negSeqSet.addAll(negMotifSeqSet);
				
				//find the union of all the seqs except d[i]
				//System.out.println("i:" + i + " motifName:" + motifName);
				HashSet<String> unionSeqSet = find_union(d, i);
				//find the difference between d[i] and the union
				HashSet<String> diffSet = new HashSet<String>(dSeqSet);
				diffSet.removeAll(unionSeqSet);
				//Percentage of how many new sequences were added
				double newAddedPer = (double)diffSet.size()/totalPosNumSeqsFasta;
				//System.out.println("addedPer:" + newAddedPer);
				if (newAddedPer < newAddedCut)
				{
					addedFlag = 1;
				}
				setSize ++;
			}
		}
		
		//find the value of third objective value; loop thru all sequences 
		int solDepthSum = 0;
		Set<String> keys = seqDepthMap.keySet();
		//check how many seqs covered by this motif solution
		int numSeqsNotHitBySol = totalPosNumSeqsFasta - keys.size();
		//System.out.println("\tNum of seqs Not hit by the sol:" + numSeqsNotHitBySol);
        for(String seqName: keys)
        {
            int diff = maxDepth - seqDepthMap.get(seqName); 
            solDepthSum += diff;
        }
        //add the penalty for seqs which were not hit by the soluiton
        solDepthSum = solDepthSum + maxDepth * numSeqsNotHitBySol;
       
		//using sets
		posCumCov = (double)posSeqSet.size()/totalPosNumSeqsFasta;
		negCumCov = (double)negSeqSet.size()/totalNegNumSeqsFasta;
		
		//find number of uncovered pos seqs and number of covered neg seqs for this solution
		//using arrays
		//int posUncovered = posSeqListHits.size() - posSeqList.size();
		//int negCovered =   negSeqList.size();
		//using sets
		int posUncovered = posSeqListHits.size() - posSeqSet.size();
		int negCovered =   negSeqSet.size();
		//cost using percentage and an alpha penalty value
		double posUncovPer = (double)posUncovered/totalPosNumSeqsFasta;
		double negCovPer = (double)negCovered/totalNegNumSeqsFasta;
		double cost = (1-alpha)*posUncovPer + (alpha)*negCovPer;
		solution.setObjective(0, cost);
		solution.setObjective(1, setSize);
		solution.setObjective(2, solDepthSum);
		
		//add a constraint such that if one of the motifs in the solution did not add enough new sequences it is penalzied
		double g = 0.0;
		
		if(negCovFlag == 1)
		{
			//g[0] = negCumCov;
			g =1000000000;
		}else
		{
			//g[0] = 0;
			g = 0.0;
		}
		double h = 0;
		//the filtering cutoff threshold
		if (addedFlag == 1)
		{
			h = 1000000000;
		}else
		{
			h = 0;
		}
		//check the set size constraint, no zeroes allowed
		double z = 0.0;
		if (setSize == 0)
		{
			z =  1000000000;
		}else
		{
			z = 0;
		}
		double m = 0;
		//check the min pos cov per motif
		if (minMotifPosCovFlag == 1)
		{
			m = 1000000000;
		}else
		{
			m = 0;
		}

		solution.setConstraint(0, z);
		solution.setConstraint(1, g);
		solution.setConstraint(2, h);
		solution.setConstraint(3, m);
	}
	
	/**
	 * 
	 */ 
	@Override
	public String getName() {
		return "Motif Select";
	}

	@Override
	public int getNumberOfConstraints() {
		return numConstraints;
	}

	@Override
	public int getNumberOfObjectives() {
		return numObjectives;
	}

	@Override
	public int getNumberOfVariables() {
		return 1;//since it is a binary string representation
	}

	@Override
	public Solution newSolution() {
		//Solution solution = new Solution(1, numObjectives);//use 1 as num of variables since using a binary string
		Solution solution = new Solution(1, numObjectives, numConstraints );//use 1 as num of variables since using a binary string and num of obj and number of constraints
		solution.setVariable(0, EncodingUtils.newBinary(binaryStringSize));
		return solution;
	}

	@Override
	public void close() {
		//do nothing
	}
}
