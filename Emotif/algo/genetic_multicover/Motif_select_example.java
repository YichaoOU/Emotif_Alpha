/* Feb-23-2015
 * First example to try motif selection using the MOEA framework
 * 
 * 
 * Next Feb-25 move the rading file stuff and info from Motif_select to here
 * 
 * 
 */
import java.io.IOException;
import java.io.InputStream;
import org.moeaframework.Executor;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Solution;
import org.moeaframework.util.Vector;
import java.io.FileInputStream;
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
import java.util.Properties;
import java.io.PrintWriter;
import java.text.*;
import java.text.DecimalFormat;
//To compile and run
//javac -cp ".:/home/rami/Documents/moea/MOEAFramework-2.4/lib/*" Motif_select_example.java 
//notice the . above
//To run
//java -cp ".:/home/rami/Documents/moea/MOEAFramework-2.4/lib/*" Motif_select_example test ./config_1.properties

/**
 * 
 * 
 */ 
public class Motif_select_example
{	
	/**
	 * Starts the example
	 * 
	 * @param args the command line arguments
	 * @throws IOException if an I/O error occurred
	 */
	public static void main(String[] args) throws IOException 
	{
		//check if motif file name provided
		// Check how many arguments were passed in
		if(args.length == 0)
		{
			System.out.println("Provide a config file and a job ID");
			System.exit(0);
		}
		String jobId = args[0];
		System.out.println("Job id:" + jobId);
		String configFile = args[1];
		
		//make a directory for results
		new File(jobId).mkdir();
		//make an array which stores names of files to be copied
		ArrayList<String> fileList = new ArrayList<String>(); 
		
		//make a config file for job parameters
		String confFileName = jobId + "_conf_file";
		fileList.add(confFileName);
		PrintWriter confFile = new PrintWriter(confFileName, "UTF-8");

		//make a results file
		String resFileName  = jobId + "_results.csv";
		fileList.add(resFileName);
		PrintWriter resFile = new PrintWriter(resFileName, "UTF-8");
		
		//read the config file
		Properties prop = new Properties();
		InputStream input = null;
		try 
		{
			input = new FileInputStream(configFile);
			// load a properties file
			prop.load(input);
			//loop thru the properties obj and write it to the conf file
			for(String key : prop.stringPropertyNames()) 
			{
			  String value = prop.getProperty(key);
			  confFile.println(key + "=" + value);
			}
		}catch (IOException ex) 
		{
			ex.printStackTrace();
		}
		// make an object for the data file
		InputData inData = new InputData();
		//get the total number of seqs in the fasta file
		inData.totalPosNumSeqsFasta = inData.getFastaSeqCount(prop.getProperty("posFastaFile"));
		inData.totalNegNumSeqsFasta = inData.getFastaSeqCount(prop.getProperty("negFastaFile"));
		//read the pos hit file
		inData.readHitFile(prop.getProperty("posMotifHitFile"),1);
		//read the neg hit file
		inData.readHitFile(prop.getProperty("negMotifHitFile"),-1);
		
		//set the other parameters for inData
		inData.setNumEvals(Integer.parseInt(prop.getProperty("numEvals")));
		
		inData.setMaxNegCovMotif(Double.parseDouble(prop.getProperty("maxNegCovPerMotif")));
		inData.setFilterCutoff(Double.parseDouble(prop.getProperty("filterCutoff")));
		inData.setminPosCov_threshold(Double.parseDouble(prop.getProperty("minPosCov_threshold")));
		inData.setAlpha(Double.parseDouble(prop.getProperty("penaltyValue")));
		
		
		
		//read the input file and set the object varaibles
		//inData.readFile(motifFileName);
		//set number of objectives
		inData.setNumObjectives(3);
		confFile.println("Number of objectives=" + inData.numObjectives);
		inData.setNumConstraints(4);
		//inData.setNumConstraints(2);
		confFile.println("Number of constraints=" + inData.numConstraints);
		
		
		//call the algthm
		String algthmName = "NSGAII";
		int numEvaluations = inData.numEvals;
		confFile.println("Algorithm=" + algthmName);
		confFile.println("Number of evaluations=" + numEvaluations);
		NondominatedPopulation result = new Executor()
			.withProblemClass(Motif_select.class, inData)
			.withAlgorithm(algthmName)
			.withMaxEvaluations(numEvaluations)
			.run();
		System.out.println("Results size:" + result.size());
		confFile.println("Results size=" + result.size());
		
		
		//map between motif name and list of solutions it occurs in
		Map< String, ArrayList<Integer> > motifSolMap = new HashMap< String, ArrayList<Integer> > ();
		ArrayList<Integer> tmpSolList = new ArrayList<Integer>();
		//file to write the number of solutions each motif occurs in
		String motifSolNumFileName = jobId + "_motif_num_sols.csv";
		fileList.add(motifSolNumFileName);
		PrintWriter motifSolFile = new PrintWriter(motifSolNumFileName, "UTF-8");
		motifSolFile.println("#Number of solutions:" + result.size() + "\n");
		motifSolFile.println("#Motif,num_sol_per,pos_cov,neg_cov,sols");
		//file to write the solution stats
		String solStatFileName = jobId + "_solution_stats.csv";
		fileList.add(solStatFileName);
		PrintWriter solStatFile = new PrintWriter(solStatFileName, "UTF-8");
		solStatFile.println("#Solution,pos_cov,set_size,neg_cov,Cov_cost,Depth_cost");
		int solId = 0;
		// print the results
		for (int i = 0; i < result.size(); i++) 
		{
			Solution solution = result.get(i);
			double[] objectives = solution.getObjectives();
					
			// negate objectives to return them to their maximized form
			//objectives = Vector.negate(objectives);	
			System.out.println("\nSolution " + (i+1) + ":");
			//System.out.println("Obj 1: " + -1*(objectives[0]));
			System.out.println("Obj 1: " + objectives[0]);
			System.out.println("Obj 2: " + objectives[1]);
			System.out.println("Obj 3: " + objectives[2]);
			System.out.println("Binary String: " + solution.getVariable(0));
			String binaryString = solution.getVariable(0).toString();
			//check the motifs of the binary string
			int motifId = 0;
			ArrayList<String> posSeqList = new ArrayList<String>();
			ArrayList<String> negSeqList = new ArrayList<String>();
			//loop thru the binary string
			double posSetCumCov = 0;
			double negSetCumCov = 0;
			int setSize = 0;
			resFile.println("\nSolution:" + (i+1) );
			resFile.println("#Motif,pos_cov,pos_cum_cov,neg_cov,neg_cum_cov, sol_list");
			solId = i+1;
			double negMotifCov = 0;
			for (int j = 0; j < binaryString.length(); j++)
			{
				char c = binaryString.charAt(j); 
				int digit = Character.getNumericValue(c);
				if (digit == 1)
				{
					String motifName = inData.idNameMap.get(motifId);
					double posMotifCov = (double)inData.posMotifMap.get(motifName).size()/inData.totalPosNumSeqsFasta;
					for (String seqName : inData.posMotifMap.get(motifName))
					{
						if (posSeqList.contains(seqName)){continue;}
						posSeqList.add(seqName);
					}
					
					if (inData.negMotifMap.containsKey(motifName))
					{
						negMotifCov = (double)inData.negMotifMap.get(motifName).size()/inData.totalNegNumSeqsFasta;
						for (String seqName : inData.negMotifMap.get(motifName))
						{
							if (negSeqList.contains(seqName)){continue;}
							negSeqList.add(seqName);
						}
					}else
					{
						negMotifCov = 0;
					}
					
					setSize ++;
					
					double posCumCov = (double)posSeqList.size()/inData.totalPosNumSeqsFasta;
					double negCumCov = (double)negSeqList.size()/inData.totalNegNumSeqsFasta; 
					//System.out.println("motifName:" + motifName + " id:" + motifId + " posCov:" + posMotifCov + " posCumCov:" + posCumCov +
					//" negCov:" + negMotifCov + " negCumCov:" + negCumCov + " tmpPosCumCov:" + tmpPosCumCov);
					resFile.println(motifName + "," + posMotifCov + "," + posCumCov + "," + negMotifCov + "," + negCumCov);
					if (! motifSolMap.containsKey(motifName))
					{
						motifSolMap.put(motifName, new ArrayList<Integer>());
					}
					motifSolMap.get(motifName).add(solId);
					//tmpSolList = motifSolMap.get(motifName);
					//tmpSolList.add(solId);
				}
				motifId ++;       
			}
			//print stats for this selection
			posSetCumCov = (double)posSeqList.size()/inData.totalPosNumSeqsFasta;
			double tmpPosCov = (double)posSeqList.size()/inData.posSeqListHits.size();
			negSetCumCov = (double)negSeqList.size()/inData.totalNegNumSeqsFasta; 
			double tmpNegCov = (double)negSeqList.size()/inData.negSeqListHits.size();
			System.out.println("posCumCov:" + posSetCumCov + " negCumCov:" + negSetCumCov + " tmpPos:" + tmpPosCov + " tmpneg:" + tmpNegCov);
			double ratio = (posSetCumCov/setSize);
			int posUncovered = inData.posSeqListHits.size() - posSeqList.size();
			int negCovered =   negSeqList.size();
			DecimalFormat df = new DecimalFormat("###.#");
			DecimalFormat df_2 = new DecimalFormat(".##");
			posSetCumCov = 100*posSetCumCov;
			negSetCumCov = 100*negSetCumCov;
			solStatFile.println(solId + "," + df.format(posSetCumCov) + "," + setSize + "," + df.format(negSetCumCov) + "," + df_2.format(objectives[0]) + "," +  objectives[2]);
			//if(i == 10){break;}
		}
		
		
		//loop thru the motifSol map and write to file
		for(String motifName: motifSolMap.keySet())
		{
			ArrayList<Integer> solList = motifSolMap.get(motifName);
			//System.out.println(motifName);
			String listString = "";
			for (Integer sol : solList) 
			{
				//System.out.println("\t" + sol);
				listString += String.valueOf(sol) + "-";
			}
			//System.out.println(listString);
			double solPer = (double)solList.size()/result.size();
			DecimalFormat df = new DecimalFormat(".##");
			double posMotifCov = (double)inData.posMotifMap.get(motifName).size()/inData.totalPosNumSeqsFasta;
			double negMotifCov = 0;
			if (inData.negMotifMap.containsKey(motifName))
			{
				negMotifCov = (double)inData.negMotifMap.get(motifName).size()/inData.totalNegNumSeqsFasta;
			}else
			{
				negMotifCov = 0;
			}
			
			motifSolFile.println(motifName + "," + df.format(solPer) + "," + posMotifCov + "," + negMotifCov + "," + listString);
		}
		
		//close the files
		solStatFile.close();
		confFile.close();
		resFile.close();
		motifSolFile.close();
		//move files to results folder
		for (String fileName : fileList) 
		{
			File mvfile =new File(fileName);
    	    mvfile.renameTo(new File(jobId + "/" + mvfile.getName()));
		} 
	}//end of main
	

}


/**
 *class for representing the pos and neg file hits after reading it 
 */
class InputData
{
	/**
	 * How many objectives to be used 
	 */ 
	public int numObjectives;
	/**
	 * How many constraints
	 */ 
	public int numConstraints;
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
	public ArrayList<String> posSeqListHits = new ArrayList<String>(); 
	/**
	* List to store total number of neg seqs hits
	*/
	public ArrayList<String> negSeqListHits = new ArrayList<String>();
	
	/**
	* Total number of seqs in the original input fasta file (foreground)
	*/  
	public int totalPosNumSeqsFasta;
	/**
	* Total number of seqs in the original input fasta file (background)
	*/  
	public int totalNegNumSeqsFasta;
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
	/**
	* The filtering value between the motifs
	*/  
	public double filterCutoff;
	
	
	// The min pos cov for each motif
	public double minPosCov_threshold;
	
	
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
	 * Empty constructor
	 */ 
	public InputData()
	{
	}
	
	//Methods
	/**
	 * Set number of objectives
	 */ 
	 public void setNumObjectives(int k)
	 {
		numObjectives = k;
	 }
	 /**
	  * Set number of constraints
	  */ 
	 public void setNumConstraints(int k)
	 {
		 numConstraints = k;
	 } 
	 /**
	  * Set number of evaluations
	  */ 
	 public void setNumEvals(int k)
	 {
		 numEvals = k;
	 } 
	 /**
	  * Set alpha value
	  */ 
	 public void setMaxNegCovMotif(double k)
	 {
		 maxNegCovMotif = k;
	 } 
	 
	 /**
	  * Set alpha value
	  */ 
	 public void setFilterCutoff(double k)
	 {
		 filterCutoff = k;
	 } 
	 
	 // new added yichao li 2-29-2016
	 public void setminPosCov_threshold(double k)
	 {
		 minPosCov_threshold = k;
	 } 
	 
	 /**
	  * Set alpha value
	  */ 
	 public void setAlpha(double k)
	 {
		 alpha = k;
	 } 
	 /**
	  * Read the motif hit file and return a map between the motif name and the seqs it hits
	  * @filrType whether we are reading a pos or neg motif hit file, 1 is pos and -1 is neg
	  */
	 public void readHitFile(String motifFileName, int fileType) throws IOException
	 {
		 Pattern motifSignal = Pattern.compile(">(.*)");
		 //the total number of motifs as found by in the file
		 int numMotifs = 0;
		 //A map between motif name and binary Id
		 Map< String, Integer > inNameIdMap = new HashMap< String, Integer > ();
		 //A map between binary ID and motif name 
		 Map< Integer, String > inIdNameMap = new HashMap< Integer, String> ();
		 //the map between motif name and list of seqs it hits
		 Map< String, ArrayList<String> > motifMap = new HashMap< String, ArrayList<String> > ();
		 String motifId = "";
		 //init the seq list per each motif
		 ArrayList<String> seqList = new ArrayList<String>();
		 //init the total seq list
		 ArrayList<String> allSeqList = new ArrayList<String>();
		 try
		 {
			BufferedReader reader = new BufferedReader(new FileReader(motifFileName));
			String line;
			//loop thru the lines
			int motifCount = 0;
			while ((line = reader.readLine()) != null)
			{
			  //skip empty lines
			  if (line.isEmpty()){continue;}
			  //System.out.println(line);
			  Matcher motifMatcher = motifSignal.matcher(line);
			  if (motifMatcher.find()) 
			  {
				  //add seq lists to the motif
				  if (seqList.size() != 0)
				  {
					  motifMap.put(motifId, new ArrayList<String>(seqList));
					  seqList.clear();
				  }
				  motifId = motifMatcher.group(1);
				  inNameIdMap.put(motifId, motifCount);
				  inIdNameMap.put(motifCount, motifId);
				  motifCount ++;
				  numMotifs ++;
				  continue;
			  }
			  seqList.add(line);
			  if (! (allSeqList.contains(line)) )
			  {
				allSeqList.add(line);
			  }
			}
			 
			//finish reading
			reader.close();
		 }catch (Exception e)
	     {
			System.err.format("Exception occurred trying to read '%s'.", motifFileName);
			e.printStackTrace();
	     }
	     //populate the map for last motif 
		 motifMap.put(motifId, seqList);
		 
		 //fill the object variables
		 if (fileType == 1)
		 {
			totalNumMotifs = numMotifs;
			//copy the seq list
			for(String item: allSeqList) posSeqListHits.add(item); 
			//copy the motif map
			for(String motifName: motifMap.keySet())
			{
				posMotifMap.put(motifName, motifMap.get(motifName));
			}
			//copy the IDs and names
			for(String motifName: inNameIdMap.keySet())
			{
				nameIdMap.put(motifName, inNameIdMap.get(motifName));
			}
			for(int id: inIdNameMap.keySet())
			{
				idNameMap.put(id, inIdNameMap.get(id));
			}
			 
		 }
		 if (fileType == -1)
		 {
			 //copy the seq list
			 for(String item: allSeqList) negSeqListHits.add(item); 
			 //copy the motif map
			 for(String motifName: motifMap.keySet())
			 {
				 negMotifMap.put(motifName, motifMap.get(motifName));
			 }
		 } 
	 } 
	 
	 
	/**
	 * Read the input file which is a map between motif names and the pos and neg seqs it occurs in
	 */
	public void readFile(String motifFileName) throws IOException 
	{
		System.out.println(motifFileName);
		//totalPosNumSeqsFasta = 1000;
		//create the Pattern objects
		Pattern motifSignal = Pattern.compile(">(.*)");
		Pattern posSeqSignal = Pattern.compile("positive_seqs");
		Pattern negSeqSignal = Pattern.compile("negative_seqs");
		
		totalNumMotifs = 0;
		String motifId = "";
		int posFlag = 0;
		int negFlag = 0;
		//init the lists
		posSeqListHits = new ArrayList<String>();
		negSeqListHits = new ArrayList<String>();
		//the tmp array lists
		ArrayList<String> posSeqList = new ArrayList<String>();
		ArrayList<String> negSeqList = new ArrayList<String>();
		try
	    {
			BufferedReader reader = new BufferedReader(new FileReader(motifFileName));
			String line;
			//loop thru the lines
			int motifCount = 0;
			while ((line = reader.readLine()) != null)
			{
			  //skip empty lines
			  if (line.isEmpty()){continue;}
			  //System.out.println(line);
			  Matcher motifMatcher = motifSignal.matcher(line);
			  if (motifMatcher.find()) 
			  {
				  //System.out.println(line);
				  //add the neg and pos seq lists to the motif
				  if (posSeqList.size() != 0)
				  {
					  posMotifMap.put(motifId, new ArrayList<String>(posSeqList));
					  posSeqList.clear();
					  //if no neg hits just init the empty list for it
					  if (negSeqList.size() == 0)
					  {
						  negMotifMap.put(motifId, new ArrayList<String>()); 
					  }
				  }
				  if (negSeqList.size() !=0)
				  {
					  negMotifMap.put(motifId, new ArrayList<String>(negSeqList));
					  negSeqList.clear();
				  }
				  //System.out.println(motifMatcher.group(1));
				  motifId = motifMatcher.group(1);
				  nameIdMap.put(motifId, motifCount);
				  idNameMap.put(motifCount, motifId);
				  motifCount ++;
				  totalNumMotifs++;
				  continue;
			  }
			  Matcher posMatch = posSeqSignal.matcher(line);
			  if (posMatch.find())
			  {
				  //System.out.println(line);
				  posFlag = 1;
				  negFlag = 0;
				  continue;
			  }
			  Matcher negMatch = negSeqSignal.matcher(line);
			  if(negMatch.find())
			  {
				  //System.out.println(line);
				  negFlag = 1;
				  posFlag = 0;
				  continue;
			  }
			  //populate the map
			  if(posFlag == 1)
			  {
				posSeqList.add(line);
				if (! (posSeqListHits.contains(line)) )
				{
					posSeqListHits.add(line);
				}
			  }
			  if(negFlag == 1)
			  {
				negSeqList.add(line);
				if( ! (negSeqListHits.contains(line)) )
				{
					negSeqListHits.add(line);
				}
			  }
			}
			
			//finish reading
			reader.close();
	    }
	    catch (Exception e)
	    {
			System.err.format("Exception occurred trying to read '%s'.", motifFileName);
			e.printStackTrace();
	    }
	    // add the last motif
	    if (posSeqList.size() !=0)
		{
			posMotifMap.put(motifId, posSeqList);
			//posSeqList.clear();
		}
		if (negSeqList.size() !=0)
		{
			negMotifMap.put(motifId, negSeqList);
			//negSeqList.clear();
		}
	    //checking
	    System.out.println("num of motifs:" + totalNumMotifs);
	}
	
	
	/**
	 * Count number of seqs in the fasta file
	 */ 
	 public int getFastaSeqCount(String fastaFileName)
	 {
		Pattern motifSignal = Pattern.compile(">(.*)");
		int totalNumSeqs = 0;
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(fastaFileName));
			String line;
			//loop thru the lines
			int motifCount = 0;
			while ((line = reader.readLine()) != null)
			{
				//skip empty lines
				if (line.isEmpty()){continue;}
				Matcher motifMatcher = motifSignal.matcher(line);
				if (motifMatcher.find()) 
				{
					totalNumSeqs ++;
				}
			}
		}
		catch (Exception e)
		{
			System.err.format("Exception occurred trying to read '%s'.", fastaFileName);
			e.printStackTrace();
		}
		System.out.println("Total num seqs fasta:" + totalNumSeqs);
		//totalPosNumSeqsFasta = totalNumSeqs;
		return totalNumSeqs;
	 }
}


