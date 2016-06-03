package wekaproject;
import java.io.*;
import java.util.*;

import ca.pfv.spmf.algorithms.sequentialpatterns.BIDE_and_prefixspan_with_strings.AlgoPrefixSpan_with_Strings;
import ca.pfv.spmf.input.sequence_database_list_strings.SequenceDatabase;
import dataPreprocessing.SAXTransformation;
import dataPreprocessing.SAXTransformation_Testing;
import getAttribute.GetAttr;
import ruleGeneration.RuleEvaluation;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.NominalPrediction;

import weka.classifiers.functions.LibSVM;
import weka.classifiers.trees.J48;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.converters.CSVLoader;


import transferToSDB.T2SDB;
import weka.core.converters.ArffSaver;
public class wekaTest {
	static HashSet<List<String>> powerSet = new HashSet<List<String>>();
	
	public static BufferedReader readDataFile(String filename) {
		BufferedReader inputReader = null;
 
		try {
			inputReader = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException ex) {
			System.err.println("File not found: " + filename);
		}
 
		return inputReader;
	}
 
	public static Evaluation classify(Classifier model,
			Instances trainingSet, Instances testingSet) throws Exception {
		Evaluation evaluation = new Evaluation(trainingSet);
 
		model.buildClassifier(trainingSet);
		evaluation.evaluateModel(model, testingSet);
 
		return evaluation;
	}
 
	public static double calculateAccuracy(FastVector predictions) {
		double correct = 0;
 
		for (int i = 0; i < predictions.size(); i++) {
			NominalPrediction np = (NominalPrediction) predictions.elementAt(i);
			if (np.predicted() == np.actual()) {
				correct++;
			}
		}
 
		return 100 * correct / predictions.size();
	}
 
	public static void run(int period, List<String> para_list, String preprocessing_path, String output_path) throws Exception {

		BufferedReader datafile = readDataFile(preprocessing_path + "weka_training_" + period + "_" + para_list +".arff");
 
		Instances data = new Instances(datafile);
		//System.out.println(data.numAttributes() - 1);
		data.setClassIndex(data.numAttributes() - 1);
		
		
		int trainSize = (int) Math.round(data.numInstances() * 0.8);
		int testSize = data.numInstances() - trainSize;
		Instances train = new Instances(data, 0, trainSize);
		//System.out.println(train);
		Instances test = new Instances(data, trainSize, testSize);
		
		
		// Use a set of classifiers
		Classifier[] models = { 
				new J48(), // a decision tree			
//				new Logistic(),	
				new LibSVM(),
		};
 
		// Run for each model
		for (int j = 0; j < models.length; j++) {
			    
			//SVM MODULE, SET KERNEL
            if (j == 1) {
            	
            	
            	try {        	
            		 //LINEAR
            	    String options = ( "-K 0" );
            	    String[] optionsArray = options.split( " " );
            	    models[j].setOptions(optionsArray);                        	
    		        Evaluation validation = classify(models[j], train, test);
    		        FastVector predictions = new FastVector();
    		        predictions.appendElements(validation.predictions());
    		        		     
    		        double percentage  = validation.correct()/(double)(validation.incorrect() + validation.correct());
		            if (percentage < 0.8) continue;
    		        
                	File fout = new File(output_path + "svm_liner_"+ period + "_" + para_list +".arff");                	
             	    FileOutputStream fos = new FileOutputStream(fout);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);            	
            	    
            	    
            	
            	   
    		        osw.write(validation.toSummaryString("\nResults:SVM(LINEAR)\n======\n", true)); 
    		        osw.write("\r\n");
    		        osw.write(validation.toClassDetailsString());
    		        //System.out.println(validation.toSummaryString("\nResults:SVM(LINEAR)\n======\n", true));
    		        //System.out.println(validation.toClassDetailsString());
    		        osw.close();
            	}catch (IOException e) {
    	        	System.out.println("[ERROR] I/O Exception.");
    	            e.printStackTrace();  	
    	        }   
    		    
            	try {        	
            		 //POLY
    		        String options = ( "-K 1" );
            	    String[] optionsArray = options.split( " " );    		   
            	    models[j].setOptions(optionsArray);      
            	    Evaluation validation = classify(models[j], train, test);
            	    validation = classify(models[j], train, test);
            	    FastVector predictions = new FastVector();
    		        predictions.appendElements(validation.predictions());
    		        double percentage  = validation.correct()/(double)(validation.incorrect() + validation.correct());
		            if (percentage < 0.8) continue;
            		
                	File fout = new File(output_path + "svm_poly_" + period + "_" + para_list +".arff");                	
             	    FileOutputStream fos = new FileOutputStream(fout);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);   
                    
    		    
    		    
    		        
    		        osw.write(validation.toSummaryString("\nResults:SVM(PLOY)\n======\n", true)); 
    		        osw.write("\r\n");
    		        osw.write(validation.toClassDetailsString());
    		        osw.close();
            	}catch (IOException e) {
    	        	System.out.println("[ERROR] I/O Exception.");
    	            e.printStackTrace();  	
    	        }   
            	
            } else {
                try {        
                	// Collect every group of predictions for current model in a FastVector
			        FastVector predictions = new FastVector();
		            Evaluation validation = classify(models[j], train, test); 
		            predictions.appendElements(validation.predictions());
		            double percentage  = validation.correct()/(double)(validation.incorrect() + validation.correct());
		            if (percentage < 0.8) continue;
		            
                    File fout = new File(output_path + models[j].getClass().getSimpleName() + "_" + period + "_" + para_list +".arff");                	
             	    FileOutputStream fos = new FileOutputStream(fout);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);   			        
 
		    // Uncomment to see the summary for each training-testing pair.
		    //System.out.println(models[j].toString());
			

			// Calculate overall accuracy of current classifier on all splits
			//double accuracy = calculateAccuracy(predictions);
			
			// Print current classifier's name and accuracy in a complicated,
			// but nice-looking way.
			//System.out.println("Accuracy of " + models[j].getClass().getSimpleName() + ": "
			//		+ String.format("%.2f%%", accuracy)
			//		+ "\n---------------------------------");
		            osw.write(validation.toSummaryString("\nResults:"+models[j].getClass().getSimpleName() + "\n======\n", true)); 
    		        osw.write("\r\n");
    		        osw.write(validation.toClassDetailsString());
		            //System.out.println(validation.toSummaryString("\nResults:"+models[j].getClass().getSimpleName() + "\n======\n", true));
		            //System.out.println(validation.toClassDetailsString());
		            osw.close();
                }catch (IOException e) {
    	        	System.out.println("[ERROR] I/O Exception.");
    	            e.printStackTrace();  	
    	        }   		           
		    
            }
		}
	}
	
	private static void buildPowerSet(List<String> list, int count) {
	   
		    powerSet.add(list);
		 
	 	    for(int i=0; i<list.size(); i++)
	 	    {
	 	        List<String> temp = new ArrayList<String>(list);
	 	        temp.remove(i);
	 	        buildPowerSet(temp, temp.size());
	 	    }
	}
	
	static ArrayList<ArrayList<String>> readCSV(String fullpath) throws FileNotFoundException{
        ArrayList<ArrayList<String>> records = new ArrayList<>();
	    File inputFile = new File(fullpath);
	    Scanner scl = new Scanner(inputFile);
	    while(scl.hasNextLine()){
		    ArrayList<String> newRecord = new ArrayList<>();
		    String[] tokens = scl.nextLine().split(",");
		    for(String token : tokens){
			    newRecord.add(token);
		    }
		    records.add(newRecord);
	    }
	    scl.close();		
	    return records; 
    }
	
	static void writeCSV(String path, String filename, ArrayList<ArrayList<String>> records) throws IOException{
		FileWriter outputFW = new FileWriter(path + filename);
		for(int i=0;i<records.size();i++){
			ArrayList<String> record = records.get(i);
			StringBuilder recordSB = new StringBuilder();
			for(String col : record) recordSB.append(col).append(',');
			recordSB.deleteCharAt(recordSB.length()-1);
			outputFW.write(recordSB.toString());
			if(i < records.size()-1) outputFW.write("\r\n");
		}
		outputFW.close();
	}	
	
	
	 public static ArrayList<ArrayList<String>> read_text_weka(String filename) throws FileNotFoundException {
	     ArrayList<ArrayList<String>> result = new ArrayList<>();    	
	     Scanner sc = new Scanner(new File(filename));
	     int i = 1;
	     while(sc.hasNextLine()){
		     String[] tokens = sc.nextLine().split(", ");  
			 ArrayList<String> temp = new ArrayList<>();  
			 if (i == 1) {		    
			     for (String s : tokens) {
			         temp.add(s);
			     }   
			     i--;
			     result.add(temp); 
			 } else {
			     for (String s : tokens) {
			         temp.add(s);
			     }   
			     result.add(temp); 
			 }
		  }
		  return result;			
	}
	 
	public static void main(String[] args) throws Exception {		
		/**參數設定**/		
		int N = 10;
		int Original_Level = 1;
		int Original_Relative = 1;
		int Original_Data = 1;
		int MA_Relative = 1;
		int MA_N = 0;
        int MA_Diff = 1;
		int user_defined_class = 0;
        int minsup = 220;
        //double minconf = 0.94;
        if (args.length < 4) {
		    System.out.println("Please input: (1) data_path  (2) preprocessing_path  (3) output_path  (4) periods"); 	
		}
        
		String data_path = args[0];
		String preprocessing_path = args[1];
		String output_path = args[2];
		//選MA BIAS的週期
		int period = Integer.parseInt(args[3]); 
		
		
		ArrayList<String> parameter = new ArrayList<>();
		parameter.add("B_N_C_" + period);
		parameter.add("B_N_S_" + period);
		parameter.add("B_N_R_" + period);
		parameter.add("B_N_T_" + period);
		if (MA_N == 1) {
		    parameter.add("M_N_C_" + period);
		    parameter.add("M_N_S_" + period);
		    parameter.add("M_N_R_" + period);
		    parameter.add("M_N_T_" + period);
		} 
		if (MA_Diff == 1){		
		    parameter.add("D_N_C_" + period);
		    parameter.add("D_N_S_" + period);
		    parameter.add("D_N_R_" + period);
		    parameter.add("D_N_T_" + period);
		}
		
		if (MA_Relative == 1) {
		    parameter.add("M_R_C_" + period);
		    parameter.add("M_R_S_" + period);
		    parameter.add("M_R_R_" + period);
		    parameter.add("M_R_T_" + period);
		}
		
		
		
		/**擷取類別**/    
		String path = data_path;    	    
	    ArrayList<ArrayList<String>> records = readCSV(path);
	    HashMap<Integer, String> feature_target = new HashMap<>();	    
	    if (user_defined_class == 1) {
	    	feature_target = GetAttr.featureExtraction_target_user_defined(records);
	    } else {
	    	feature_target = GetAttr.featureExtraction_target(records);
	    }  
		
	    
	    
		/**Sequential Pattern Mining**/
	    
	    /*
	    //(1)對原始之料使用SAX
		//離散化 SAX (Training Data)
		SAXTransformation.start("SAXTransformation_config_petro_subset1_2010.txt");
		
		
		
	    //轉成sequence (Training Data)
		String path_after_discrete_train = "petro_subset1_2010_rate_after_sax_training.csv";
		T2SDB t = new T2SDB();
		int  SDB_Training_Size = t.translate_training_sliding_window(N, path_after_discrete_train,  feature_target, "SDB(Training).txt");
		//System.out.println("Train size " + SDB_Training_Size);
		
		
		
		
		//離散化 SAX (Testing Data)		
	    SAXTransformation_Testing.start("petro_subset1_breakpoints_2010.txt");
	    
	    //轉成sequence (Testing Data)
	    String path_after_discrete_test = "petro_subset1_2010_rate_after_sax_testing.csv";
	    int  SDB_Testing_Size = t.translate_testing_sliding_window(N, path_after_discrete_test, "SDB(Testing).txt");	
	    //System.out.println("Test size " + SDB_Testing_Size);
	    */
	    
	    
	    
	    //(2)
	    //先取BIAS與MA
		GetAttr.featureExtraction_N("transformed_petro_subset1_feature.csv", records, feature_target, N,  period);	
		
		
		//再對bias值進行sax 
		SAXTransformation.start("SAXTransformation_config_petro_subset1_2010.txt");
		SAXTransformation_Testing.start("petro_subset1_breakpoints_2010.txt");
		System.out.println("Done for SAX!");
		
		//轉成Sequence
		String path_after_discrete = "transformed_petro_subset1_feature_for_sax_training.csv";
		T2SDB t = new T2SDB();
		int SDB_Training_Size = t.translate_training_sliding_window(N, path_after_discrete,  feature_target, "SDB(Training).txt");
		System.out.println("Done for Sequence(Training)!");
		String path_of_testing_file = "transformed_petro_subset1_feature_for_sax_testing.csv";
        int SDB_Testing_Size = t.translate_testing_sliding_window(N, path_of_testing_file, "SDB(Testing).txt");		
        System.out.println("Done for Sequence(Testing)!");
		
        
	    //Sequential Pattern Mining
	    SequenceDatabase sequenceDatabase = new SequenceDatabase();
	    sequenceDatabase.loadFile("SDB(Training).txt");	    
	    AlgoPrefixSpan_with_Strings algo = new AlgoPrefixSpan_with_Strings(); 
	    algo.runAlgorithm(sequenceDatabase, "sequential_patterns.txt", minsup);
//	    algo.printStatistics(sequenceDatabase.size());
	    System.out.println("Done for Mining!");	    
	    
	   
	    
	    //產生Rule
	    //int rule_size = RuleEvaluation.start("RuleEvaluation_config.txt", minconf, minsup, N, SDB_Training_Size);
	    //讀取Sequence
	    ArrayList<ArrayList<ArrayList<String>>> sequences = ReadSDB_for_sequence("sequential_patterns.txt");
//	    for (ArrayList<ArrayList<String>> sequence : sequences) {
//	    	System.out.println(sequence);
//	    }
	    System.out.println("Sequences size: " + sequences.size());
	    
	    
	    int debug = 0;
	    if (debug == 0) {
	    /**產生Sequential Feature*/	    
	    HashMap<Integer, ArrayList<Integer>> SF = GetAttr.sequential_feture(records, sequences, ReadSDB_for_testing("SDB(Testing).txt"), Read_Training_Data("SDB(Training).txt"));	    
	    System.out.println("Done for Rule!");	
	    //for (int index : SF.keySet()) {
	    //	System.out.println(index);
	    //}
	    
	    
	    
	    
	    
	    
	    
	   
	    buildPowerSet(parameter, parameter.size());
	    System.out.println("Down for build powerset");
	    
		for (List<String> para_list : powerSet) {		
			if (para_list.isEmpty()) continue;
			
    	    GetAttr.featureExtraction_weka(Original_Relative, Original_Data, preprocessing_path + "weka_"  + period + "_" + para_list +".csv" , records, feature_target, period, para_list);  
    	    //System.out.println(para_list);
    	    /**Translate To SDB**/
    	    /**1.Training Data**/
    	    
    	    T2SDB t2sdb = new T2SDB();   
    	    
    	    t2sdb.translate_training_sliding_window_weka_including_level_new(N, preprocessing_path + "weka_"  + period + "_" + para_list +".csv", feature_target, preprocessing_path+"weka_training_" + period + "_" + para_list +".txt", Original_Level, records, records.get(0).size()-1, SF);
    	    
    	    try {
                ArrayList<ArrayList<String>> txt_training = read_text_weka(preprocessing_path+"weka_training_" + period + "_" + para_list +".txt");  
                try {
    		        writeCSV("", preprocessing_path + "weka_training_" + period + "_" + para_list +".csv", txt_training);
    		    } catch (IOException e) {
   			        System.out.println("[ERROR] I/O Exception.");
    			    e.printStackTrace();
   		        }  
            } catch (FileNotFoundException e) {
               
            }
    	    
    	    
    	    // load CSV
    	    CSVLoader loader = new CSVLoader();
    	    loader.setSource(new File(preprocessing_path + "weka_training_" + period + "_" + para_list+".csv"));
    	    Instances data1 = loader.getDataSet();
    	    // save ARFF
    	    ArffSaver saver = new ArffSaver();
    	    saver.setInstances(data1);
    	    saver.setFile(new File(preprocessing_path + "weka_training_" + period + "_" + para_list +".arff"));
    	    //saver.setDestination(new File(args[1]));
    	    saver.writeBatch();    	    
    	    run(period, para_list, preprocessing_path, output_path);
            
		}
	    }
		//Clear
		powerSet = new HashSet<List<String>>();
		
	}
	
	private static void sequential_feture(HashMap<ArrayList<ArrayList<String>>, ArrayList<Double>> readRules,
			HashMap<Integer, ArrayList<ArrayList<String>>> readSDB_for_testing,
			HashMap<Integer, ArrayList<ArrayList<String>>> read_Training_Data) {
		// TODO Auto-generated method stub
		
	}

	//讀取SDB(Training).txt
	static HashMap<Integer, ArrayList<ArrayList<String>>> Read_Training_Data(String filename) throws FileNotFoundException{
	    HashMap<Integer, ArrayList<ArrayList<String>>> result = new HashMap<>();
	    int index = 1;        
	    Scanner sc = new Scanner(new File(filename));        
	    while(sc.hasNextLine()) {        
	        ArrayList<ArrayList<String>> itemsets = new ArrayList<>();
	        String[] tokens = sc.nextLine().split(" -1 -2");
	        String[] tokens_next = tokens[0].split(" -1 ");
	        for (String s : tokens_next) {
	            ArrayList<String> itemset = new ArrayList<>();
	            String[] tokens_next_next = s.split(" ");
	            for (String ss : tokens_next_next) {
	                itemset.add(ss);
	            }
	            itemsets.add(itemset);
	        }
	        result.put(index, itemsets);
	        index = index + 1;
	    }            
	    /*
	    //debug
	    for (Integer i : result.keySet()) {
		    System.out.println(i + " " + result.get(i));
		    
		}*/
	    //System.out.println(result.size());
	    sc.close();
	    return result;
	        
	}
	
	
	//讀取SDB(Testing).txt
	static HashMap<Integer, ArrayList<ArrayList<String>>> ReadSDB_for_testing(String filename) throws FileNotFoundException{
        HashMap<Integer, ArrayList<ArrayList<String>>> result = new HashMap<>();
        int index = 1;        
        Scanner sc = new Scanner(new File(filename));        
        while(sc.hasNextLine()) {        
            ArrayList<ArrayList<String>> itemsets = new ArrayList<>();
         
            String[] tokens = sc.nextLine().split(" -1 -2");
            String[] tokens_next = tokens[0].split(" -1 ");
            for (String s : tokens_next) {
                ArrayList<String> itemset = new ArrayList<>();
                String[] tokens_next_next = s.split(" ");
                for (String ss : tokens_next_next) {
                    itemset.add(ss);
                }
                itemsets.add(itemset);
            }
            result.put(index, itemsets);
            index = index + 1;
        }            
        /*
        //debug
        for (Integer i : result.keySet()) {
	        System.out.println(i + " " + result.get(i));
	    
	    }*/
        //System.out.println(result.size());
        sc.close();
        return result;
        
    }
	//讀取Sequence
		static ArrayList<ArrayList<ArrayList<String>>> ReadSDB_for_sequence(String filename) throws FileNotFoundException{
	        ArrayList<ArrayList<ArrayList<String>>> result = new ArrayList<>();     
	        Scanner sc = new Scanner(new File(filename));        
	        while(sc.hasNextLine()) {        
	            ArrayList<ArrayList<String>> itemsets = new ArrayList<>();
	         
	            String[] tokens = sc.nextLine().split("-1  #SUP:");
	            String[] tokens_next = tokens[0].split(" -1 ");
	            for (String s : tokens_next) {
	                ArrayList<String> itemset = new ArrayList<>();
	                String[] tokens_next_next = s.split(" ");
	                for (String ss : tokens_next_next) {
	                    itemset.add(ss);
	                }
	                itemsets.add(itemset);
	            }
	            result.add(itemsets);
	        }            
	        /*
	        //debug
	        for (Integer i : result.keySet()) {
		        System.out.println(i + " " + result.get(i));
		    
		    }*/
	        //System.out.println(result.size());
	        sc.close();
	        return result;
	        
	    }	
		
	//Read rule file
    static HashMap<ArrayList<ArrayList<String>>, ArrayList<Double>> readRules(String filename) throws FileNotFoundException{
	        
		HashMap<ArrayList<ArrayList<String>>, ArrayList<Double>> result = new HashMap<>();
				
		Scanner sc = new Scanner(new File(filename));
		while(sc.hasNextLine()){
		
		    ArrayList<ArrayList<String>> itemsets = new ArrayList<>();
		    ArrayList<Double> list = new ArrayList<>();
			String[] tokens = sc.nextLine().split("\t:\t");
			//For sup, confidence
			String[] number = tokens[1].split(",\t");
			for (String s : number) {
			    double n = Double.parseDouble(s);
			    list.add(n);
			}
			
			//For items
			String[] tokens_next = tokens[0].split(" -> ");
			String[] tokens_next_next = tokens_next[0].split(" -1 ");
			
			//tokens_next[1] : Rise/Down
			ArrayList<String> itemset_next = new ArrayList<>();
			itemset_next.add(tokens_next[1]);
			
			for(String s : tokens_next_next) {
			    String[] tokens_next_next_next =  s.split(" ");
			    ArrayList<String> itemset = new ArrayList<>();   
			    for(String ss : tokens_next_next_next) {
			        itemset.add(ss);    			    
			    }
			    itemsets.add(itemset);
            }
			itemsets.add(itemset_next);		
			result.put(itemsets, list);
			
		}
		/*
		//debug
		for (ArrayList<ArrayList<String>> key : result.keySet()) {
		    System.out.println(key + " " + result.get(key));
		
		}*/
		
		sc.close();
		return result;	
    }	
	
}
