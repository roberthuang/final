package getAttribute;

import java.io.*;
import java.util.*;





public class GetAttr {
	private static HashMap<Integer, Double> temp_sl = new HashMap<>();
	private static HashMap<Integer, Double> temp_ll = new HashMap<>();
	
	public static HashMap<Integer, String> FS_oil(int att_index, ArrayList<ArrayList<String>> records) {
		 HashMap<Integer, String> result = new HashMap<>(); 	    
	     int col = att_index; 
	     int rise_number = 0;
	     int down_number = 0;
	     for (int i = 1; i < records.size(); i++ ) {       
	            if (i == 1) {	                
	                continue;
	            }
	            
	            if (Double.parseDouble(records.get(i).get(col))- Double.parseDouble(records.get(i-1).get(col)) > 0 ) {
	    	    	result.put(i, "R");     
	    	    	rise_number++;
	    	    } else {
	    	    	result.put(i, "D");  
	    	    	down_number++;
	    	    }	
        } 
	    for (int i = 1; i < records.size(); i++) {
	    	if (result.get(i) == null) {
	    		if (rise_number > down_number) {
	    			result.put(i, "R");	    			
	    		} else {
	    			result.put(i, "D");	    
	    		}
	    	}
	    }
	     
	    return result;		
	}
	
	//得到Sequence
	public static ArrayList<String> getsequence(ArrayList<String> sequence) {
		ArrayList<String> result = new ArrayList<>();	
		for (int i = 0; i < sequence.size(); i++) {
			result.add(sequence.get(i));
		}
		return result;
	}
	
	//得到Rule的前項
	public static ArrayList<ArrayList<String>> get_prefix(ArrayList<ArrayList<String>> rule) {
		ArrayList<ArrayList<String>> result = new ArrayList<>();	
		for (int i = 0; i < rule.size()-1; i++) {
			ArrayList<String> temp = new ArrayList<>();
			temp = getsequence(rule.get(i));
			result.add(temp);
		}
	    return result;
	}
	
	public static HashMap<Integer, ArrayList<Integer>> sequential_feture(ArrayList<ArrayList<String>> records, final ArrayList<ArrayList<ArrayList<String>>> sequences, HashMap<Integer, ArrayList<ArrayList<String>>> SDB_for_testing, HashMap<Integer, ArrayList<ArrayList<String>>> SDB_for_training) {
		//刪除Conflict rules				
		/*
		//對規則做排序:1.confidence 2.support 3.length
		Collections.sort(rule_set_before_top_k, new Comparator<ArrayList<ArrayList<String>>>() {
			@Override
			public int compare(ArrayList<ArrayList<String>> rule1, ArrayList<ArrayList<String>> rule2) {
				double conf_1 = rules.get(rule1).get(1);
				double conf_2 = rules.get(rule2).get(1);
				if (conf_1 > conf_2) {
					return -1;
				} else if (conf_1 == conf_2) {
					double sup_1 = rules.get(rule1).get(0);
					double sup_2 = rules.get(rule2).get(0);
					if (sup_1 > sup_2) {
						return -1;
					} else if (sup_1 == sup_2) {
						int len_1 = rules.get(rule1).size()-1;
						int len_2 = rules.get(rule2).size()-1;
						return len_1 - len_2;
					} else {
						return 1;
					}
				} else {
					return 1;
				}
			}	
		});		
		//for (ArrayList<ArrayList<String>> rule : rule_set_before_top_k) {			
		//	System.out.println("Conf: "+ rules.get(rule).get(1));		
		//}
		int count = 0;
		ArrayList<ArrayList<ArrayList<String>>> rule_set = new ArrayList<>();
		for (ArrayList<ArrayList<String>> rule : rule_set_before_top_k) {
			if (count < top_k) {
//				System.out.println(rules.get(rule).get(1));
		        rule_set.add(rule);	
		        count++;
			}
		}*/
		
		
		
		
		
		/*
		for (int i = 0 ; i < rule_set.size(); i++) {
    	    boolean same = false;
    	    for (int j = i+1; j < rule_set.size(); j++) {
    	        ArrayList<ArrayList<String>> temp1 = new ArrayList<>();
    		    for (int k1 = 0; k1 < rule_set.get(i).size()-1; k1++) {
    		        temp1.add(rule_set.get(i).get(k1));
    	        }    
    	        String str1 = rule_set.get(i).get(rule_set.get(i).size()-1).get(0);
    	        ArrayList<ArrayList<String>> temp2 = new ArrayList<>();
    		    for (int k1 = 0; k1 < rule_set.get(j).size()-1; k1++) {
    		        temp2.add(rule_set.get(j).get(k1));
    	        }
    	        String str2 = rule_set.get(j).get(rule_set.get(j).size()-1).get(0);
    	        if ((temp1.equals(temp2)) && (!str1.equals(str2))) {
    	        	//PRINT DUPLICATES
    	        	//int index_1 = rules_all_index.get(rule_set.get(i));  
    	        	//int index_2 = rules_all_index.get(rule_set.get(j));
    	        	//osw.write(index_1 + "    " + index_2 + "\r\n");
    	        	//osw.write(rules.get(rule_set.get(i)).get(1) + "    " + rules.get(rule_set.get(j)).get(1) + "\r\n");
    		        same = true;
    		        rule_set.remove(j--);		    		        	
    		        break;
    	        } 
    	    }    
    	    if (same) {
    	        //System.out.println(i);
    	    	
    	    		rule_set.remove(i--);	 
    	    	
       
    	    }
    	}	*/
		
		
		
		
		
		HashMap<Integer, ArrayList<Integer>> result = new HashMap<>();    	
		int i = 0;		
//		System.out.println("SDB_for_training size: " + SDB_for_training.size());
//		System.out.println("SDB_for_testing size: " + SDB_for_testing.size());
		for (i = 1; i <= SDB_for_training.size(); i++) {
			ArrayList<Integer> match = new ArrayList<>();
			//欲檢查的sequence
			ArrayList<ArrayList<String>> sequence = SDB_for_training.get(i);
		    for (ArrayList<ArrayList<String>> pattern : sequences) {		    		
		    	//System.out.println(rule);
		    	//System.out.println(prefix_of_rule);
		    	//看每個Sequence是否包含了Rule's prefix
		    	int size = 0;
                int current = 0;
                for (int i_1 = 0; i_1 < sequence.size(); i_1++) {                	
                    for (int j = current; j < pattern.size(); j++) {                                         
                        if (sequence.get(i_1).containsAll(pattern.get(j))) {    
                            current = j;
                            current++;
                            size++;
                        }   
                        break;
                    }                                                            
           
                }        
                //有包含Rule's prefix
                if (size == pattern.size()) {
                    match.add(1);            	                      	
                } else {
                	match.add(0);
                }
		    }
		    result.put(i, match);
		}
		int training_data = (int)((records.size()- 1)*0.8);
        for (int j = 1; j <= SDB_for_testing.size(); j++) {
        	ArrayList<Integer> match = new ArrayList<>();
			//欲檢查的sequence
			ArrayList<ArrayList<String>> sequence = SDB_for_testing.get(j);
		    for (ArrayList<ArrayList<String>> pattern : sequences) {
		    	//System.out.println("s:  " + sequence);
		    	//System.out.println("r:  " + prefix_of_rule);
		    	//看每個Sequence是否包含了Rule's prefix
		    	int size = 0;
                int current = 0;
                for (int i_1 = 0; i_1 < sequence.size(); i_1++) {                	
                    for (int j_1 = current; j_1 < pattern.size(); j_1++) {                                         
                        if (sequence.get(i_1).containsAll(pattern.get(j_1))) {    
                            current = j_1;
                            current++;
                            size++;
                        }   
                        break;
                    }                                                            
           
                }        
                //有包含Rule's prefix
                if (size == pattern.size()) {
                	//System.out.println("Yes");
                    match.add(1);            	                      	
                } else {
                	match.add(0);
                }	
		    }
		    result.put(training_data + j, match);
		}
//        System.out.println("feature size: " +result.size());
		return result;
	}
	public static HashMap<Integer, Double> BIAS_N(int length, int att_index, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, Double> result = new HashMap<>();
    	int col = att_index;   
    	for (int i = 1; i < records.size(); i++) {
    		double bias;
    	    if (i <= length-1) {
    	    	
    	    } else {
    	    	double sum_t = 0;
    	    	if (i - length + 1 >= 1) {
    	    		for (int j = i; j >= i-length+1; j--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(j).get(col));
                    } 	    	    		
    	    	}
    	    	sum_t = sum_t / (double)length;
    	    	bias = (Double.parseDouble(records.get(i).get(att_index)) - sum_t)/(double) sum_t;
    	    	result.put(i, bias);  	    	
    	    }    		
    	}
    	
    	double average =  0.0;
    	for (int i = 1; i <= result.size(); i++) {
    		if (result.get(i) != null) {
    			average += result.get(i);
    		}
    	}
    	
    	average /= (double) result.size();
    	for (int i = 1; i < result.size(); i++) {
    		if (result.get(i) == null) {
    			result.put(i, average);
    		}
    	}

    	return result;
    }	
	//length是要選定的週期
	public static HashMap<Integer, String> Move_Average(int length, String att, int att_index, ArrayList<ArrayList<String>> records) {
        HashMap<Integer, String> result = new HashMap<>();    
        //The column of Target
        int col = att_index;                                                                                                                            
        for (int i = 1; i < records.size(); i++ ) {       
            if (i <= length) {
                result.put(i, "MA"+ att.charAt(0) + length + "_1");     
                continue;
            }
            
            double sum_t = 0;
            double sum_t_1 = 0;
            if (i -length + 1 >= 1) {         
                for (int p_1 = i; p_1 >= i-length+1; p_1--) {                
                    sum_t = sum_t + Double.parseDouble(records.get(p_1).get(col));
                } 
                     
                int j = i - 1;
                if (j - length + 1 >=1) {
                    
                    for (int p_2 = j; p_2 >= j-length+1; p_2--) {
                       
                        sum_t_1 = sum_t_1 + Double.parseDouble(records.get(p_2).get(col));
                    }
                }
            }          
            
            //Rise or Down
            double MA = sum_t/length - sum_t_1/length;     
            if (MA >= 0) {
                //System.out.println("i: " + i + " " + MA);
                result.put(i, "MA" + att.charAt(0) + length + "_1");    
            } else {
                //System.out.println("i: " + i + " " + MA);
                result.put(i, "MA" + att.charAt(0) + length + "_0"); 
            }              
        }            
        return result;
    }    	
	//擷取Numerical BIAS
	public static void featureExtraction_N(String output_filename, ArrayList<ArrayList<String>> records, HashMap<Integer, String> feature_target, int window_size, int period_for_MA_BIAS) {		
    	ArrayList<ArrayList<String>> result = new ArrayList<>();
    	HashMap<Integer, Double> BIAS_N_4 = BIAS_N(period_for_MA_BIAS, 4,records);
    	HashMap<Integer, Double> BIAS_N_3= BIAS_N(period_for_MA_BIAS, 3,records);
    	HashMap<Integer, Double> BIAS_N_2 = BIAS_N(period_for_MA_BIAS, 2,records);
    	HashMap<Integer, Double> BIAS_N_1 = BIAS_N(period_for_MA_BIAS, 1,records);    	
    	HashMap<Integer, String> FT_but = feature(4, records);
    	
    	HashMap<Integer, String> MA_4 = Move_Average(period_for_MA_BIAS, records.get(0).get(4), 4, records);
    	HashMap<Integer, String> MA_3 = Move_Average(period_for_MA_BIAS, records.get(0).get(3), 3, records);
    	HashMap<Integer, String> MA_2 = Move_Average(period_for_MA_BIAS, records.get(0).get(2), 2, records);
    	HashMap<Integer, String> MA_1 = Move_Average(period_for_MA_BIAS, records.get(0).get(1), 1, records);
    	
		int training_data_size = (int) ((records.size()-1)*0.8);

		for (int i = 0; i < records.size(); i++) {		
			ArrayList<String> temp = new ArrayList<>();
			//Add time
			temp.add(records.get(i).get(0));
			if(i == 0) {
				temp.add("BIAS_N_1");
				temp.add("BIAS_N_2");
				temp.add("BIAS_N_3");
				temp.add("BIAS_N_4");	
				temp.add("MA_1");
				temp.add("MA_2");
				temp.add("MA_3");
				temp.add("MA_4");				
			} else {
				temp.add(String.valueOf(BIAS_N_1.get(i)));
				temp.add(String.valueOf(BIAS_N_2.get(i)));
				temp.add(String.valueOf(BIAS_N_3.get(i)));
				temp.add(String.valueOf(BIAS_N_4.get(i)));
				temp.add(MA_1.get(i));
				temp.add(MA_2.get(i));
				temp.add(MA_3.get(i));
				temp.add(MA_4.get(i));
			}	
			//temp.add(records.get(i).get(records.get(i).size()-1));	
			result.add(temp);
		}		
		try {
		writeCSV("", output_filename,result);
		} catch (IOException e) {
			System.out.println("[ERROR] I/O Exception.");
			e.printStackTrace();
		}
	}
	
	
	
	
	
	public static HashMap<Integer, String> feature(int att_index, ArrayList<ArrayList<String>> records) {
		 HashMap<Integer, String> result = new HashMap<>(); 	    
	     int col = att_index; 
	     int rise_number = 0;
	     int down_number = 0;
	     for (int i = 1; i < records.size(); i++ ) {       
	            if (i == 1) {	                
	                continue;
	            }
	            
	            if (Double.parseDouble(records.get(i).get(col))- Double.parseDouble(records.get(i-1).get(col)) > 0 ) {
	    	    	result.put(i, "R");     
	    	    	rise_number++;
	    	    } else {
	    	    	result.put(i, "D");  
	    	    	down_number++;
	    	    }	
       } 
	    for (int i = 1; i < records.size(); i++) {
	    	if (result.get(i) == null) {
	    		if (rise_number > down_number) {
	    			result.put(i, "R");	    			
	    		} else {
	    			result.put(i, "D");	    
	    		}
	    	}
	    }
	     
	    return result;		
	}
	public static HashMap<Integer, String> feature_categories(int att_index, ArrayList<ArrayList<String>> records, int rise, int down) {
		 HashMap<Integer, String> result = new HashMap<>(); 	    
	     int col = att_index; 
	     for (int i = 1; i < records.size(); i++ ) {       
	            if (i == 1) {
	                result.put(i, Integer.toString(rise));     
	                continue;
	            }
	            
	            if (Double.parseDouble(records.get(i).get(col))- Double.parseDouble(records.get(i-1).get(col)) > 0 ) {
	    	    	result.put(i, Integer.toString(rise));     
	    	    } else {
	    	    	result.put(i, Integer.toString(down));  
	    	    }	
        }       	        	
	    return result;		
	}
	
	
	
	public static HashMap<Integer, Double> feature2_weka(int att_index, ArrayList<ArrayList<String>> records) {
		 HashMap<Integer, Double> result = new HashMap<>(); 	    
	     int col = att_index; 
	     for (int i = 1; i < records.size(); i++ ) {       
	            if (i == 1) {
	            	
	            } else {	            
	                double min = Double.parseDouble(records.get(i).get(col))- Double.parseDouble(records.get(i-1).get(col));
	                result.put(i, min);
	            }
         }    
	     
	     double average =  0.0;
	    	for (int i = 1; i < result.size(); i++) {
	    		if (result.get(i) != null) {
	    			average += result.get(i);
	    		}
	    	}
	    	
	    	average /= (double) result.size();
	    	for (int i = 1; i < result.size(); i++) {
	    		if (result.get(i)== null) {
	    			result.put(i, average);
	    		}
	    	}
	     return result;		
	}
	
	public static HashMap<Integer, String> feature2(int att_index, ArrayList<ArrayList<String>> records) {
		 HashMap<Integer, String> result = new HashMap<>(); 	    
	     int col = att_index; 
	     for (int i = 1; i < records.size(); i++ ) {       
	            if (i == 1) {
	            	result.put(i, records.get(0).get(col) + "_R");      
	                continue;
	            }
	            
	            if (Double.parseDouble(records.get(i).get(col))- Double.parseDouble(records.get(i-1).get(col)) >= 0 ) {
	    	    	result.put(i,  records.get(0).get(col) + "_R");     
	    	    } else {
	    	    	result.put(i,  records.get(0).get(col) + "_D");  
	    	    }	
       }       	        	
	    return result;		
	}
	
	
	public static HashMap<Integer, String> match_source_target(HashMap<Integer, String> s, HashMap<Integer, String> t, int sou, int tar) {
		HashMap<Integer, String> result = new HashMap<>(); 
	    for (int i = 1;i <= t.size(); i++) {
	        if (s.get(i) == t.get(i)) {
	        	result.put(i, "Same" + "_" + sou + "_" + tar);
	        } else {
	        	result.put(i, "Diff" + "_" + sou + "_" + tar);
	        }
	    }	    
	    return result;
	}
	
	
	public static HashMap<Integer, String> match_source_target_categories(HashMap<Integer, String> s, HashMap<Integer, String> t, int rise, int down) {
		HashMap<Integer, String> result = new HashMap<>(); 
	    for (int i = 1;i <= t.size(); i++) {
	        if (s.get(i) == t.get(i)) {
	        	result.put(i, Integer.toString(rise));
	        } else {
	        	result.put(i, Integer.toString(down));
	        }
	    }	    
	    return result;
	}
	
	public static HashMap<Integer, String> match_source_target_technical(HashMap<Integer, String> s, HashMap<Integer, String> t, int sou, int tar, String str) {
		HashMap<Integer, String> result = new HashMap<>(); 
	    for (int i = 1;i <= t.size(); i++) {
	    	char s1 = s.get(i).charAt(s.get(i).length()-1);
	    	char t1 = t.get(i).charAt(t.get(i).length()-1);
	        if (s1 == t1) {
	        	result.put(i, "Same" + "_" + sou + "_" + tar + "_" + str);
	        } else {
	        	result.put(i, "Diff" + "_" + sou + "_" + tar + "_" + str);
	        }
	    }	    
	    return result;
	}
	
	
	public static HashMap<Integer, String> Move_Average_same(int length1, int length2, String att, int att_index, ArrayList<ArrayList<String>> records) {
		HashMap<Integer, String> result = new HashMap<>();
		int col = att_index;   
		for (int i = 1; i < records.size(); i++ ) {       
	           if (i <= length1) {
	               result.put(i, "MAa"+ att.charAt(0) + length1 + "_1");     
	               continue;
	           }
	           double sum_t1 = 0;
	           double sum_t2= 0;
	           if (i - length1 + 1 >= 1) { 
	        	   for (int p_1 = i; p_1 >= i-length1+1; p_1--) {                
	                    sum_t1 = sum_t1+ Double.parseDouble(records.get(p_1).get(col));
	        	   }
	           }
	           if (i - length2 + 1 >= 1) { 
	        	   for (int p_1 = i; p_1 >= i-length2+1; p_1--) {                
	                    sum_t2 = sum_t2 + Double.parseDouble(records.get(p_1).get(col));
	        	   }
	           }
	           double MA = sum_t1/length1 - sum_t2/length2;
	           if (MA >= 0) {	                
	                result.put(i, "MAa" + att.charAt(0) + length1 + "_1");    
	           } else {	                
	                result.put(i, "MAa" + att.charAt(0) + length1 + "_0"); 
	           }       	           
		}
		return result;
	}
	
	//有差值Moving Average
	public static HashMap<Integer, Double> Move_Average_Diff_Numeric(int period, int att_index, ArrayList<ArrayList<String>> records) {
        HashMap<Integer, Double> result = new HashMap<>();    
        //The column of Target
        int col = att_index;                                                                                                                            
        for (int i = 1; i < records.size(); i++ ) {       
            if (i < period) {
              
            } else {
            
                double sum_t = 0;
                double sum_t_1 = 0;
                if (i - period + 1 >= 1) {         
                    for (int p_1 = i; p_1 >= i-period+1; p_1--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(p_1).get(col));
                    } 
                     
                    int j = i - 1;
                    if (j - period + 1 >=1) {
                    
                        for (int p_2 = j; p_2 >= j-period+1; p_2--) {
                       
                            sum_t_1 = sum_t_1 + Double.parseDouble(records.get(p_2).get(col));
                        }
                    }
                }          
            
                double MA = sum_t/period - sum_t_1/period;     
                result.put(i, MA);
            
            }
        }  
        
        double average =  0.0;
    	for (int i = 1; i < result.size(); i++) {
    		if (result.get(i) != null) {
    			average += result.get(i);
    		}
    	}
    	
    	average /= (double) result.size();
    	for (int i = 1; i < records.size(); i++) {
    		if (result.get(i) == null) {
    			result.put(i, average);
    		}
    	}
        return result;
    }    
	
	public static HashMap<Integer, Double> Move_Average_Numeric(int period, int att_index, ArrayList<ArrayList<String>> records) {
        HashMap<Integer, Double> result = new HashMap<>();                                                                                                                               
        for (int i = 1; i < records.size(); i++ ) {       
            if (i < period) {
               
            } else {
                double sum_t = 0;
                if (i - period + 1 >= 1) {         
                    for (int p_1 = i; p_1 >= i-period+1; p_1--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(p_1).get(att_index));
                    } 
                     
                
                }          
                double MA = sum_t/period;     
                result.put(i, MA);
            }
        }  
        
        double average =  0.0;
    	for (int i = 1; i < result.size(); i++) {
    		if (result.get(i) != null) {
    			average += result.get(i);
    		}
    	}
//    	System.out.println(average);
    	average /= (double) result.size();
    	for (int i = 1; i < records.size(); i++) {
    		if (result.get(i) == null) {
    			result.put(i, average);
    		}
    	}
        return result;
    }    
	
    
	 
    public static HashMap<Integer, String> BIAS(int length, int att_index, double threshold, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, String> result = new HashMap<>();
    	int col = att_index;   
    	int rise_number = 0;
    	int down_number = 0;
    	for (int i = 1; i < records.size(); i++) {
    		double bias;
    	    if (i <= length-1) {
    	    	//result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_1");
    	    } else {
    	    	double sum_t = 0;
    	    	if (i - length + 1 >= 1) {
    	    		for (int j = i; j >= i-length+1; j--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(j).get(col));
                    } 	    	    		
    	    	}
    	    	sum_t = sum_t / (double)length;
    	    	bias = (Double.parseDouble(records.get(i).get(att_index)) - sum_t)/(double) sum_t;
    	    	if (bias >= threshold) {
    	    		result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_1");	
    	    		 rise_number++;
    	    	} else {
    	    		result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_0");	
    	    		down_number++;
    	    	}    	    	
    	    }
    		
    	}
    	for (int i = 1; i < records.size(); i++) {
    		if (result.get(i) == null) {
    			if (rise_number > down_number) {
    				result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_1");
    			} else {
    				result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_0");
    			}
    		}
    	}
    	
    	
    	return result;
    }
    
    public static HashMap<Integer, String> BIAS_categories(int length, int att_index, double threshold, ArrayList<ArrayList<String>> records, int rise, int down) {
    	HashMap<Integer, String> result = new HashMap<>();
    	int col = att_index;   
    	int rise_number = 0;
    	int down_number = 0;
    	for (int i = 1; i < records.size(); i++) {
    		double bias;
    	    if (i <= length-1) {
    	    	//result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_1");
    	    } else {
    	    	double sum_t = 0;
    	    	if (i - length + 1 >= 1) {
    	    		for (int j = i; j >= i-length+1; j--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(j).get(col));
                    } 	    	    		
    	    	}
    	    	sum_t = sum_t / (double)length;
    	    	bias = (Double.parseDouble(records.get(i).get(att_index)) - sum_t)/(double) sum_t;
    	    	if (bias >= threshold) {
    	    		result.put(i, Integer.toString(i));	
    	    		 rise_number++;
    	    	} else {
    	    		result.put(i, "BIAS_" + records.get(0).get(att_index).charAt(0) + "_" + length + "_" + threshold + "_0");	
    	    		down_number++;
    	    	}    	    	
    	    }
    		
    	}
    	for (int i = 1; i < records.size(); i++) {
    		if (result.get(i) == null) {
    			if (rise_number > down_number) {
    				result.put(i, Integer.toString(rise));
    			} else {
    				result.put(i, Integer.toString(down));
    			}
    		}
    	}
    	
    	
    	return result;
    }
    
    public static HashMap<Integer, Double> BIAS_Numeric(int period, int att_index, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, Double> result = new HashMap<>();
    	
    	for (int i = 1; i < records.size(); i++) {
    		double bias;
    	    if (i <= period-1) {
    	    	
    	    } else {
    	    	double sum_t = 0;
    	    	if (i - period + 1 >= 1) {
    	    		for (int j = i; j >= i-period+1; j--) {                
                        sum_t = sum_t + Double.parseDouble(records.get(j).get(att_index));
                    } 	    	    		
    	    	}
    	    	sum_t = sum_t / (double) period;
    	    	bias = (Double.parseDouble(records.get(i).get(att_index)) - sum_t)/(double) sum_t;
    	    	result.put(i, bias);  	    	
    	    }
    		
    	}
    	
    	double average =  0.0;
    	for (int i = 1; i < result.size(); i++) {
    		if (result.get(i) != null) {
    			average += result.get(i);
    		}
    	}
    	
    	average /= (double) result.size();
    	for (int i = 1; i < records.size(); i++) {
    		if (result.get(i) == null) {
    			result.put(i, average);
    		}
    	}
    	return result;
    }
    
    public static HashMap<Integer, Double> origianl_relative(int att_index, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, Double> result = new HashMap<>();
    	double average = 0;
    	for (int i = 1; i < records.size(); i++) {
    	    //Empty
    		if (i == 1)	{
    	    	
    	    } else {
    	    	double pre = Double.parseDouble(records.get(i-1).get(att_index));
    	    	double now = Double.parseDouble(records.get(i).get(att_index));
    	    	
    	    	double relative = (now - pre)/ (double) pre;
    	    	result.put(i, relative);	    	
    	    	average += relative;
    	    }
    		
    	}
    	
    	average /= result.keySet().size();
    	for (int i = 1; i < records.size(); i++) {
    	    if (result.get(i) == null) {
    	    	result.put(i, average);
    	    }    		
    	}
    	return result;
    }
    
    public static HashMap<Integer, Double> MA_relative(HashMap<Integer, Double> MA) {
    	HashMap<Integer, Double> result = new HashMap<>();
    	double average = 0;
    	for (int i = 1; i <= MA.size(); i++) {
    	    //Empty
    		if (i == 1)	{
    	    	
    	    } else {
    	    	double pre = MA.get(i-1);
    	    	double now = MA.get(i);
    	    	
    	    	double relative = (now - pre)/ (double) pre;
    	    	result.put(i, relative);	    	
    	    	average += relative;
    	    }
    		
    	}
    	
    	average /= result.keySet().size();
    	for (int i = 1; i <= MA.size(); i++) {
    	    if (result.get(i) == null) {
    	    	result.put(i, average);
    	    }    		
    	}
    	return result;
    }
    
    
    
    
    //weka
    public static void featureExtraction_weka(int Original_Relative, int Original_Data,String output_filename, ArrayList<ArrayList<String>> records, HashMap<Integer, String> feature_target, int period, List<String> para_list) {		
    	
		ArrayList<ArrayList<String>> result = new ArrayList<>();
		
		//Original Data Relative
		HashMap<Integer, Double> R_C = origianl_relative(1, records);
		HashMap<Integer, Double> R_S = origianl_relative(2, records);
		HashMap<Integer, Double> R_R = origianl_relative(3, records);
		HashMap<Integer, Double> R_T = origianl_relative(4, records);
		
		
		//BIAS
		HashMap<Integer, Double> B_N_C_2 = BIAS_Numeric(2,1, records);
		HashMap<Integer, Double> B_N_S_2 = BIAS_Numeric(2,2, records);
		HashMap<Integer, Double> B_N_R_2 = BIAS_Numeric(2,3, records);
	  	HashMap<Integer, Double> B_N_T_2 = BIAS_Numeric(2,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_3 = BIAS_Numeric(3,1, records);
		HashMap<Integer, Double> B_N_S_3 = BIAS_Numeric(3,2, records);
		HashMap<Integer, Double> B_N_R_3 = BIAS_Numeric(3,3, records);
	  	HashMap<Integer, Double> B_N_T_3 = BIAS_Numeric(3,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_4 = BIAS_Numeric(4,1, records);
		HashMap<Integer, Double> B_N_S_4 = BIAS_Numeric(4,2, records);
		HashMap<Integer, Double> B_N_R_4 = BIAS_Numeric(4,3, records);
	  	HashMap<Integer, Double> B_N_T_4 = BIAS_Numeric(4,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_5 = BIAS_Numeric(5,1, records);
		HashMap<Integer, Double> B_N_S_5 = BIAS_Numeric(5,2, records);
		HashMap<Integer, Double> B_N_R_5 = BIAS_Numeric(5,3, records);
	  	HashMap<Integer, Double> B_N_T_5 = BIAS_Numeric(5,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_6 = BIAS_Numeric(6,1, records);
		HashMap<Integer, Double> B_N_S_6 = BIAS_Numeric(6,2, records);
		HashMap<Integer, Double> B_N_R_6 = BIAS_Numeric(6,3, records);
	  	HashMap<Integer, Double> B_N_T_6 = BIAS_Numeric(6,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_7 = BIAS_Numeric(7,1, records);
		HashMap<Integer, Double> B_N_S_7 = BIAS_Numeric(7,2, records);
		HashMap<Integer, Double> B_N_R_7 = BIAS_Numeric(7,3, records);
	  	HashMap<Integer, Double> B_N_T_7 = BIAS_Numeric(7,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_8 = BIAS_Numeric(8,1, records);
		HashMap<Integer, Double> B_N_S_8 = BIAS_Numeric(8,2, records);
		HashMap<Integer, Double> B_N_R_8 = BIAS_Numeric(8,3, records);
	  	HashMap<Integer, Double> B_N_T_8 = BIAS_Numeric(8,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_9 = BIAS_Numeric(9,1, records);
		HashMap<Integer, Double> B_N_S_9 = BIAS_Numeric(9,2, records);
		HashMap<Integer, Double> B_N_R_9 = BIAS_Numeric(9,3, records);
	  	HashMap<Integer, Double> B_N_T_9 = BIAS_Numeric(9,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_10 = BIAS_Numeric(10,1, records);
		HashMap<Integer, Double> B_N_S_10 = BIAS_Numeric(10,2, records);
		HashMap<Integer, Double> B_N_R_10 = BIAS_Numeric(10,3, records);
	  	HashMap<Integer, Double> B_N_T_10 = BIAS_Numeric(10,4, records);
	  	
	  	
	  	HashMap<Integer, Double> B_N_C_11 = BIAS_Numeric(11,1, records);
		HashMap<Integer, Double> B_N_S_11 = BIAS_Numeric(11,2, records);
		HashMap<Integer, Double> B_N_R_11 = BIAS_Numeric(11,3, records);
	  	HashMap<Integer, Double> B_N_T_11 = BIAS_Numeric(11,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_12 = BIAS_Numeric(12,1, records);
		HashMap<Integer, Double> B_N_S_12 = BIAS_Numeric(12,2, records);
		HashMap<Integer, Double> B_N_R_12 = BIAS_Numeric(12,3, records);
	  	HashMap<Integer, Double> B_N_T_12 = BIAS_Numeric(12,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_13 = BIAS_Numeric(13,1, records);
		HashMap<Integer, Double> B_N_S_13 = BIAS_Numeric(13,2, records);
		HashMap<Integer, Double> B_N_R_13 = BIAS_Numeric(13,3, records);
	  	HashMap<Integer, Double> B_N_T_13 = BIAS_Numeric(13,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_14 = BIAS_Numeric(14,1, records);
		HashMap<Integer, Double> B_N_S_14 = BIAS_Numeric(14,2, records);
		HashMap<Integer, Double> B_N_R_14 = BIAS_Numeric(14,3, records);
	  	HashMap<Integer, Double> B_N_T_14 = BIAS_Numeric(14,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_15 = BIAS_Numeric(15,1, records);
		HashMap<Integer, Double> B_N_S_15 = BIAS_Numeric(15,2, records);
		HashMap<Integer, Double> B_N_R_15 = BIAS_Numeric(15,3, records);
	  	HashMap<Integer, Double> B_N_T_15 = BIAS_Numeric(15,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_16 = BIAS_Numeric(16,1, records);
		HashMap<Integer, Double> B_N_S_16 = BIAS_Numeric(16,2, records);
		HashMap<Integer, Double> B_N_R_16 = BIAS_Numeric(16,3, records);
	  	HashMap<Integer, Double> B_N_T_16 = BIAS_Numeric(16,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_17 = BIAS_Numeric(17,1, records);
		HashMap<Integer, Double> B_N_S_17 = BIAS_Numeric(17,2, records);
		HashMap<Integer, Double> B_N_R_17 = BIAS_Numeric(17,3, records);
	  	HashMap<Integer, Double> B_N_T_17 = BIAS_Numeric(17,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_18 = BIAS_Numeric(18,1, records);
		HashMap<Integer, Double> B_N_S_18 = BIAS_Numeric(18,2, records);
		HashMap<Integer, Double> B_N_R_18 = BIAS_Numeric(18,3, records);
	  	HashMap<Integer, Double> B_N_T_18 = BIAS_Numeric(18,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_19 = BIAS_Numeric(19,1, records);
		HashMap<Integer, Double> B_N_S_19 = BIAS_Numeric(19,2, records);
		HashMap<Integer, Double> B_N_R_19 = BIAS_Numeric(19,3, records);
	  	HashMap<Integer, Double> B_N_T_19 = BIAS_Numeric(19,4, records);
	  	
	  	HashMap<Integer, Double> B_N_C_20 = BIAS_Numeric(20,1, records);
		HashMap<Integer, Double> B_N_S_20 = BIAS_Numeric(20,2, records);
		HashMap<Integer, Double> B_N_R_20 = BIAS_Numeric(20,3, records);
	  	HashMap<Integer, Double> B_N_T_20 = BIAS_Numeric(20,4, records);
	  	
	  	
	  	//MA
	  	HashMap<Integer, Double> M_N_C_2 = Move_Average_Numeric(2,1, records);
	  	HashMap<Integer, Double> M_N_S_2 = Move_Average_Numeric(2,2, records);
	  	HashMap<Integer, Double> M_N_R_2 = Move_Average_Numeric(2,3, records);
	  	HashMap<Integer, Double> M_N_T_2 = Move_Average_Numeric(2,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_3 = Move_Average_Numeric(3,1, records);
	  	HashMap<Integer, Double> M_N_S_3 = Move_Average_Numeric(3,2, records);
	  	HashMap<Integer, Double> M_N_R_3 = Move_Average_Numeric(3,3, records);
	  	HashMap<Integer, Double> M_N_T_3 = Move_Average_Numeric(3,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_4 = Move_Average_Numeric(4,1, records);
	  	HashMap<Integer, Double> M_N_S_4 = Move_Average_Numeric(4,2, records);
	  	HashMap<Integer, Double> M_N_R_4 = Move_Average_Numeric(4,3, records);
	  	HashMap<Integer, Double> M_N_T_4 = Move_Average_Numeric(4,4, records);
	  	
		HashMap<Integer, Double> M_N_C_5 = Move_Average_Numeric(5,1, records);
	  	HashMap<Integer, Double> M_N_S_5 = Move_Average_Numeric(5,2, records);
	  	HashMap<Integer, Double> M_N_R_5 = Move_Average_Numeric(5,3, records);
	  	HashMap<Integer, Double> M_N_T_5 = Move_Average_Numeric(5,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_6 = Move_Average_Numeric(6,1, records);
	  	HashMap<Integer, Double> M_N_S_6 = Move_Average_Numeric(6,2, records);
	  	HashMap<Integer, Double> M_N_R_6 = Move_Average_Numeric(6,3, records);
	  	HashMap<Integer, Double> M_N_T_6 = Move_Average_Numeric(6,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_7 = Move_Average_Numeric(7,1, records);
	  	HashMap<Integer, Double> M_N_S_7 = Move_Average_Numeric(7,2, records);
	  	HashMap<Integer, Double> M_N_R_7 = Move_Average_Numeric(7,3, records);
	  	HashMap<Integer, Double> M_N_T_7 = Move_Average_Numeric(7,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_8 = Move_Average_Numeric(8,1, records);
	  	HashMap<Integer, Double> M_N_S_8 = Move_Average_Numeric(8,2, records);
	  	HashMap<Integer, Double> M_N_R_8 = Move_Average_Numeric(8,3, records);
	  	HashMap<Integer, Double> M_N_T_8 = Move_Average_Numeric(8,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_9 = Move_Average_Numeric(9,1, records);
	  	HashMap<Integer, Double> M_N_S_9 = Move_Average_Numeric(9,2, records);
	  	HashMap<Integer, Double> M_N_R_9 = Move_Average_Numeric(9,3, records);
	  	HashMap<Integer, Double> M_N_T_9 = Move_Average_Numeric(9,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_10 = Move_Average_Numeric(10,1, records);
	  	HashMap<Integer, Double> M_N_S_10 = Move_Average_Numeric(10,2, records);
	  	HashMap<Integer, Double> M_N_R_10 = Move_Average_Numeric(10,3, records);
	  	HashMap<Integer, Double> M_N_T_10 = Move_Average_Numeric(10,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_11 = Move_Average_Numeric(11,1, records);
	  	HashMap<Integer, Double> M_N_S_11 = Move_Average_Numeric(11,2, records);
	  	HashMap<Integer, Double> M_N_R_11 = Move_Average_Numeric(11,3, records);
	  	HashMap<Integer, Double> M_N_T_11 = Move_Average_Numeric(11,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_12 = Move_Average_Numeric(12,1, records);
	  	HashMap<Integer, Double> M_N_S_12 = Move_Average_Numeric(12,2, records);
	  	HashMap<Integer, Double> M_N_R_12 = Move_Average_Numeric(12,3, records);
	  	HashMap<Integer, Double> M_N_T_12 = Move_Average_Numeric(12,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_13 = Move_Average_Numeric(13,1, records);
	  	HashMap<Integer, Double> M_N_S_13 = Move_Average_Numeric(13,2, records);
	  	HashMap<Integer, Double> M_N_R_13 = Move_Average_Numeric(13,3, records);
	  	HashMap<Integer, Double> M_N_T_13 = Move_Average_Numeric(13,4, records);
	  	
		HashMap<Integer, Double> M_N_C_14 = Move_Average_Numeric(14,1, records);
	  	HashMap<Integer, Double> M_N_S_14 = Move_Average_Numeric(14,2, records);
	  	HashMap<Integer, Double> M_N_R_14 = Move_Average_Numeric(14,3, records);
	  	HashMap<Integer, Double> M_N_T_14 = Move_Average_Numeric(14,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_15 = Move_Average_Numeric(15,1, records);
	  	HashMap<Integer, Double> M_N_S_15 = Move_Average_Numeric(15,2, records);
	  	HashMap<Integer, Double> M_N_R_15 = Move_Average_Numeric(15,3, records);
	  	HashMap<Integer, Double> M_N_T_15 = Move_Average_Numeric(15,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_16 = Move_Average_Numeric(16,1, records);
	  	HashMap<Integer, Double> M_N_S_16 = Move_Average_Numeric(16,2, records);
	  	HashMap<Integer, Double> M_N_R_16 = Move_Average_Numeric(16,3, records);
	  	HashMap<Integer, Double> M_N_T_16 = Move_Average_Numeric(16,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_17 = Move_Average_Numeric(17,1, records);
	  	HashMap<Integer, Double> M_N_S_17 = Move_Average_Numeric(17,2, records);
	  	HashMap<Integer, Double> M_N_R_17 = Move_Average_Numeric(17,3, records);
	  	HashMap<Integer, Double> M_N_T_17 = Move_Average_Numeric(17,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_18 = Move_Average_Numeric(18,1, records);
	  	HashMap<Integer, Double> M_N_S_18 = Move_Average_Numeric(18,2, records);
	  	HashMap<Integer, Double> M_N_R_18 = Move_Average_Numeric(18,3, records);
	  	HashMap<Integer, Double> M_N_T_18 = Move_Average_Numeric(18,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_19 = Move_Average_Numeric(19,1, records);
	  	HashMap<Integer, Double> M_N_S_19 = Move_Average_Numeric(19,2, records);
	  	HashMap<Integer, Double> M_N_R_19 = Move_Average_Numeric(19,3, records);
	  	HashMap<Integer, Double> M_N_T_19 = Move_Average_Numeric(19,4, records);
	  	
	  	HashMap<Integer, Double> M_N_C_20 = Move_Average_Numeric(20,1, records);
	  	HashMap<Integer, Double> M_N_S_20 = Move_Average_Numeric(20,2, records);
	  	HashMap<Integer, Double> M_N_R_20 = Move_Average_Numeric(20,3, records);
	  	HashMap<Integer, Double> M_N_T_20 = Move_Average_Numeric(20,4, records);
	  	
	  	//MA's Relative
	  	HashMap<Integer, Double> M_R_C_2 = MA_relative(M_N_C_2);
	  	HashMap<Integer, Double> M_R_S_2 = MA_relative(M_N_S_2);
	  	HashMap<Integer, Double> M_R_R_2 = MA_relative(M_N_R_2);
	  	HashMap<Integer, Double> M_R_T_2 = MA_relative(M_N_T_2);
	  	
		HashMap<Integer, Double> M_R_C_3 = MA_relative(M_N_C_3);
	  	HashMap<Integer, Double> M_R_S_3 = MA_relative(M_N_S_3);
	  	HashMap<Integer, Double> M_R_R_3 = MA_relative(M_N_R_3);
	  	HashMap<Integer, Double> M_R_T_3 = MA_relative(M_N_T_3);
	  	
	  	HashMap<Integer, Double> M_R_C_4 = MA_relative(M_N_C_4);
	  	HashMap<Integer, Double> M_R_S_4 = MA_relative(M_N_S_4);
	  	HashMap<Integer, Double> M_R_R_4 = MA_relative(M_N_R_4);
	  	HashMap<Integer, Double> M_R_T_4 = MA_relative(M_N_T_4);
	  	
	 	HashMap<Integer, Double> M_R_C_5 = MA_relative(M_N_C_5);
	  	HashMap<Integer, Double> M_R_S_5 = MA_relative(M_N_S_5);
	  	HashMap<Integer, Double> M_R_R_5 = MA_relative(M_N_R_5);
	  	HashMap<Integer, Double> M_R_T_5 = MA_relative(M_N_T_5);
	  	
	  	HashMap<Integer, Double> M_R_C_6 = MA_relative(M_N_C_6);
	  	HashMap<Integer, Double> M_R_S_6 = MA_relative(M_N_S_6);
	  	HashMap<Integer, Double> M_R_R_6 = MA_relative(M_N_R_6);
	  	HashMap<Integer, Double> M_R_T_6 = MA_relative(M_N_T_6);
	  	
	  	HashMap<Integer, Double> M_R_C_7 = MA_relative(M_N_C_7);
	  	HashMap<Integer, Double> M_R_S_7 = MA_relative(M_N_S_7);
	  	HashMap<Integer, Double> M_R_R_7 = MA_relative(M_N_R_7);
	  	HashMap<Integer, Double> M_R_T_7 = MA_relative(M_N_T_7);
	  	
	  	HashMap<Integer, Double> M_R_C_8 = MA_relative(M_N_C_8);
	  	HashMap<Integer, Double> M_R_S_8 = MA_relative(M_N_S_8);
	  	HashMap<Integer, Double> M_R_R_8 = MA_relative(M_N_R_8);
	  	HashMap<Integer, Double> M_R_T_8 = MA_relative(M_N_T_8);
	  	
	  	HashMap<Integer, Double> M_R_C_9 = MA_relative(M_N_C_9);
	  	HashMap<Integer, Double> M_R_S_9 = MA_relative(M_N_S_9);
	  	HashMap<Integer, Double> M_R_R_9 = MA_relative(M_N_R_9);
	  	HashMap<Integer, Double> M_R_T_9 = MA_relative(M_N_T_9);
	  	
	  	HashMap<Integer, Double> M_R_C_10 = MA_relative(M_N_C_10);
	  	HashMap<Integer, Double> M_R_S_10 = MA_relative(M_N_S_10);
	  	HashMap<Integer, Double> M_R_R_10 = MA_relative(M_N_R_10);
	  	HashMap<Integer, Double> M_R_T_10 = MA_relative(M_N_T_10);
	  	
	  	HashMap<Integer, Double> M_R_C_11 = MA_relative(M_N_C_11);
	  	HashMap<Integer, Double> M_R_S_11 = MA_relative(M_N_S_11);
	  	HashMap<Integer, Double> M_R_R_11 = MA_relative(M_N_R_11);
	  	HashMap<Integer, Double> M_R_T_11 = MA_relative(M_N_T_11);
	  	
		HashMap<Integer, Double> M_R_C_12 = MA_relative(M_N_C_12);
	  	HashMap<Integer, Double> M_R_S_12 = MA_relative(M_N_S_12);
	  	HashMap<Integer, Double> M_R_R_12 = MA_relative(M_N_R_12);
	  	HashMap<Integer, Double> M_R_T_12 = MA_relative(M_N_T_12);
	  	
	  	HashMap<Integer, Double> M_R_C_13 = MA_relative(M_N_C_13);
	  	HashMap<Integer, Double> M_R_S_13 = MA_relative(M_N_S_13);
	  	HashMap<Integer, Double> M_R_R_13 = MA_relative(M_N_R_13);
	  	HashMap<Integer, Double> M_R_T_13 = MA_relative(M_N_T_13);
	  	
	 	HashMap<Integer, Double> M_R_C_14 = MA_relative(M_N_C_14);
	  	HashMap<Integer, Double> M_R_S_14 = MA_relative(M_N_S_14);
	  	HashMap<Integer, Double> M_R_R_14 = MA_relative(M_N_R_14);
	  	HashMap<Integer, Double> M_R_T_14 = MA_relative(M_N_T_14);
	  	
	  	HashMap<Integer, Double> M_R_C_15 = MA_relative(M_N_C_15);
	  	HashMap<Integer, Double> M_R_S_15 = MA_relative(M_N_S_15);
	  	HashMap<Integer, Double> M_R_R_15 = MA_relative(M_N_R_15);
	  	HashMap<Integer, Double> M_R_T_15 = MA_relative(M_N_T_15);
	  	
	  	HashMap<Integer, Double> M_R_C_16 = MA_relative(M_N_C_16);
	  	HashMap<Integer, Double> M_R_S_16 = MA_relative(M_N_S_16);
	  	HashMap<Integer, Double> M_R_R_16 = MA_relative(M_N_R_16);
	  	HashMap<Integer, Double> M_R_T_16 = MA_relative(M_N_T_16);
	  	
	  	HashMap<Integer, Double> M_R_C_17 = MA_relative(M_N_C_17);
	  	HashMap<Integer, Double> M_R_S_17 = MA_relative(M_N_S_17);
	  	HashMap<Integer, Double> M_R_R_17 = MA_relative(M_N_R_17);
	  	HashMap<Integer, Double> M_R_T_17 = MA_relative(M_N_T_17);
	  	
	  	HashMap<Integer, Double> M_R_C_18 = MA_relative(M_N_C_18);
	  	HashMap<Integer, Double> M_R_S_18 = MA_relative(M_N_S_18);
	  	HashMap<Integer, Double> M_R_R_18 = MA_relative(M_N_R_18);
	  	HashMap<Integer, Double> M_R_T_18 = MA_relative(M_N_T_18);
	  	
	  	HashMap<Integer, Double> M_R_C_19 = MA_relative(M_N_C_19);
	  	HashMap<Integer, Double> M_R_S_19 = MA_relative(M_N_S_19);
	  	HashMap<Integer, Double> M_R_R_19 = MA_relative(M_N_R_19);
	  	HashMap<Integer, Double> M_R_T_19 = MA_relative(M_N_T_19);
	  	
	  	HashMap<Integer, Double> M_R_C_20 = MA_relative(M_N_C_20);
	  	HashMap<Integer, Double> M_R_S_20 = MA_relative(M_N_S_20);
	  	HashMap<Integer, Double> M_R_R_20 = MA_relative(M_N_R_20);
	  	HashMap<Integer, Double> M_R_T_20 = MA_relative(M_N_T_20);
	  	
	    //MA_Diff
	  	HashMap<Integer, Double> D_N_C_2 = Move_Average_Diff_Numeric(2,1, records);
	  	HashMap<Integer, Double> D_N_S_2 = Move_Average_Diff_Numeric(2,2, records);
	  	HashMap<Integer, Double> D_N_R_2 = Move_Average_Diff_Numeric(2,3, records);
	  	HashMap<Integer, Double> D_N_T_2 = Move_Average_Diff_Numeric(2,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_3 = Move_Average_Diff_Numeric(3,1, records);
	  	HashMap<Integer, Double> D_N_S_3 = Move_Average_Diff_Numeric(3,2, records);
	  	HashMap<Integer, Double> D_N_R_3 = Move_Average_Diff_Numeric(3,3, records);
	  	HashMap<Integer, Double> D_N_T_3 = Move_Average_Diff_Numeric(3,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_4 = Move_Average_Diff_Numeric(4,1, records);
	  	HashMap<Integer, Double> D_N_S_4 = Move_Average_Diff_Numeric(4,2, records);
	  	HashMap<Integer, Double> D_N_R_4 = Move_Average_Diff_Numeric(4,3, records);
	  	HashMap<Integer, Double> D_N_T_4 = Move_Average_Diff_Numeric(4,4, records);
	  	
		HashMap<Integer, Double> D_N_C_5 = Move_Average_Diff_Numeric(5,1, records);
	  	HashMap<Integer, Double> D_N_S_5 = Move_Average_Diff_Numeric(5,2, records);
	  	HashMap<Integer, Double> D_N_R_5 = Move_Average_Diff_Numeric(5,3, records);
	  	HashMap<Integer, Double> D_N_T_5 = Move_Average_Diff_Numeric(5,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_6 = Move_Average_Diff_Numeric(6,1, records);
	  	HashMap<Integer, Double> D_N_S_6 = Move_Average_Diff_Numeric(6,2, records);
	  	HashMap<Integer, Double> D_N_R_6 = Move_Average_Diff_Numeric(6,3, records);
	  	HashMap<Integer, Double> D_N_T_6 = Move_Average_Diff_Numeric(6,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_7 = Move_Average_Diff_Numeric(7,1, records);
	  	HashMap<Integer, Double> D_N_S_7 = Move_Average_Diff_Numeric(7,2, records);
	  	HashMap<Integer, Double> D_N_R_7 = Move_Average_Diff_Numeric(7,3, records);
	  	HashMap<Integer, Double> D_N_T_7 = Move_Average_Diff_Numeric(7,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_8 = Move_Average_Diff_Numeric(8,1, records);
	  	HashMap<Integer, Double> D_N_S_8 = Move_Average_Diff_Numeric(8,2, records);
	  	HashMap<Integer, Double> D_N_R_8 = Move_Average_Diff_Numeric(8,3, records);
	  	HashMap<Integer, Double> D_N_T_8 = Move_Average_Diff_Numeric(8,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_9 = Move_Average_Diff_Numeric(9,1, records);
	  	HashMap<Integer, Double> D_N_S_9 = Move_Average_Diff_Numeric(9,2, records);
	  	HashMap<Integer, Double> D_N_R_9 = Move_Average_Diff_Numeric(9,3, records);
	  	HashMap<Integer, Double> D_N_T_9 = Move_Average_Diff_Numeric(9,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_10 = Move_Average_Diff_Numeric(10,1, records);
	  	HashMap<Integer, Double> D_N_S_10 = Move_Average_Diff_Numeric(10,2, records);
	  	HashMap<Integer, Double> D_N_R_10 = Move_Average_Diff_Numeric(10,3, records);
	  	HashMap<Integer, Double> D_N_T_10 = Move_Average_Diff_Numeric(10,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_11 = Move_Average_Diff_Numeric(11,1, records);
	  	HashMap<Integer, Double> D_N_S_11 = Move_Average_Diff_Numeric(11,2, records);
	  	HashMap<Integer, Double> D_N_R_11 = Move_Average_Diff_Numeric(11,3, records);
	  	HashMap<Integer, Double> D_N_T_11 = Move_Average_Diff_Numeric(11,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_12 = Move_Average_Diff_Numeric(12,1, records);
	  	HashMap<Integer, Double> D_N_S_12 = Move_Average_Diff_Numeric(12,2, records);
	  	HashMap<Integer, Double> D_N_R_12 = Move_Average_Diff_Numeric(12,3, records);
	  	HashMap<Integer, Double> D_N_T_12 = Move_Average_Diff_Numeric(12,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_13 = Move_Average_Diff_Numeric(13,1, records);
	  	HashMap<Integer, Double> D_N_S_13 = Move_Average_Diff_Numeric(13,2, records);
	  	HashMap<Integer, Double> D_N_R_13 = Move_Average_Diff_Numeric(13,3, records);
	  	HashMap<Integer, Double> D_N_T_13 = Move_Average_Diff_Numeric(13,4, records);
	  	
		HashMap<Integer, Double> D_N_C_14 = Move_Average_Diff_Numeric(14,1, records);
	  	HashMap<Integer, Double> D_N_S_14 = Move_Average_Diff_Numeric(14,2, records);
	  	HashMap<Integer, Double> D_N_R_14 = Move_Average_Diff_Numeric(14,3, records);
	  	HashMap<Integer, Double> D_N_T_14 = Move_Average_Diff_Numeric(14,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_15 = Move_Average_Diff_Numeric(15,1, records);
	  	HashMap<Integer, Double> D_N_S_15 = Move_Average_Diff_Numeric(15,2, records);
	  	HashMap<Integer, Double> D_N_R_15 = Move_Average_Diff_Numeric(15,3, records);
	  	HashMap<Integer, Double> D_N_T_15 = Move_Average_Diff_Numeric(15,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_16 = Move_Average_Diff_Numeric(16,1, records);
	  	HashMap<Integer, Double> D_N_S_16 = Move_Average_Diff_Numeric(16,2, records);
	  	HashMap<Integer, Double> D_N_R_16 = Move_Average_Diff_Numeric(16,3, records);
	  	HashMap<Integer, Double> D_N_T_16 = Move_Average_Diff_Numeric(16,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_17 = Move_Average_Diff_Numeric(17,1, records);
	  	HashMap<Integer, Double> D_N_S_17 = Move_Average_Diff_Numeric(17,2, records);
	  	HashMap<Integer, Double> D_N_R_17 = Move_Average_Diff_Numeric(17,3, records);
	  	HashMap<Integer, Double> D_N_T_17 = Move_Average_Diff_Numeric(17,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_18 = Move_Average_Diff_Numeric(18,1, records);
	  	HashMap<Integer, Double> D_N_S_18 = Move_Average_Diff_Numeric(18,2, records);
	  	HashMap<Integer, Double> D_N_R_18 = Move_Average_Diff_Numeric(18,3, records);
	  	HashMap<Integer, Double> D_N_T_18 = Move_Average_Diff_Numeric(18,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_19 = Move_Average_Diff_Numeric(19,1, records);
	  	HashMap<Integer, Double> D_N_S_19 = Move_Average_Diff_Numeric(19,2, records);
	  	HashMap<Integer, Double> D_N_R_19 = Move_Average_Diff_Numeric(19,3, records);
	  	HashMap<Integer, Double> D_N_T_19 = Move_Average_Diff_Numeric(19,4, records);
	  	
	  	HashMap<Integer, Double> D_N_C_20 = Move_Average_Diff_Numeric(20,1, records);
	  	HashMap<Integer, Double> D_N_S_20 = Move_Average_Diff_Numeric(20,2, records);
	  	HashMap<Integer, Double> D_N_R_20 = Move_Average_Diff_Numeric(20,3, records);
	  	HashMap<Integer, Double> D_N_T_20 = Move_Average_Diff_Numeric(20,4, records);
	  	
		for (int i = 0; i < records.size(); i++) {		
			ArrayList<String> temp = new ArrayList<>();
			//Add time
			temp.add(records.get(i).get(0));
			if(i == 0) {
			   if (Original_Data == 1) {
			       //Original Data
                   for (int j = 1; j < records.get(0).size(); j++) {
                       temp.add(records.get(i).get(j));
			       }
			   }
			   if (Original_Relative == 1) {
                   //Original Data Relative
                   temp.add("R_C");
                   temp.add("R_S");
                   temp.add("R_R");
                   temp.add("R_T");
			   }
			   
               //MA and BIAS
               for (int k = 0; k < para_list.size();k++) {
            	   temp.add(para_list.get(k));   
               }
                
               temp.add("Target");
			} else {
				if (Original_Data == 1) {
				    //Original Data
				    temp.add(records.get(i).get(1));
				    temp.add(records.get(i).get(2));
				    temp.add(records.get(i).get(3));
				    temp.add(records.get(i).get(4));
				}
				if (Original_Relative == 1) {
				    //Original Data Relative
				    temp.add(String.valueOf(R_C.get(i)));
	                temp.add(String.valueOf(R_S.get(i)));
	                temp.add(String.valueOf(R_R.get(i)));
	                temp.add(String.valueOf(R_T.get(i)));
				}
				
				//MA and BIAS
				for (int k = 0; k < para_list.size(); k++) {
 					switch(para_list.get(k)) {
 					    case "B_N_C_2":
					        temp.add(String.valueOf(B_N_C_2.get(i)));
					    	break;
 					    case "B_N_S_2":
					        temp.add(String.valueOf(B_N_S_2.get(i)));
					    	break;
 					    case "B_N_R_2":
					        temp.add(String.valueOf(B_N_R_2.get(i)));
					    	break;
 					    case "B_N_T_2":
					        temp.add(String.valueOf(B_N_T_2.get(i)));
					    	break;
					    	
 					    case "B_N_C_3":
					        temp.add(String.valueOf(B_N_C_3.get(i)));
					    	break;
					    case "B_N_S_3":
					        temp.add(String.valueOf(B_N_S_3.get(i)));
					    	break;
					    case "B_N_R_3":
					        temp.add(String.valueOf(B_N_R_3.get(i)));
					    	break;
					    case "B_N_T_3":
					        temp.add(String.valueOf(B_N_T_3.get(i)));
					    	break;
					    	
					    	
					    case "B_N_C_4":
					        temp.add(String.valueOf(B_N_C_4.get(i)));
					    	break;
					    case "B_N_S_4":
					        temp.add(String.valueOf(B_N_S_4.get(i)));
					    	break;
					    case "B_N_R_4":
					        temp.add(String.valueOf(B_N_R_4.get(i)));
					    	break;
					    case "B_N_T_4":
					        temp.add(String.valueOf(B_N_T_4.get(i)));
					    	break;
					    	
					    	
					    case "B_N_C_5":
					        temp.add(String.valueOf(B_N_C_5.get(i)));
					    	break;
					    case "B_N_S_5":
					        temp.add(String.valueOf(B_N_S_5.get(i)));
					    	break;
					    case "B_N_R_5":
					        temp.add(String.valueOf(B_N_R_5.get(i)));
					    	break;
					    case "B_N_T_5":
					        temp.add(String.valueOf(B_N_T_5.get(i)));
					    	break;
					    	
					    case "B_N_C_6":
					        temp.add(String.valueOf(B_N_C_6.get(i)));
					    	break;
					    case "B_N_S_6":
					        temp.add(String.valueOf(B_N_S_6.get(i)));
					    	break;
					    case "B_N_R_6":
					        temp.add(String.valueOf(B_N_R_6.get(i)));
					    	break;
					    case "B_N_T_6":
					        temp.add(String.valueOf(B_N_T_6.get(i)));
					    	break;
					    	
					    case "B_N_C_7":
					        temp.add(String.valueOf(B_N_C_7.get(i)));
					    	break;
					    case "B_N_S_7":
					        temp.add(String.valueOf(B_N_S_7.get(i)));
					    	break;
					    case "B_N_R_7":
					        temp.add(String.valueOf(B_N_R_7.get(i)));
					    	break;
					    case "B_N_T_7":
					        temp.add(String.valueOf(B_N_T_7.get(i)));
					    	break;
					    	
					    case "B_N_C_8":
					        temp.add(String.valueOf(B_N_C_8.get(i)));
					    	break;
					    case "B_N_S_8":
					        temp.add(String.valueOf(B_N_S_8.get(i)));
					    	break;
					    case "B_N_R_8":
					        temp.add(String.valueOf(B_N_R_8.get(i)));
					    	break;
					    case "B_N_T_8":
					        temp.add(String.valueOf(B_N_T_8.get(i)));
					    	break;
					    	
					    case "B_N_C_9":
					        temp.add(String.valueOf(B_N_C_9.get(i)));
					    	break;
					    case "B_N_S_9":
					        temp.add(String.valueOf(B_N_S_9.get(i)));
					    	break;
					    case "B_N_R_9":
					        temp.add(String.valueOf(B_N_R_9.get(i)));
					    	break;
					    case "B_N_T_9":
					        temp.add(String.valueOf(B_N_T_9.get(i)));
					    	break;  	
					   
					    case "B_N_C_10":
					        temp.add(String.valueOf(B_N_C_10.get(i)));
					    	break;
					    case "B_N_S_10":
					        temp.add(String.valueOf(B_N_S_10.get(i)));
					    	break;
					    case "B_N_R_10":
					        temp.add(String.valueOf(B_N_R_10.get(i)));
					    	break;
					    case "B_N_T_10":
					        temp.add(String.valueOf(B_N_T_10.get(i)));
					    	break; 
					    	
					    case "B_N_C_11":
					        temp.add(String.valueOf(B_N_C_11.get(i)));
					    	break;
 					    case "B_N_S_11":
					        temp.add(String.valueOf(B_N_S_11.get(i)));
					    	break;
 					    case "B_N_R_11":
					        temp.add(String.valueOf(B_N_R_11.get(i)));
					    	break;
 					    case "B_N_T_11":
					        temp.add(String.valueOf(B_N_T_11.get(i)));
					    	break;
					    	
 					    case "B_N_C_12":
					        temp.add(String.valueOf(B_N_C_12.get(i)));
					    	break;
					    case "B_N_S_12":
					        temp.add(String.valueOf(B_N_S_12.get(i)));
					    	break;
					    case "B_N_R_12":
					        temp.add(String.valueOf(B_N_R_12.get(i)));
					    	break;
					    case "B_N_T_12":
					        temp.add(String.valueOf(B_N_T_12.get(i)));
					    	break;
					    	
					    	
					    case "B_N_C_13":
					        temp.add(String.valueOf(B_N_C_13.get(i)));
					    	break;
					    case "B_N_S_13":
					        temp.add(String.valueOf(B_N_S_13.get(i)));
					    	break;
					    case "B_N_R_13":
					        temp.add(String.valueOf(B_N_R_13.get(i)));
					    	break;
					    case "B_N_T_13":
					        temp.add(String.valueOf(B_N_T_13.get(i)));
					    	break;
					    	
					    	
					    case "B_N_C_14":
					        temp.add(String.valueOf(B_N_C_14.get(i)));
					    	break;
					    case "B_N_S_14":
					        temp.add(String.valueOf(B_N_S_14.get(i)));
					    	break;
					    case "B_N_R_14":
					        temp.add(String.valueOf(B_N_R_14.get(i)));
					    	break;
					    case "B_N_T_14":
					        temp.add(String.valueOf(B_N_T_14.get(i)));
					    	break;
					    	
					    case "B_N_C_15":
					        temp.add(String.valueOf(B_N_C_15.get(i)));
					    	break;
					    case "B_N_S_15":
					        temp.add(String.valueOf(B_N_S_15.get(i)));
					    	break;
					    case "B_N_R_15":
					        temp.add(String.valueOf(B_N_R_15.get(i)));
					    	break;
					    case "B_N_T_15":
					        temp.add(String.valueOf(B_N_T_15.get(i)));
					    	break;
					    	
					    case "B_N_C_16":
					        temp.add(String.valueOf(B_N_C_16.get(i)));
					    	break;
					    case "B_N_S_16":
					        temp.add(String.valueOf(B_N_S_16.get(i)));
					    	break;
					    case "B_N_R_16":
					        temp.add(String.valueOf(B_N_R_16.get(i)));
					    	break;
					    case "B_N_T_16":
					        temp.add(String.valueOf(B_N_T_16.get(i)));
					    	break;
					    	
					    case "B_N_C_17":
					        temp.add(String.valueOf(B_N_C_17.get(i)));
					    	break;
					    case "B_N_S_17":
					        temp.add(String.valueOf(B_N_S_17.get(i)));
					    	break;
					    case "B_N_R_17":
					        temp.add(String.valueOf(B_N_R_17.get(i)));
					    	break;
					    case "B_N_T_17":
					        temp.add(String.valueOf(B_N_T_17.get(i)));
					    	break;
					    	
					    case "B_N_C_18":
					        temp.add(String.valueOf(B_N_C_18.get(i)));
					    	break;
					    case "B_N_S_18":
					        temp.add(String.valueOf(B_N_S_18.get(i)));
					    	break;
					    case "B_N_R_18":
					        temp.add(String.valueOf(B_N_R_18.get(i)));
					    	break;
					    case "B_N_T_18":
					        temp.add(String.valueOf(B_N_T_18.get(i)));
					    	break;  	
					   
					    case "B_N_C_19":
					        temp.add(String.valueOf(B_N_C_19.get(i)));
					    	break;
					    case "B_N_S_19":
					        temp.add(String.valueOf(B_N_S_19.get(i)));
					    	break;
					    case "B_N_R_19":
					        temp.add(String.valueOf(B_N_R_19.get(i)));
					    	break;
					    case "B_N_T_19":
					        temp.add(String.valueOf(B_N_T_19.get(i)));
					    	break; 
					    	

					    case "B_N_C_20":
					        temp.add(String.valueOf(B_N_C_20.get(i)));
					    	break;
					    case "B_N_S_20":
					        temp.add(String.valueOf(B_N_S_20.get(i)));
					    	break;
					    case "B_N_R_20":
					        temp.add(String.valueOf(B_N_R_20.get(i)));
					    	break;
					    case "B_N_T_20":
					        temp.add(String.valueOf(B_N_T_20.get(i)));
					    	break;  	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    	
					    case "M_N_C_2":
					        temp.add(String.valueOf(M_N_C_2.get(i)));
					    	break;
					    case "M_N_S_2":
					        temp.add(String.valueOf(M_N_S_2.get(i)));
					    	break;
					    case "M_N_R_2":
					        temp.add(String.valueOf(M_N_R_2.get(i)));
					    	break;
					    case "M_N_T_2":
					        temp.add(String.valueOf(M_N_T_2.get(i)));
					    	break; 
					    	
					    case "M_N_C_3":
					        temp.add(String.valueOf(M_N_C_3.get(i)));
					    	break;
					    case "M_N_S_3":
					        temp.add(String.valueOf(M_N_S_3.get(i)));
					    	break;
					    case "M_N_R_3":
					        temp.add(String.valueOf(M_N_R_3.get(i)));
					    	break;
					    case "M_N_T_3":
					        temp.add(String.valueOf(M_N_T_3.get(i)));
					    	break; 
					  
					    case "M_N_C_4":
					        temp.add(String.valueOf(M_N_C_4.get(i)));
					    	break;
					    case "M_N_S_4":
					        temp.add(String.valueOf(M_N_S_4.get(i)));
					    	break;
					    case "M_N_R_4":
					        temp.add(String.valueOf(M_N_R_4.get(i)));
					    	break;
					    case "M_N_T_4":
					        temp.add(String.valueOf(M_N_T_4.get(i)));
					    	break; 
					    	
					    case "M_N_C_5":
					        temp.add(String.valueOf(M_N_C_5.get(i)));
					    	break;
					    case "M_N_S_5":
					        temp.add(String.valueOf(M_N_S_5.get(i)));
					    	break;
					    case "M_N_R_5":
					        temp.add(String.valueOf(M_N_R_5.get(i)));
					    	break;
					    case "M_N_T_5":
					        temp.add(String.valueOf(M_N_T_5.get(i)));
					    	break; 
					
					    case "M_N_C_6":
					        temp.add(String.valueOf(M_N_C_6.get(i)));
					    	break;
					    case "M_N_S_6":
					        temp.add(String.valueOf(M_N_S_6.get(i)));
					    	break;
					    case "M_N_R_6":
					        temp.add(String.valueOf(M_N_R_6.get(i)));
					    	break;
					    case "M_N_T_6":
					        temp.add(String.valueOf(M_N_T_6.get(i)));
					    	break; 
					    	 	
					    case "M_N_C_7":
					        temp.add(String.valueOf(M_N_C_7.get(i)));
					    	break;
					    case "M_N_S_7":
					        temp.add(String.valueOf(M_N_S_7.get(i)));
					    	break;
					    case "M_N_R_7":
					        temp.add(String.valueOf(M_N_R_7.get(i)));
					    	break;
					    case "M_N_T_7":
					        temp.add(String.valueOf(M_N_T_7.get(i)));
					    	break;  	
				
					    case "M_N_C_8":
					        temp.add(String.valueOf(M_N_C_8.get(i)));
					    	break;
					    case "M_N_S_8":
					        temp.add(String.valueOf(M_N_S_8.get(i)));
					    	break;
					    case "M_N_R_8":
					        temp.add(String.valueOf(M_N_R_8.get(i)));
					    	break;
					    case "M_N_T_8":
					        temp.add(String.valueOf(M_N_T_8.get(i)));
					    	break; 
					    	
					    case "M_N_C_9":
					        temp.add(String.valueOf(M_N_C_9.get(i)));
					    	break;
					    case "M_N_S_9":
					        temp.add(String.valueOf(M_N_S_9.get(i)));
					    	break;
					    case "M_N_R_9":
					        temp.add(String.valueOf(M_N_R_9.get(i)));
					    	break;
					    case "M_N_T_9":
					        temp.add(String.valueOf(M_N_T_9.get(i)));
					    	break; 
					    	
					    case "M_N_C_10":
					        temp.add(String.valueOf(M_N_C_10.get(i)));
					    	break;
					    case "M_N_S_10":
					        temp.add(String.valueOf(M_N_S_10.get(i)));
					    	break;
					    case "M_N_R_10":
					        temp.add(String.valueOf(M_N_R_10.get(i)));
					    	break;
					    case "M_N_T_10":
					        temp.add(String.valueOf(M_N_T_10.get(i)));
					    	break; 
					    	
					  
					    case "M_N_C_11":
					        temp.add(String.valueOf(M_N_C_11.get(i)));
					    	break;
					    case "M_N_S_11":
					        temp.add(String.valueOf(M_N_S_11.get(i)));
					    	break;
					    case "M_N_R_11":
					        temp.add(String.valueOf(M_N_R_11.get(i)));
					    	break;
					    case "M_N_T_11":
					        temp.add(String.valueOf(M_N_T_2.get(i)));
					    	break; 
					    	
					    case "M_N_C_12":
					        temp.add(String.valueOf(M_N_C_12.get(i)));
					    	break;
					    case "M_N_S_12":
					        temp.add(String.valueOf(M_N_S_12.get(i)));
					    	break;
					    case "M_N_R_12":
					        temp.add(String.valueOf(M_N_R_12.get(i)));
					    	break;
					    case "M_N_T_12":
					        temp.add(String.valueOf(M_N_T_12.get(i)));
					    	break; 
					  
					    case "M_N_C_13":
					        temp.add(String.valueOf(M_N_C_13.get(i)));
					    	break;
					    case "M_N_S_13":
					        temp.add(String.valueOf(M_N_S_13.get(i)));
					    	break;
					    case "M_N_R_13":
					        temp.add(String.valueOf(M_N_R_13.get(i)));
					    	break;
					    case "M_N_T_13":
					        temp.add(String.valueOf(M_N_T_13.get(i)));
					    	break; 
					    	
					    case "M_N_C_14":
					        temp.add(String.valueOf(M_N_C_14.get(i)));
					    	break;
					    case "M_N_S_14":
					        temp.add(String.valueOf(M_N_S_14.get(i)));
					    	break;
					    case "M_N_R_14":
					        temp.add(String.valueOf(M_N_R_14.get(i)));
					    	break;
					    case "M_N_T_14":
					        temp.add(String.valueOf(M_N_T_14.get(i)));
					    	break; 
					
					    case "M_N_C_15":
					        temp.add(String.valueOf(M_N_C_15.get(i)));
					    	break;
					    case "M_N_S_15":
					        temp.add(String.valueOf(M_N_S_15.get(i)));
					    	break;
					    case "M_N_R_15":
					        temp.add(String.valueOf(M_N_R_15.get(i)));
					    	break;
					    case "M_N_T_15":
					        temp.add(String.valueOf(M_N_T_15.get(i)));
					    	break; 
					    	 	
					    case "M_N_C_16":
					        temp.add(String.valueOf(M_N_C_16.get(i)));
					    	break;
					    case "M_N_S_16":
					        temp.add(String.valueOf(M_N_S_16.get(i)));
					    	break;
					    case "M_N_R_16":
					        temp.add(String.valueOf(M_N_R_16.get(i)));
					    	break;
					    case "M_N_T_16":
					        temp.add(String.valueOf(M_N_T_16.get(i)));
					    	break;  	
				
					    case "M_N_C_17":
					        temp.add(String.valueOf(M_N_C_17.get(i)));
					    	break;
					    case "M_N_S_17":
					        temp.add(String.valueOf(M_N_S_17.get(i)));
					    	break;
					    case "M_N_R_17":
					        temp.add(String.valueOf(M_N_R_17.get(i)));
					    	break;
					    case "M_N_T_17":
					        temp.add(String.valueOf(M_N_T_17.get(i)));
					    	break; 
					    	
					    case "M_N_C_18":
					        temp.add(String.valueOf(M_N_C_18.get(i)));
					    	break;
					    case "M_N_S_18":
					        temp.add(String.valueOf(M_N_S_18.get(i)));
					    	break;
					    case "M_N_R_18":
					        temp.add(String.valueOf(M_N_R_18.get(i)));
					    	break;
					    case "M_N_T_18":
					        temp.add(String.valueOf(M_N_T_18.get(i)));
					    	break; 
					    	
					    case "M_N_C_19":
					        temp.add(String.valueOf(M_N_C_19.get(i)));
					    	break;
					    case "M_N_S_19":
					        temp.add(String.valueOf(M_N_S_19.get(i)));
					    	break;
					    case "M_N_R_19":
					        temp.add(String.valueOf(M_N_R_19.get(i)));
					    	break;
					    case "M_N_T_19":
					        temp.add(String.valueOf(M_N_T_19.get(i)));
					    	break;  	
					    	
					    case "M_N_C_20":
					        temp.add(String.valueOf(M_N_C_20.get(i)));
					    	break;
					    case "M_N_S_20":
					        temp.add(String.valueOf(M_N_S_20.get(i)));
					    	break;
					    case "M_N_R_20":
					        temp.add(String.valueOf(M_N_R_20.get(i)));
					    	break;
					    case "M_N_T_20":
					        temp.add(String.valueOf(M_N_T_20.get(i)));
					    	break;  	
					    	 	
					    	
					    	
					    	
					    case "D_N_C_2":
					        temp.add(String.valueOf(D_N_C_2.get(i)));
					    	break;
					    case "D_N_S_2":
					        temp.add(String.valueOf(D_N_S_2.get(i)));
					    	break;
					    case "D_N_R_2":
					        temp.add(String.valueOf(D_N_R_2.get(i)));
					    	break;
					    case "D_N_T_2":
					        temp.add(String.valueOf(D_N_T_2.get(i)));
					    	break; 
					    	
					    case "D_N_C_3":
					        temp.add(String.valueOf(D_N_C_3.get(i)));
					    	break;
					    case "D_N_S_3":
					        temp.add(String.valueOf(D_N_S_3.get(i)));
					    	break;
					    case "D_N_R_3":
					        temp.add(String.valueOf(D_N_R_3.get(i)));
					    	break;
					    case "D_N_T_3":
					        temp.add(String.valueOf(D_N_T_3.get(i)));
					    	break; 
					  
					    case "D_N_C_4":
					        temp.add(String.valueOf(D_N_C_4.get(i)));
					    	break;
					    case "D_N_S_4":
					        temp.add(String.valueOf(D_N_S_4.get(i)));
					    	break;
					    case "D_N_R_4":
					        temp.add(String.valueOf(D_N_R_4.get(i)));
					    	break;
					    case "D_N_T_4":
					        temp.add(String.valueOf(D_N_T_4.get(i)));
					    	break; 
					    	
					    case "D_N_C_5":
					        temp.add(String.valueOf(D_N_C_5.get(i)));
					    	break;
					    case "D_N_S_5":
					        temp.add(String.valueOf(D_N_S_5.get(i)));
					    	break;
					    case "D_N_R_5":
					        temp.add(String.valueOf(D_N_R_5.get(i)));
					    	break;
					    case "D_N_T_5":
					        temp.add(String.valueOf(D_N_T_5.get(i)));
					    	break; 
					
					    case "D_N_C_6":
					        temp.add(String.valueOf(D_N_C_6.get(i)));
					    	break;
					    case "D_N_S_6":
					        temp.add(String.valueOf(D_N_S_6.get(i)));
					    	break;
					    case "D_N_R_6":
					        temp.add(String.valueOf(D_N_R_6.get(i)));
					    	break;
					    case "D_N_T_6":
					        temp.add(String.valueOf(D_N_T_6.get(i)));
					    	break; 
					    	 	
					    case "D_N_C_7":
					        temp.add(String.valueOf(D_N_C_7.get(i)));
					    	break;
					    case "D_N_S_7":
					        temp.add(String.valueOf(D_N_S_7.get(i)));
					    	break;
					    case "D_N_R_7":
					        temp.add(String.valueOf(D_N_R_7.get(i)));
					    	break;
					    case "D_N_T_7":
					        temp.add(String.valueOf(D_N_T_7.get(i)));
					    	break;  	
				
					    case "D_N_C_8":
					        temp.add(String.valueOf(D_N_C_8.get(i)));
					    	break;
					    case "D_N_S_8":
					        temp.add(String.valueOf(D_N_S_8.get(i)));
					    	break;
					    case "D_N_R_8":
					        temp.add(String.valueOf(D_N_R_8.get(i)));
					    	break;
					    case "D_N_T_8":
					        temp.add(String.valueOf(D_N_T_8.get(i)));
					    	break; 
					    	
					    case "D_N_C_9":
					        temp.add(String.valueOf(D_N_C_9.get(i)));
					    	break;
					    case "D_N_S_9":
					        temp.add(String.valueOf(D_N_S_9.get(i)));
					    	break;
					    case "D_N_R_9":
					        temp.add(String.valueOf(D_N_R_9.get(i)));
					    	break;
					    case "D_N_T_9":
					        temp.add(String.valueOf(D_N_T_9.get(i)));
					    	break; 
					    	
					    case "D_N_C_10":
					        temp.add(String.valueOf(D_N_C_10.get(i)));
					    	break;
					    case "D_N_S_10":
					        temp.add(String.valueOf(D_N_S_10.get(i)));
					    	break;
					    case "D_N_R_10":
					        temp.add(String.valueOf(D_N_R_10.get(i)));
					    	break;
					    case "D_N_T_10":
					        temp.add(String.valueOf(D_N_T_10.get(i)));
					    	break; 
					   
					    	
					    case "D_N_C_11":
					        temp.add(String.valueOf(D_N_C_11.get(i)));
					    	break;
					    case "D_N_S_11":
					        temp.add(String.valueOf(D_N_S_11.get(i)));
					    	break;
					    case "D_N_R_11":
					        temp.add(String.valueOf(D_N_R_11.get(i)));
					    	break;
					    case "D_N_T_11":
					        temp.add(String.valueOf(D_N_T_11.get(i)));
					    	break; 
					    	
					    case "D_N_C_12":
					        temp.add(String.valueOf(D_N_C_12.get(i)));
					    	break;
					    case "D_N_S_12":
					        temp.add(String.valueOf(D_N_S_12.get(i)));
					    	break;
					    case "D_N_R_12":
					        temp.add(String.valueOf(D_N_R_12.get(i)));
					    	break;
					    case "D_N_T_12":
					        temp.add(String.valueOf(D_N_T_12.get(i)));
					    	break; 
					  
					    case "D_N_C_13":
					        temp.add(String.valueOf(D_N_C_13.get(i)));
					    	break;
					    case "D_N_S_13":
					        temp.add(String.valueOf(D_N_S_13.get(i)));
					    	break;
					    case "D_N_R_13":
					        temp.add(String.valueOf(D_N_R_13.get(i)));
					    	break;
					    case "D_N_T_13":
					        temp.add(String.valueOf(D_N_T_13.get(i)));
					    	break; 
					    	
					    case "D_N_C_14":
					        temp.add(String.valueOf(D_N_C_14.get(i)));
					    	break;
					    case "D_N_S_14":
					        temp.add(String.valueOf(D_N_S_14.get(i)));
					    	break;
					    case "D_N_R_14":
					        temp.add(String.valueOf(D_N_R_14.get(i)));
					    	break;
					    case "D_N_T_14":
					        temp.add(String.valueOf(D_N_T_14.get(i)));
					    	break; 
					
					    case "D_N_C_15":
					        temp.add(String.valueOf(D_N_C_15.get(i)));
					    	break;
					    case "D_N_S_15":
					        temp.add(String.valueOf(D_N_S_15.get(i)));
					    	break;
					    case "D_N_R_15":
					        temp.add(String.valueOf(D_N_R_15.get(i)));
					    	break;
					    case "D_N_T_15":
					        temp.add(String.valueOf(D_N_T_15.get(i)));
					    	break; 
					    	 	
					    case "D_N_C_16":
					        temp.add(String.valueOf(D_N_C_16.get(i)));
					    	break;
					    case "D_N_S_16":
					        temp.add(String.valueOf(D_N_S_16.get(i)));
					    	break;
					    case "D_N_R_16":
					        temp.add(String.valueOf(D_N_R_16.get(i)));
					    	break;
					    case "D_N_T_16":
					        temp.add(String.valueOf(D_N_T_16.get(i)));
					    	break;  	
				
					    case "D_N_C_17":
					        temp.add(String.valueOf(D_N_C_17.get(i)));
					    	break;
					    case "D_N_S_17":
					        temp.add(String.valueOf(D_N_S_17.get(i)));
					    	break;
					    case "D_N_R_17":
					        temp.add(String.valueOf(D_N_R_17.get(i)));
					    	break;
					    case "D_N_T_17":
					        temp.add(String.valueOf(D_N_T_17.get(i)));
					    	break; 
					    	
					    case "D_N_C_18":
					        temp.add(String.valueOf(D_N_C_18.get(i)));
					    	break;
					    case "D_N_S_18":
					        temp.add(String.valueOf(D_N_S_18.get(i)));
					    	break;
					    case "D_N_R_18":
					        temp.add(String.valueOf(D_N_R_18.get(i)));
					    	break;
					    case "D_N_T_18":
					        temp.add(String.valueOf(D_N_T_18.get(i)));
					    	break; 
					    	
					    case "D_N_C_19":
					        temp.add(String.valueOf(D_N_C_19.get(i)));
					    	break;
					    case "D_N_S_19":
					        temp.add(String.valueOf(D_N_S_19.get(i)));
					    	break;
					    case "D_N_R_19":
					        temp.add(String.valueOf(D_N_R_19.get(i)));
					    	break;
					    case "D_N_T_19":
					        temp.add(String.valueOf(D_N_T_19.get(i)));
					    	break;  	
					  
					    case "D_N_C_20":
					        temp.add(String.valueOf(D_N_C_20.get(i)));
					    	break;
					    case "D_N_S_20":
					        temp.add(String.valueOf(D_N_S_20.get(i)));
					    	break;
					    case "D_N_R_20":
					        temp.add(String.valueOf(D_N_R_20.get(i)));
					    	break;
					    case "D_N_T_20":
					        temp.add(String.valueOf(D_N_T_20.get(i)));
					    	break;  	
					    	
					    //MA Relative
					    case "M_R_C_2":
					        temp.add(String.valueOf(M_R_C_2.get(i)));
					    	break;
					    case "M_R_S_2":
					        temp.add(String.valueOf(M_R_S_2.get(i)));
					    	break;
					    case "M_R_R_2":
					        temp.add(String.valueOf(M_R_R_2.get(i)));
					    	break;
					    case "M_R_T_2":
					        temp.add(String.valueOf(M_R_T_2.get(i)));
					    	break; 
					    	
					    case "M_R_C_3":
					        temp.add(String.valueOf(M_R_C_3.get(i)));
					    	break;
					    case "M_R_S_3":
					        temp.add(String.valueOf(M_R_S_3.get(i)));
					    	break;
					    case "M_R_R_3":
					        temp.add(String.valueOf(M_R_R_3.get(i)));
					    	break;
					    case "M_R_T_3":
					        temp.add(String.valueOf(M_R_T_3.get(i)));
					    	break; 
					  
					    case "M_R_C_4":
					        temp.add(String.valueOf(M_R_C_4.get(i)));
					    	break;
					    case "M_R_S_4":
					        temp.add(String.valueOf(M_R_S_4.get(i)));
					    	break;
					    case "M_R_R_4":
					        temp.add(String.valueOf(M_R_R_4.get(i)));
					    	break;
					    case "M_R_T_4":
					        temp.add(String.valueOf(M_R_T_4.get(i)));
					    	break; 
					    	
					    case "M_R_C_5":
					        temp.add(String.valueOf(M_R_C_5.get(i)));
					    	break;
					    case "M_R_S_5":
					        temp.add(String.valueOf(M_R_S_5.get(i)));
					    	break;
					    case "M_R_R_5":
					        temp.add(String.valueOf(M_R_R_5.get(i)));
					    	break;
					    case "M_R_T_5":
					        temp.add(String.valueOf(M_R_T_5.get(i)));
					    	break; 
					
					    case "M_R_C_6":
					        temp.add(String.valueOf(M_R_C_6.get(i)));
					    	break;
					    case "M_R_S_6":
					        temp.add(String.valueOf(M_R_S_6.get(i)));
					    	break;
					    case "M_R_R_6":
					        temp.add(String.valueOf(M_R_R_6.get(i)));
					    	break;
					    case "M_R_T_6":
					        temp.add(String.valueOf(M_R_T_6.get(i)));
					    	break; 
					    	 	
					    case "M_R_C_7":
					        temp.add(String.valueOf(M_R_C_7.get(i)));
					    	break;
					    case "M_R_S_7":
					        temp.add(String.valueOf(M_R_S_7.get(i)));
					    	break;
					    case "M_R_R_7":
					        temp.add(String.valueOf(M_R_R_7.get(i)));
					    	break;
					    case "M_R_T_7":
					        temp.add(String.valueOf(M_R_T_7.get(i)));
					    	break; 
				
					    case "M_R_C_8":
					        temp.add(String.valueOf(M_R_C_8.get(i)));
					    	break;
					    case "M_R_S_8":
					        temp.add(String.valueOf(M_R_S_8.get(i)));
					    	break;
					    case "M_R_R_8":
					        temp.add(String.valueOf(M_R_R_8.get(i)));
					    	break;
					    case "M_R_T_8":
					        temp.add(String.valueOf(M_R_T_8.get(i)));
					    	break; 
					    	
					    case "M_R_C_9":
					        temp.add(String.valueOf(M_R_C_9.get(i)));
					    	break;
					    case "M_R_S_9":
					        temp.add(String.valueOf(M_R_S_9.get(i)));
					    	break;
					    case "M_R_R_9":
					        temp.add(String.valueOf(M_R_R_9.get(i)));
					    	break;
					    case "M_R_T_9":
					        temp.add(String.valueOf(M_R_T_9.get(i)));
					    	break; 
					    	
					    case "M_R_C_10":
					        temp.add(String.valueOf(M_R_C_10.get(i)));
					    	break;
					    case "M_R_S_10":
					        temp.add(String.valueOf(M_R_S_10.get(i)));
					    	break;
					    case "M_R_R_10":
					        temp.add(String.valueOf(M_R_R_10.get(i)));
					    	break;
					    case "M_R_T_10":
					        temp.add(String.valueOf(M_R_T_10.get(i)));
					    	break; 
					    	
					    case "M_R_C_11":
					        temp.add(String.valueOf(M_R_C_11.get(i)));
					    	break;
					    case "M_R_S_11":
					        temp.add(String.valueOf(M_R_S_11.get(i)));
					    	break;
					    case "M_R_R_11":
					        temp.add(String.valueOf(M_R_R_11.get(i)));
					    	break;
					    case "M_R_T_11":
					        temp.add(String.valueOf(M_R_T_11.get(i)));
					    	break; 
					    	
					    case "M_R_C_12":
					        temp.add(String.valueOf(M_R_C_12.get(i)));
					    	break;
					    case "M_R_S_12":
					        temp.add(String.valueOf(M_R_S_12.get(i)));
					    	break;
					    case "M_R_R_12":
					        temp.add(String.valueOf(M_R_R_12.get(i)));
					    	break;
					    case "M_R_T_12":
					        temp.add(String.valueOf(M_R_T_12.get(i)));
					    	break; 
					  
					    case "M_R_C_13":
					        temp.add(String.valueOf(M_R_C_4.get(i)));
					    	break;
					    case "M_R_S_13":
					        temp.add(String.valueOf(M_R_S_13.get(i)));
					    	break;
					    case "M_R_R_13":
					        temp.add(String.valueOf(M_R_R_13.get(i)));
					    	break;
					    case "M_R_T_13":
					        temp.add(String.valueOf(M_R_T_13.get(i)));
					    	break; 
					    	
					    case "M_R_C_14":
					        temp.add(String.valueOf(M_R_C_14.get(i)));
					    	break;
					    case "M_R_S_14":
					        temp.add(String.valueOf(M_R_S_14.get(i)));
					    	break;
					    case "M_R_R_14":
					        temp.add(String.valueOf(M_R_R_14.get(i)));
					    	break;
					    case "M_R_T_14":
					        temp.add(String.valueOf(M_R_T_14.get(i)));
					    	break; 
					
					    case "M_R_C_15":
					        temp.add(String.valueOf(M_R_C_15.get(i)));
					    	break;
					    case "M_R_S_15":
					        temp.add(String.valueOf(M_R_S_15.get(i)));
					    	break;
					    case "M_R_R_15":
					        temp.add(String.valueOf(M_R_R_15.get(i)));
					    	break;
					    case "M_R_T_15":
					        temp.add(String.valueOf(M_R_T_15.get(i)));
					    	break; 
					    	 	
					    case "M_R_C_16":
					        temp.add(String.valueOf(M_R_C_16.get(i)));
					    	break;
					    case "M_R_S_16":
					        temp.add(String.valueOf(M_R_S_16.get(i)));
					    	break;
					    case "M_R_R_16":
					        temp.add(String.valueOf(M_R_R_16.get(i)));
					    	break;
					    case "M_R_T_16":
					        temp.add(String.valueOf(M_R_T_16.get(i)));
					    	break; 
				
					    case "M_R_C_17":
					        temp.add(String.valueOf(M_R_C_17.get(i)));
					    	break;
					    case "M_R_S_17":
					        temp.add(String.valueOf(M_R_S_17.get(i)));
					    	break;
					    case "M_R_R_17":
					        temp.add(String.valueOf(M_R_R_17.get(i)));
					    	break;
					    case "M_R_T_17":
					        temp.add(String.valueOf(M_R_T_17.get(i)));
					    	break; 
					    	
					    case "M_R_C_18":
					        temp.add(String.valueOf(M_R_C_18.get(i)));
					    	break;
					    case "M_R_S_18":
					        temp.add(String.valueOf(M_R_S_18.get(i)));
					    	break;
					    case "M_R_R_18":
					        temp.add(String.valueOf(M_R_R_18.get(i)));
					    	break;
					    case "M_R_T_18":
					        temp.add(String.valueOf(M_R_T_18.get(i)));
					    	break; 
					    	
					    case "M_R_C_19":
					        temp.add(String.valueOf(M_R_C_19.get(i)));
					    	break;
					    case "M_R_S_19":
					        temp.add(String.valueOf(M_R_S_19.get(i)));
					    	break;
					    case "M_R_R_19":
					        temp.add(String.valueOf(M_R_R_19.get(i)));
					    	break;
					    case "M_R_T_19":
					        temp.add(String.valueOf(M_R_T_19.get(i)));
					    	break; 
					    	   	
					     	
					    case "M_R_C_20":
					        temp.add(String.valueOf(M_R_C_20.get(i)));
					    	break;
					    case "M_R_S_20":
					        temp.add(String.valueOf(M_R_S_20.get(i)));
					    	break;
					    case "M_R_R_20":
					        temp.add(String.valueOf(M_R_R_20.get(i)));
					    	break;
					    case "M_R_T_20":
					        temp.add(String.valueOf(M_R_T_20.get(i)));
					    	break;   	
					    	
 					}
 					
				}
				
                  
				temp.add(feature_target.get(i));		
			}	
			//temp.add(records.get(i).get(records.get(i).size()-1));	
			result.add(temp);
		}		
		try {
		writeCSV("", output_filename,result);
		} catch (IOException e) {
			System.out.println("[ERROR] I/O Exception.");
			e.printStackTrace();
		}
	}
	
    public static HashMap<Integer, String> featureExtraction_target(ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, String> result = new HashMap<>();
    	
    	int index_of_target_att = records.get(0).size()-1;
    	for (int i = 1; i < records.size(); i++) {
    	    if (i==1) {
    	    	result.put(i, "Rise"); 
    	    	continue;
    	    }
//    	    System.out.println(i);
//    	    System.out.println(Double.parseDouble(records.get(i).get(index_of_target_att)));
    	   
    	    if (Double.parseDouble(records.get(i).get(index_of_target_att))- Double.parseDouble(records.get(i-1).get(index_of_target_att)) >= 0 ) {
    	    	result.put(i, "Rise");     
    	    } else {
    	    	result.put(i, "Down");  
    	    }	
    	}    	  
    	return result;  
    	
    }
    
    public static HashMap<Integer, String> featureExtraction_target_user_defined(ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, String> result = new HashMap<>();
    	
    	int index_of_target_att = records.get(0).size()-1;
    	for (int i = 1; i < records.size(); i++) {
    	    if (i==1) {  
    	    	result.put(i, "None_Rise");
    	    	continue;
    	    }
    	    
    	    double price_now = Double.parseDouble(records.get(i).get(index_of_target_att));
    	    double price_pre = Double.parseDouble(records.get(i-1).get(index_of_target_att));
    	    if (0 < (price_now-price_pre) && (price_now-price_pre) < 30) {
    	    	result.put(i, "None_Rise");     
    	    } else if (30 <= (price_now-price_pre) && (price_now-price_pre) <= 100) {
    	    	result.put(i, "Rise");  
    	    } else if ((price_now-price_pre) > 100) {
    	    	result.put(i, "Most_Rise");  
    	    } else if (-30 < (price_now-price_pre) && (price_now-price_pre) <=0) {
    	    	result.put(i, "None_Down");  
    	    } else if (-100 <= (price_now-price_pre) && (price_now-price_pre) <= -30) {
    	    	result.put(i, "Down");
    	    } else if ((price_now-price_pre) < -100) {
    	    	result.put(i, "Most_Down");
    	    }
    	}    	  
    	return result;  
    	
    }
           
    public static HashMap<Integer, String> MACD(int tl, int sl, int ll, String att, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, String> result = new HashMap<>(); 
    	for (int i = 1; i < records.size(); i++) {
    	    double MACD = DIF(i, sl, ll, records) - DEM(i, sl, ll, tl, records);        	
    		if (MACD < 0) {
    			result.put(i, "MACD_" + att.charAt(0) + sl + ll+"_0");
    		} else {
    			result.put(i, "MACD_" + att.charAt(0) + sl + ll+"_1");			
    		}
    	}
    	return result;
    } 
    
    public static HashMap<Integer, Double> MACD_weka(int tl, int sl, int ll, String att, ArrayList<ArrayList<String>> records) {
    	HashMap<Integer, Double> result = new HashMap<>(); 
    	for (int i = 1; i < records.size(); i++) {
    	    double MACD = DIF(i, sl, ll, records) - DEM(i, sl, ll, tl, records);        	
    		result.put(i, MACD);
    	}
    	return result;
    } 
     
    public static double EMA(int t, int l, ArrayList<ArrayList<String>> records, String s) {
    	if (t == 0) {  
    		return 0.0;
    	}
    	int col = 2;
    	double alpha = 2/(double)(l+1);
    	double p = Double.parseDouble(records.get(t).get(col));
        if (s.equals("sl")) {
        	temp_sl.put(0, 0.0);
            temp_sl.put(t, temp_sl.get(t-1) + alpha*(p - temp_sl.get(t-1)));
            return temp_sl.get(t);
        } else {  	
        	temp_ll.put(0, 0.0); 
            temp_ll.put(t, temp_ll.get(t-1) + alpha*(p - temp_ll.get(t-1)));
            return temp_ll.get(t);
        }  
    	
    }
    
    public static double DIF(int t, int sl, int ll, ArrayList<ArrayList<String>> records) {
        return EMA(t, sl, records, "sl") - EMA(t, ll, records, "ll"); 	
    }
    
    public static double DEM(int t, int sl, int ll, int tl, ArrayList<ArrayList<String>> records) {
        return 	(DIF(t, sl, ll, records) + DIF(t-1, sl, ll, records))/(double) tl;
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
}
