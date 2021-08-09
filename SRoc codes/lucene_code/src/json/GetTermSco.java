package json;

import org.apache.commons.io.FileUtils; 
//import org.json.JSONException; 
//import org.json.JSONObject; 
import java.io.File; 
import java.io.IOException;
//import java.util.ArrayList;
import java.util.HashMap;

import org.json.JSONArray;

public class GetTermSco {
	public static void main(String args[]) throws IOException { 
		File file=new File("json/NewTermScoList.json"); 
		String content= FileUtils.readFileToString(file,"UTF-8"); 
		//JSONObject jsonObject=new JSONObject(content); 
		JSONArray jsonArray=new JSONArray(content);
		//System.out.println("list是："+jsonObject.getJSONArray("list"));
		//JSONArray ArrayList=jsonObject.getJSONArray("list");
		//ArrayList<HashMap<String,Double>> DocN_list = new ArrayList<HashMap<String,Double>>(); 
		for (int j = 0; j < 50; j++) {
		//JSONObject info = jsonArray.getJSONObject(j);
		System.out.println((j+1)+"次查询的扩展词表是："+jsonArray.getJSONArray(j)); 
		//System.out.println("年龄："+info.getDouble("age")); 
		//System.out.println("学到的技能："+info.getJSONArray("major")); 
		//System.out.println("国家："+info.getJSONObject("Nativeplace").getString("country"));
		JSONArray termsco=jsonArray.getJSONArray(j);
		HashMap<String,Double> pre_DocN = new HashMap<String,Double>();
		for(int t=0;t<30;t++){
			
			JSONArray array=termsco.getJSONArray(t);
			pre_DocN.put(array.getString(0),array.getDouble(1));
			//DocN_list.add(Vecter_DocN);
		}
		System.out.println(pre_DocN);
		}
		 

		
		
		} 

}
