package json;

import org.apache.commons.io.FileUtils; 
//import org.json.JSONException; 
import org.json.JSONObject; 
import java.io.File; 
import java.io.IOException; 
import org.json.JSONArray;

public class Demo {
	public static void main(String args[]) throws IOException { 
		File file=new File("json/demo.json"); 
		String content= FileUtils.readFileToString(file,"UTF-8"); 
		//JSONObject jsonObject=new JSONObject(content); 
		JSONArray jsonArray=new JSONArray(content);
		//System.out.println("list是："+jsonObject.getJSONArray("list"));
		//JSONArray ArrayList=jsonObject.getJSONArray("list");
		JSONObject info = jsonArray.getJSONObject(0);
		System.out.println("姓名是："+info.getString("name")); 
		System.out.println("年龄："+info.getDouble("age")); 
		System.out.println("学到的技能："+info.getJSONArray("major")); 
		System.out.println("国家："+info.getJSONObject("Nativeplace").getString("country"));
		
		
		} 
	}
