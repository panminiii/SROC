 import java.io.*;

public class ExtractAveragePrecision {	

	public static void main(String[] args) throws Exception {
		String path="MAPCC/";
		FileWriter writer = new FileWriter("MAPFiles/MAP"+"name");
		for (File file : new File(path).listFiles()) {
			extractAveragePrecision(file,writer);
			writer.write("\r\n");
			System.out.println(file);
		}
		writer.close();
	}

	public static void extractAveragePrecision(File file, FileWriter writer) throws Exception {

		if (file.isDirectory() || !file.getName().endsWith(".txt")) {
			return;
		}

		LineNumberReader reader = new LineNumberReader(new FileReader(file));

//		FileWriter writer = new FileWriter("MAPFiles/MAP-" + file.getName());

		String line = null;
		writer.write(file.getName()+ "\t");
		while ( (line = reader.readLine()) != null ) {

			if (line.indexOf("Average Precision:") != -1) {
				writer.write(line.substring(line.indexOf(":")+1, line.length()).trim() + "\t");
			}
		}


		reader.close();
	}
}