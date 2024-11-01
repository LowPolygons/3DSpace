import java.util.*;
import java.io.File;  // Import the File class
import java.io.FileWriter;  // Import the FileWriter class
import java.io.IOException;  // Import the IOException class to handle errors
import java.io.FileNotFoundException;

public class reading {

    public static void main (String[] args) {
        File file = new File("output1.txt");
        String[] totalTime = new String[10000];
        String[] comp = new String[10000];
        String[] comms = new String[10000];

        int totalTimeCounter = 0;
        int totalCompCounter = 0;
        int totalCommsCounter = 0;

        try {
            Scanner scanner = new Scanner(file);
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();

                if (line.contains("Main Program Time:")) {
                    totalTime[totalTimeCounter] = line;
                    totalTimeCounter++;
                }                
                if (line.contains("Total Computation Time:")) {
                    comp[totalCompCounter] = line;
                    totalCompCounter++;
                }                
                if (line.contains("Total Communication Time:")) {
                    comms[totalCommsCounter] = line;
                    totalCommsCounter++;
                }

            }          
        } catch (FileNotFoundException error) {
            System.out.println("File not found.");
        }

        // file = new File("output2.txt");
        // try {
        //     Scanner scanner = new Scanner(file);
        //     while (scanner.hasNextLine()) {
        //         String line = scanner.nextLine();

        //         if (line.contains("Main Program Time:")) {
        //             totalTime[totalTimeCounter] = line;
        //             totalTimeCounter++;
        //         }                
        //         if (line.contains("Total Computation Time:")) {
        //             comp[totalCompCounter] = line;
        //             totalCompCounter++;
        //         }                
        //         if (line.contains("Total Communication Time:")) {
        //             comms[totalCommsCounter] = line;
        //             totalCommsCounter++;
        //         }

        //     }          
        // } catch (FileNotFoundException error) {
        //     System.out.println("File not found.");
        // }
        
        // file = new File("output3.txt");
        // try {
        //     Scanner scanner = new Scanner(file);
        //     while (scanner.hasNextLine()) {
        //         String line = scanner.nextLine();

        //         if (line.contains("Main Program Time:")) {
        //             totalTime[totalTimeCounter] = line;
        //             totalTimeCounter++;
        //         }                
        //         if (line.contains("Total Computation Time:")) {
        //             comp[totalCompCounter] = line;
        //             totalCompCounter++;
        //         }                
        //         if (line.contains("Total Communication Time:")) {
        //             comms[totalCommsCounter] = line;
        //             totalCommsCounter++;
        //         }

        //     }          
        // } catch (FileNotFoundException error) {
        //     System.out.println("File not found.");
        // }
        
        double[] dtotal = new double[totalTimeCounter];
        double[] dcomms = new double[totalTimeCounter];
        double[] dcomps = new double[totalTimeCounter];
    
        for (int i = 0; i < totalTimeCounter; i++) {
            String currTotal = totalTime[i];
            String currComs = comms[i];
            String currComp = comp[i];
            
            //System.out.println(currTotal + "," + currComs + ", " + currComp);
            double _total = Double.parseDouble(currTotal.substring(currTotal.indexOf("Time:")+5, currTotal.length()));
            double _comms = Double.parseDouble(currComs.substring(currComs.indexOf("Time:")+5, currComs.length()));
            double _comps = Double.parseDouble(currComp.substring(currComp.indexOf("Time:")+5, currComp.length()));

            dtotal[i] = _total;
            dcomms[i] = _comms;
            dcomps[i] = _comps;
        }

        double[] atotal = new double[totalTimeCounter];
        double[] acomms = new double[totalTimeCounter];
        double[] acomps = new double[totalTimeCounter];


        int numMods = 0;
        double currTotal = 0;
        double currComs = 0;
        double currComp = 0;
        
        for (int i = 0; i < totalTimeCounter; i++) {
            currTotal = currTotal + dtotal[i];
            currComs = currComs + dcomms[i];
            currComp = currComp + dcomps[i];

            if ((i+1) % 10 == 0) {
                atotal[numMods] = currTotal/10;
                acomms[numMods] = currComs/10;
                acomps[numMods] = currComp/10;

                numMods = numMods + 1;

                currTotal = 0;
                currComs = 0;
                currComp = 0;
            }
        }


        System.out.println("Average Total Times: "+totalTimeCounter/10);
        for (int i = 0; i < totalTimeCounter/10; i++) {
           // System.out.println((i+1)+": "+atotal[i]+", ");
            System.out.print(atotal[i]+", ");
        }
        
        System.out.println("\n\nAverage Total Comms: "+totalCommsCounter/10);
        for (int i = 0; i < totalCommsCounter/10; i++) {
            System.out.print(acomms[i]+", ");
        }       
        System.out.println("\n\nAverage Total Comp: "+totalCompCounter/10);
        for (int i = 0; i < totalCompCounter/10; i++) {
            System.out.print(acomps[i]+", ");
        }     
    }
}