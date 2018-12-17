package workshop.server.common;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Logger implements ILogger { 
    private static final String TIME_FORMAT = "dd.MM.yyyy HH:mm:ss";
    
    private final String mFileName;
    private final BufferedWriter mWriter;
    private static final ILogger instance_ = new Logger("log.txt");
    
    public static ILogger getInstance() {
        return instance_;
    }
    
    Logger(final String fileName) {
        mFileName = fileName;

        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(fileName, true));
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        mWriter = writer;
    }
    
    private String getCurrentTime() {
        Calendar calendar = Calendar.getInstance();
        SimpleDateFormat dateFormat = new SimpleDateFormat(TIME_FORMAT);
        return dateFormat.format(calendar.getTime());
    }

    @Override
    public synchronized boolean log(String message) {
        if(mWriter != null) {
            try {
                mWriter.write(getCurrentTime() + " >> " + message + "\n");
                mWriter.flush();
                return true;
            } catch (IOException e) {
                System.out.println(e.toString());
                return false;
            }
        } else {
            return false;
        }
    }

    @Override
    public String getLog() {
        BufferedReader br;
        StringBuilder sb = new StringBuilder();
        try {
            br = new BufferedReader(new FileReader(mFileName));
            String line;
            while ((line = br.readLine()) != null) {
                sb.append(line).append("\n");
            }
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        return sb.toString();
    }

    @Override
    public boolean clear() {   
        try {
            PrintWriter pw = new PrintWriter(mFileName);
            pw.close();
            return true;
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return false;
        }
    }   
}
