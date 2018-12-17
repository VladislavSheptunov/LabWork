package workshop.server.common;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.logging.Level;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.Test;

public class LoggerTest {    
    private ILogger logger;
    
    public LoggerTest() {
    }

    @Before
    public void setUp() {
        PrintWriter pw;
        try {
            pw = new PrintWriter("log.txt");
            pw.close();
        } catch (FileNotFoundException ex) {
            java.util.logging.Logger.getLogger(LoggerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        logger = new Logger("log.txt");
    }

    @Test
    public void testMethods() {
        testLog();

        testGetLog();
    }

    public void testLog() {
        logger.log("qwerty");
        
        String log = logger.getLog();
        
        assertTrue(log.contains("qwerty"));
    } 
    
    public void testGetLog() {
        String log = logger.getLog();
        
        assertTrue(log.contains("qwerty"));
    }
}
