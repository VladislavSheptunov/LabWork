package workshop.client.head;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import org.junit.AfterClass;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;
import workshop.client.user.UController;
import workshop.common.ISocketConnection;
import workshop.common.SocketConnection;

public class HControllerTest {
    private static final String IP = "localhost";
    private static final int PORT = 8888;
    private static final ISocketConnection socketConnection = new SocketConnection();
    private static HController hController = null;

    
    public HControllerTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        if (socketConnection.connect(IP, PORT)) {
            hController = new HController(socketConnection);
        }
    }
    
    @AfterClass
    public static void tearDownClass() {
        socketConnection.disconnect();
    }
    
    @Test
    public void testMethods() {
        testGetLog();
        
        testClearLog();
    } 
    
    public void testGetLog() {              
        String logs = hController.getLog();
        
        assertTrue(logs.contains("qwerty"));
    }
    
    public void testClearLog() {
        hController.clearLog();
        
        String logs = hController.getLog();
        
        assertTrue(logs.equals(""));
    }
}
