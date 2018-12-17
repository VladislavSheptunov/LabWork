package workshop.common;

import org.junit.AfterClass;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;
import workshop.client.user.UController;

public class SocketConnectionTest {
    private static final String IP = "localhost";
    private static final int PORT = 8888;
    private static final ISocketConnection socketConnection = new SocketConnection();

    public SocketConnectionTest() {
    }
    
    @Test
    public void testMethods() {
        testConnect();
        
        testSend();
        
        testDisconnect();
    }
    
    public void testConnect() {
        boolean result = socketConnection.connect(IP, PORT);
        
        assertTrue(result);
    }
    
    public void testDisconnect() {
        socketConnection.disconnect();
        
        assertFalse(socketConnection.isConnected());
    }

    public void testSend() {
        Request request = new Request();
        request.senderType = "qwerty";
        request.requestType = "qwerty";
        request.body = "qwerty";
        boolean result = socketConnection.send(request);
        
        assertTrue(result);
    }
}
