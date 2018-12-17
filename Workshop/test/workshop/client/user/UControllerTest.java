package workshop.client.user;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.AfterClass;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import workshop.common.ISocketConnection;
import workshop.common.SocketConnection;

public class UControllerTest {
    private static final String IP = "localhost";
    private static final int PORT = 8888;
    private static final ISocketConnection socketConnection = new SocketConnection();
    private static UController uController = null;

    
    public UControllerTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        if (socketConnection.connect(IP, PORT)) {
            uController = new UController(socketConnection);
        }
    }
    
    @AfterClass
    public static void tearDownClass() {
        socketConnection.disconnect();
    }
    
    @Before
    public void setUp() {
        Connection mConnection = null;
        try {
            Class.forName("com.mysql.fabric.jdbc.FabricMySQLDriver");
            mConnection = DriverManager.getConnection(
                    "jdbc:mysql://localhost:3306/workshop", "root", "1234");
        } catch (ClassNotFoundException | SQLException ex) {
        }
        try {
            String query = "DELETE FROM users;";
            Statement statement = mConnection.createStatement();
            statement.executeUpdate(query);

            query = "DELETE FROM orders;";
            statement = mConnection.createStatement();
            statement.executeUpdate(query);

            query = "INSERT INTO users VALUES (1, 'Andrey', '12345678');";
            statement = mConnection.createStatement();
            statement.executeUpdate(query);

            query = "INSERT INTO users VALUES (2, 'Alex', 'asdfghjk');";
            statement = mConnection.createStatement();
            statement.executeUpdate(query);

            query = "INSERT INTO orders VALUES (1, 1, '11:00-15.05.2016', 'abc', 'Ожидание автомобиля', 'абв', '123456');";
            statement = mConnection.createStatement();
            statement.executeUpdate(query);

            query = "INSERT INTO orders VALUES (2, 1, '12:00-15.05.2016', 'def', 'Нужны детали', 'где', '123456');";
            statement = mConnection.createStatement();
            statement.executeUpdate(query);
        } catch (SQLException ex) {
        }
    }
    
    @Test
    public void testSuccessRegistration() {      
        boolean johnReg = uController.registration("John", "asdfghjk");
                
        assertTrue(johnReg);
    }
    
    @Test
    public void testFailedRegistration() {      
        boolean alexReg = uController.registration("Alex", "asdfghjk");
        
        assertFalse(alexReg);
    }
    
    @Test
    public void testAuthorization() {
        boolean alexReg = uController.authorization("Alex", "asdfghjk");
        
        assertTrue(alexReg);
    }
    
    @Test
    public void testGetSchedule() {
        String[] shedule = uController.getShedule();
        
        assertNotNull(shedule); 
    }
    
    @Test
    public void testSuccessMakeOrder() {
        boolean ordersWasMade = uController.makeOrder("Alex", "abc", "111111", "13:00-15.05.2016");
        
        assertTrue(ordersWasMade); 
    }
    
    @Test
    public void testFailedMakeOrder() {
        boolean ordersWasMade = uController.makeOrder("Alex", "abc", "111111", "12:00-15.05.2016");
        
        assertFalse(ordersWasMade); 
    }
    
    @Test
    public void testGetOrders() {
        String[] orders = uController.getOrders("Andrey");
        
        assertTrue(orders[0].contains("абв") && orders[1].contains("где")); 
    }
    
    @Test
    public void testAcceptOrder() {
        uController.acceptOrder("Andrey", "12:00-15.05.2016");
        
        String[] orders = uController.getOrders("Andrey");
        
        assertTrue(orders[1].contains("Согласие клиента")); 
    }
    
    @Test
    public void testRejectOrder() {
        uController.rejectOrder("Andrey", "12:00-15.05.2016");
        
        String[] orders = uController.getOrders("Andrey");
        
        assertTrue(orders[1].contains("Отказ клиента")); 
    }
}
