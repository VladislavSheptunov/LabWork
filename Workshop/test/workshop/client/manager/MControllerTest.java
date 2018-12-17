package workshop.client.manager;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import org.junit.AfterClass;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import workshop.common.ISocketConnection;
import workshop.common.SocketConnection;

public class MControllerTest {
    private static final String IP = "localhost";
    private static final int PORT = 8888;
    private static final ISocketConnection socketConnection = new SocketConnection();
    private static MController mController = null;

    
    public MControllerTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        if (socketConnection.connect(IP, PORT)) {
            mController = new MController(socketConnection);
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
    public void testGetSchedule() {
        String[] shedule = mController.getShedule();
        
        assertNotNull(shedule); 
    }
    
    @Test
    public void testGetInfo() {
        String info = mController.getInfo("11:00-15.05.2016");
        
        assertNotNull(info.contains("абв")); 
    }
    
    @Test
    public void testSetStatus() {
        String response = mController.editStatus("11:00-15.05.2016", "В работе", "qwerty");
        
        assertTrue(response.contains("Статус был изменен")); 
    }
    
    @Test
    public void testChangeTime() {
        String response = mController.changeTime("11:00-15.05.2016", "14:00-15.05.2016");
        
        assertTrue(response.contains("Время заказа изменено"));
    }
    
    @Test
    public void testDeleteTime() {
        String response = mController.deleteRecord("12:00-15.05.2016");
        
        assertTrue(response.contains("Заказ был удален"));
    }
}