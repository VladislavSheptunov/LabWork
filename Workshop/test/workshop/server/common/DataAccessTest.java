package workshop.server.common;

import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Statement;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class DataAccessTest {
    
    private IDataAccess instance;
    
    public DataAccessTest() {
    }
    
    @Before
    public void setUp() {
        instance = DataAccess.getInstance();
        instance.connect();
        Connection mConnection = instance.getConnection();
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
    
    @After
    public void tearDown() {
        instance.disconnect();
    }
    
    /**
     * Test of connect method, of class DataAccess.
     */
    @Test
    public void testConnect() {
        boolean result = instance.connect();
        
        assertTrue(result);
    }

    /**
     * Test of disconnect method, of class DataAccess.
     */
    @Test
    public void testDisconnect() {
        instance.disconnect();
        boolean result = instance.isConnectionClosed();
        
        assertTrue(result);
    }

    /**
     * Test of addUser method with right user, of class DataAccess.
     */
    @Test
    public void testAddRightUser() {
        String username = "John";
        String password = "asdfghjk";
        
        boolean result = instance.addUser(username, password);
    
        assertTrue(result);
    }
    
        /**
     * Test of addUser method with wrong user, of class DataAccess.
     */
    @Test
    public void testAddWrongUser() {
        String username = "Alex";
        String password = "asdfghjk";
        
        boolean result = instance.addUser(username, password);
    
        assertFalse(result);
    }


    /**
     * Test of checkUser method with incorrect user, of class DataAccess.
     */
    @Test
    public void testCheckUserWithIncorrectUser() {
        String username = "John";
        String password = "asdfghjk";
        
        boolean result = instance.checkUser(username, password);
    
        assertFalse(result);
    }

    /**
     * Test of checkUser method with correct user, of class DataAccess.
     */
    @Test
    public void testCheckUserWithCorrectUser() {
        String username = "Alex";
        String password = "asdfghjk";
        
        boolean result = instance.checkUser(username, password);
    
        assertTrue(result);
    }

    /**
     * Test of getShedule method, of class DataAccess.
     */
    @Test
    public void testGetShedule() {
        String result = instance.getShedule();
        
        assertNotNull(result);
    }

    /**
     * Test of makeOrder method with correct order, of class DataAccess.
     */
    @Test
    public void testMakeOrderWithCorrectOrder() {
        String username = "Alex";
        String description = "abc";
        String phone = "111111";
        String time = "13:00-15.05.2016";

        boolean result = instance.makeOrder(username, description, phone, time);

        assertTrue(result);
    }

    /**
     * Test of makeOrder method with incorrect order, of class DataAccess.
     */
    @Test
    public void testMakeOrderWithIncorrectOrder() {
        String username = "Alex";
        String description = "abc";
        String phone = "111111";
        String time = "12:00-15.05.2016";

        boolean result = instance.makeOrder(username, description, phone, time);

        assertFalse(result);
    }

    /**
     * Test of getOrders method, of class DataAccess.
     */
    @Test
    public void testGetOrders() {
        String username = "Andrey";
        
        String result = instance.getOrders(username);
        
        assertNotNull(result);
    }

    /**
     * Test of acceptOrder method, of class DataAccess.
     */
    @Test
    public void testAcceptOrder() {
        String username = "Andrey";
        String time = "12:00-15.05.2016";
        
        boolean result = instance.acceptOrder(username, time);
    
        assertTrue(result);
    }

    /**
     * Test of rejectOrder method, of class DataAccess.
     */
    @Test
    public void testRejectOrder() {
        String username = "Andrey";
        String time = "12:00-15.05.2016";
        
        boolean result = instance.rejectOrder(username, time);
    
        assertTrue(result);
    }

    /**
     * Test of getInfo method, of class DataAccess.
     */
    @Test
    public void testGetInfo() {
        String time = "12:00-15.05.2016";
        
        String result = instance.getInfo(time);
        
        assertNotNull(result);
    }

    /**
     * Test of setStatus method, of class DataAccess.
     */
    @Test
    public void testSetStatus() {
        String time = "11:00-15.05.2016";
        String status = "В работе";
        String statusDescription = "qwerty";
        
        String result = instance.setStatus(time, status, statusDescription);
        
        assertEquals("Статус был изменен", result);
    }

    /**
     * Test of deleteRecord method, of class DataAccess.
     */
    @Test
    public void testDeleteRecord() {
        String time = "11:00-15.05.2016";
        
        String result = instance.deleteRecord(time);
        
        assertEquals("Заказ был удален", result);
    }

    /**
     * Test of changeTime method, of class DataAccess.
     */
    @Test
    public void testChangeTime() {
        String oldTime = "11:00-15.05.2016";
        String newTime = "13:00-15.05.2016";
        
        String result = instance.changeTime(oldTime, newTime);
        
        assertEquals("Время заказа изменено", result);
    }
    
}
