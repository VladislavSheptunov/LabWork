package workshop.server.common;

import java.sql.Connection;

public interface IDataAccess {
    /**
     * Conects to the database.
     * @return true if connection was successfully completed.
     */
    public boolean connect ();

    /**
     * Disconects from database.
     */
    public void disconnect ();
    
    /**
     * Check if connection is closed.
     * @return true if connection is closed.
     */
    public boolean isConnectionClosed ();
    
    /**
     * Get connection.
     * @return connection.
     */
    public Connection getConnection ();

    /**
     * Adds user with current username and password.
     * @param username is name of user.
     * @param password is password of user.
     * @return true if user was successfully added.
     */
    public boolean addUser(final String username, final String password);
    
    /**
     * Checks is user exists in database.
     * @param username is name of user.
     * @param password is password of user.
     * @return true if user exists.
     */
    public boolean checkUser(final String username, final String password);
    
    /**
     * Gets shedule from database.
     * @return string of busy days, separated by ";".
     */
    public String getShedule();
    
    /**
     * Makes order.
     * @param username is name of user.
     * @param description is description of the proplem.
     * @param phone is phine number.
     * @param time is time from shedule.
     * @return true if order was made successfully. 
     */
    public boolean makeOrder(final String username, final String description, 
            final String phone, final String time);
			
    /**
     * Gets orders of current user.
     * @param username is user's name.
     * @return string of orders, separated by ";".
     */
    public String getOrders(final String username);
	
    /**  
     * Sets status to order as accepted.  
     * @param username is name of user.  
     * @param time is time from shedule.  
     * @return true if operation was made successfully.   
     */  
    public boolean acceptOrder(final String username, final String time);  
      
    /**  
     * Sets status to order as rejected.  
     * @param username is name of user.  
     * @param time is time from shedule.  
     * @return true if operation was made successfully.   
     */  
    public boolean rejectOrder(final String username, final String time);
    
    public String getInfo(final String time);

    public String setStatus(String time, String status, String statusDescription);

    public String deleteRecord(String time);

    public String changeTime(String oldTime, String newTime);
}
