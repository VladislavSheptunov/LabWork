package workshop.client.user;

public interface IUController {
    /**
     * Registration of user with current login and password.
     * @param username is user's username.
     * @param password is user's password.
     * @return true if registration was successfully completed.
     * Returns false if user with current login exists.
     */
    public boolean registration(final String username, final String password);
    
    /**
     * Authorization of user with current login and password.
     * @param username is user's username.
     * @param password is user's password.
     * @return true if authorization was successfully completed.
     */
    public boolean authorization(final String username, final String password);
    
    /**
     * Gets shedule.
     * @return array of times when workshop is busy.
     */
    public String[] getShedule();
    
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
     * @return array of orders.
     */
    public String[] getOrders(final String username);
    
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
}
