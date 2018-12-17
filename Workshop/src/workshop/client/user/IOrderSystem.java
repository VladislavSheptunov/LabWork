package workshop.client.user;

public interface IOrderSystem {
    /**
     * Gets shedule.
     * @return list of times when workshop is busy.
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
